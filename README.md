# gizmo_milkyway

Initial conditions and run setup for an isolated, Milky Way-like galaxy in GIZMO: a live stellar disc,
bulge, dark halo and magnetised gas disc built to match the observed Milky Way, with the FIRE-3 physics
configuration, parameter file and analysis scripts to go with it.

One script, `make_mw_ics.py`, writes a complete run directory. The rest of the repository builds and
runs GIZMO on it and measures the result against Milky Way data.

## The model

| Component | Source | Mass [Msun] | Particles at the default resolution |
|---|---|---|---|
| Dark halo | Binney & Vasiliev (2023) | 9.8e11 | 9.8 million at 1e5 Msun |
| Stellar discs: three thin-disc age groups and a thick disc | Binney & Vasiliev (2023) | 2.9e10 | 2.9 million at 1e4 Msun |
| Bulge and stellar halo | Binney & Vasiliev (2023) | 1.2e10 | 1.2 million at 1e4 Msun |
| Gas disc | Binney & Vasiliev (2023) surface density, structure solved here | 1.0e10 | 1.0 million at 1e4 Msun |

- **Stars and dark matter** follow the action-based distribution functions that Binney & Vasiliev fitted
  to Gaia DR2 kinematics. They are iterated to self-consistency with AGAMA and then sampled. The circular
  speed is 237 km/s at 8.3 kpc. The parameters are in `mw_bv23.ini`, which is AGAMA's `SCM_MW.ini`.
- **Gas** is an isothermal disc in vertical hydrostatic equilibrium in the total potential, iterated
  together with the other components, and rotating at the circular speed corrected for the radial pressure
  gradient. This is the gas treatment of the Nexus framework (Tepper-García et al. 2024). The surface
  density is Sigma ~ exp(-R/5 kpc - 5 kpc/R), which peaks at 5 kpc and has a central hole.
- **Gas metallicity** is the radial profile of stars younger than 4 Gyr from Lian et al. (2023):
  -0.049 dex/kpc with [Fe/H] = -0.07 at 8.2 kpc, held constant inside 4.6 kpc and outside 14.6 kpc where
  they have no young-star measurement. The zero point was read off their Figure 1, because the paper
  quotes only the slope. Abundance ratios are solar, with Z_sun = 0.0142.
- **Magnetic field** is purely toroidal with strength proportional to gas density,
  B = 10 microgauss x rho / (0.169 Msun/pc^3). This is the field-density relation of Wibking & Krumholz
  (2023) applied to this disc's density, giving 6 microgauss at the density peak and 2 at the solar circle.

Units are pc, Msun and km/s throughout, so the time unit is 0.978 Myr.

## Requirements

- Python with numpy, scipy, h5py and matplotlib.
- [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama):
  `git clone https://github.com/GalacticDynamics-Oxford/Agama && cd Agama && python setup.py install`.
  On Flatiron's nix Python the build cannot find libintl; run
  `export LIBRARY_PATH=/mnt/sw/nix/store/zdnba5ggf37q6cxx1rszxabi8m260cfh-gettext-0.22.4/lib:$LIBRARY_PATH`
  first, and use `yes | python setup.py install --yes` to answer its prompts.
- A GIZMO source tree with the FIRE-3 physics, which is not public. This setup was developed on the
  `starforge_dev` branch of `pfhopkins/gizmo` at commit `70d7cb7a1`.
- `patches/gizmo_ic_metallicity.patch` applied to that tree: `git apply <path>/gizmo_ic_metallicity.patch`.
  Without it GIZMO ignores the metallicity in the ICs and gives all gas `InitMetallicity`.

### What the GIZMO patch does

On a fresh start, stock GIZMO overwrites every particle's metallicity with `InitMetallicity`. With the
patch, a `Metallicity` dataset in an HDF5 IC is kept, and `InitMetallicity` applies only to particle types
that have no such dataset. An all-zero dataset therefore stays zero.

- The header attribute `Flag_Metals` gives the number of columns. With `Flag_Metals = 1` the dataset is
  the total metal mass fraction and GIZMO fills the species at solar ratios; without the attribute GIZMO
  expects one column per species it tracks.
- A dataset with the wrong number of columns is ignored with a message in the log.
- GIZMO reads and writes metallicities for gas, newly formed stars and type-5 particles only, so the
  stars in the ICs (types 2 and 3) always get `InitMetallicity`.

## Quick start

```
# 1. ICs, parameter file and Config.sh (about 6 minutes on one node at the default resolution)
python make_mw_ics.py -o runs/mw --gizmo <GIZMO source> --threads 8

# 2. build GIZMO with that Config.sh
cluster/build.sh <GIZMO source> runs/mw/Config.sh
cp <GIZMO source>/GIZMO runs/mw/

# 3. run on one node, 96 MPI ranks
cd runs/mw && <repo>/cluster/run.sh ./GIZMO mw.params 96
```

`--gizmo` copies the cooling tables GIZMO needs (`TREECOOL`, `spcool_tables/`) into the run directory.
To resume a run that stopped, add restart flag 1: `cluster/run.sh ./GIZMO mw.params 96 1`.

The scripts in `cluster/` load modules and set `SYSTYPE` for the Genoa nodes of Flatiron's Rusty cluster.
On another machine, change the `module` lines and `SYSTYPE` in `build.sh` and `run.sh`, and the
`#SBATCH` directives in `run_grid.sbatch`.

## What the generator writes

| File | Content |
|---|---|
| `mw_ics.hdf5` | The ICs. Particle types: 0 gas, 1 dark matter, 2 stellar discs, 3 bulge and stellar halo. Gas carries `InternalEnergy`, `MagneticField` in gauss and `Metallicity`. The header records the ParticleID range of each model component as `IDRange_<component>`. |
| `mw.params` | GIZMO parameter file. |
| `Config.sh` | Compile-time flags: `FIRE_PHYSICS_DEFAULTS=3` plus MHD, Spitzer conduction and Braginskii viscosity. |
| `mw_ics_rotcurve.txt` | Circular speed against radius, total and for each potential component. |
| `mw_ics_gas.txt` | Gas surface density, mid-plane density, scale height and rotation speed against radius. |
| `mw_ics_potential.ini` | The converged potential, readable by AGAMA. |

The Config lists the MHD flags one by one and does not use `FIRE_MHD`. That flag also turns on
`MHD_B_SET_IN_PARAMS`, which replaces the field in the ICs with a uniform one.

## Generator options

| Option | Default | Meaning |
|---|---|---|
| `-o` | `.` | Run directory to write into. |
| `--setup` | `fire3` | `fire3`, or `isothermal` for the equilibrium test below. |
| `--gizmo` | none | GIZMO source directory, to copy the cooling tables from. |
| `--mgas`, `--mstar`, `--mdm` | 1e4, 1e4, 1e5 | Particle masses in Msun for gas, stars and dark matter. |
| `--cs` | 10 | Isothermal sound speed of the gas disc in km/s. It sets the disc thickness and the initial temperature. |
| `--b0`, `--b-rho` | 10, 0.169 | Field strength in microgauss at the reference gas density in Msun/pc^3. `--b0 0` gives a run without magnetic fields. |
| `--gas-metallicity` | `lian23` | `lian23` writes the radial profile to the ICs; `uniform` leaves the gas at `InitMetallicity`. |
| `--metallicity` | 1.0 | `InitMetallicity` in solar units, used for the stars in the ICs. |
| `--stellar-age` | 10 | `InitStellarAge` in Gyr; see Limits. |
| `--eps-ff` | 1.0 | `SfEffPerFreeFall`, the star formation efficiency per free-fall time. |
| `--tmax`, `--dtsnap` | 1000, 10 | Run time and snapshot interval in Myr. |
| `--maxmem` | 12000 | `MaxMemSize` in MB per MPI rank. The default suits 96 ranks on a 1.5 TB node. |
| `--gas-rmax` | 30 | Outer edge of the gas disc in kpc. |
| `--threads` | 4 | Threads for the model iterations. More than about 8 does not help, because the distribution functions are Python callbacks. |
| `--iterations`, `--tol`, `--seed` | 8, 0.01, 42 | Self-consistency iteration limit, convergence threshold and random seed. |

Softenings are set from the particle masses and scale as m^(1/3). At the default masses they are 10 pc
for the stars in the ICs, 50 pc for dark matter, 32 pc for newly formed stars and a 2.2 pc floor for the
adaptive gas softening. With these values the stellar softening is about twice the local particle
spacing inside 100 pc and a fifth of it at the solar circle.

## Cost

The default resolution is 14.9 million particles. On one 96-core AMD Genoa node the FIRE-3 MHD run
advances about 0.44 Myr per wall-clock minute over its first 20 Myr, so 150 Myr takes roughly 6 hours.
That was measured at `SfEffPerFreeFall = 0.3`. Lower values leave more dense gas and are expected to be
slower; this has not been measured yet.

## Checking the ICs

`--setup isothermal` writes a gravity plus isothermal-hydro setup with no cooling or star formation, so
any evolution is a departure from equilibrium:

```
python make_mw_ics.py -o runs/eqtest --setup isothermal --b0 0
python analysis/check_equilibrium.py runs/eqtest/output
```

In a 130 Myr test at the default resolution, without a magnetic field, the stellar disc's surface
density, thickness and velocity dispersions drifted by 1 to 2 per cent (median over 2 to 16 kpc). The gas
surface density drifted by 1 per cent and its thickness by 6 per cent, and the gas radial velocity
dispersion stayed near 1.3 km/s against a sound speed of 10 km/s.

## Analysis

All scripts are in `analysis/` and take a run directory or a snapshot; each prints its usage in its
docstring. The Milky Way comparison data in `analysis/data/` are digitised from Kennicutt & Evans (2012),
Figure 7.

| Script | Output |
|---|---|
| `check_equilibrium.py <run>/output` | Radial profiles of gas and stellar disc at every snapshot, and their drift from the first. |
| `fmol_series.py <run> [<run> ...]` | Table of molecular and atomic surface density and molecular fraction at 7.5 to 9 kpc against time. |
| `plot_gas_profiles.py <run>` | Atomic and molecular surface density profiles at every snapshot, against the Milky Way. |
| `plot_profiles.py <snapshot>` | Atomic gas, molecular gas and star formation rate profiles of one snapshot, against the Milky Way. |
| `plot_sigma_sfr.py <run>` | Star formation rate surface density profile at every snapshot. |
| `plot_sfr.py <run>` | Star formation history. |
| `plot_ism.py <snapshot>` | Star formation rate profile and gas phase diagram of one snapshot. |

## Calibrating the star formation efficiency

With the FIRE-3 default `SfEffPerFreeFall = 1`, dense gas is consumed too quickly at this resolution. A
165 Myr run at that value, with uniform solar metallicity and softer gravity than the current defaults,
settled at a molecular surface density of 0.06 Msun/pc^2 at the solar circle and a molecular fraction of
0.6 per cent. The Milky Way has about 1.4 Msun/pc^2 and 12 per cent (Kennicutt & Evans 2012).

The intended fix is to lower `SfEffPerFreeFall` until the solar-circle molecular fraction matches. A
survey over 1, 0.3, 0.1, 0.03 and 0.01 was started in November 2026; no value is recommended yet.

To run such a survey, generate and build one run, then clone it per value and submit one node each:

```
cluster/make_grid.sh runs/mw SfEffPerFreeFall 1 0.3 0.1 0.03 0.01
sbatch --array=0-4 cluster/run_grid.sbatch runs/mw-SfEffPerFreeFall{1,0.3,0.1,0.03,0.01}
python analysis/fmol_series.py runs/mw-SfEffPerFreeFall*
```

- `make_grid.sh` works for any parameter in `mw.params`. The run directory must already contain the
  binary as `GIZMO`.
- Submitting `run_grid.sbatch` again resumes each run from its restart files.
- Judge each run on its 100 to 150 Myr average. In the reference run the solar-circle molecular fraction
  spiked at 20 Myr as the disc first collapsed, fell to a trough from 40 to 80 Myr, and was flat from
  about 100 Myr.

## Limits

- The model is axisymmetric. A bar and spiral arms have to grow from the live disc.
- There is no hot halo gas and no initial turbulence. The disc starts on thermal pressure alone and
  collapses vertically once cooling is on, until feedback takes over. The magnetic pressure, about
  20 per cent of thermal at the density peak, is not part of the hydrostatic solution.
- GIZMO gives the stars in the ICs ages uniform between zero and twice `InitStellarAge`, whichever
  component they belong to. The default of 10 Gyr makes them equivalent to 2 Msun/yr of ongoing star
  formation, but spreads the young stars over the bulge and thick disc as well as the thin disc.
- The stars in the ICs all have one metallicity, `InitMetallicity`.
- The gas disc of Binney & Vasiliev has no gas inside about 1 kpc, so there is no central molecular zone.

## References

- Binney & Vasiliev 2023, MNRAS 520, 1832: the stellar and dark matter model.
- Vasiliev 2019, MNRAS 482, 1525: AGAMA. The distribution-function code in `make_mw_ics.py` and
  `mw_bv23.ini` are taken from its Milky Way example.
- Tepper-García et al. 2024, MNRAS 535, arXiv:2406.00342: Nexus, the gas disc method.
- Lian et al. 2023, Nature Astronomy, arXiv:2306.14100: the metallicity profile.
- Wibking & Krumholz 2023, MNRAS 526, 5065: the magnetic field.
- Kennicutt & Evans 2012, ARA&A 50, 531: the Milky Way gas and star formation profiles.
