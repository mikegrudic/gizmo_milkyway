# gizmo_milkyway

Initial conditions and run setup for an isolated Milky Way in GIZMO: live stellar disc, bulge, dark halo
and magnetised gas disc, with a FIRE-3 configuration, parameter file and analysis scripts.

## Model

| Component | Mass [Msun] | Default particle mass [Msun] |
|---|---|---|
| Dark halo | 9.8e11 | 1e5 |
| Stellar discs (three thin-disc age groups, thick disc) | 2.9e10 | 1e4 |
| Bulge and stellar halo | 1.2e10 | 1e4 |
| Gas disc | 1.0e10 | 1e4 |

- **Stars and dark matter**: the action-based distribution functions of Binney & Vasiliev (2023), fitted
  to Gaia DR2, made self-consistent and sampled with AGAMA. Parameters are in `mw_bv23.ini`.
- **Gas**: an isothermal disc in vertical hydrostatic equilibrium in the total potential, with rotation
  corrected for the radial pressure gradient, following Tepper-García et al. (2024). Its surface density
  is Sigma ~ exp(-R/5 kpc - 5 kpc/R).
- **Gas metallicity**: the radial profile of stars younger than 4 Gyr from Lian et al. (2023),
  -0.049 dex/kpc with [Fe/H] = -0.07 at 8.2 kpc, constant inside 4.6 kpc and outside 14.6 kpc.
- **Magnetic field**: toroidal, with B = 10 microgauss x rho / (0.169 Msun/pc^3), after Wibking &
  Krumholz (2023).

Units are pc, Msun and km/s.

## Requirements

- Python with numpy, scipy, h5py and matplotlib.
- [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama).
- A GIZMO source tree with the FIRE-3 physics (developed on `starforge_dev` at `70d7cb7a1`).
- `patches/gizmo_ic_metallicity.patch` applied to it with `git apply`. With the patch, GIZMO keeps a
  `Metallicity` dataset found in the ICs and uses `InitMetallicity` only for particle types that have
  none. Without it, all gas gets `InitMetallicity`.

## Usage

```
# ICs, parameter file and Config.sh
python make_mw_ics.py -o runs/mw --gizmo <GIZMO source> --threads 8

# build GIZMO with that Config.sh
cluster/build.sh <GIZMO source> runs/mw/Config.sh
cp <GIZMO source>/GIZMO runs/mw/

# run on one node with 96 MPI ranks; append 1 to resume from restart files
cd runs/mw && <repo>/cluster/run.sh ./GIZMO mw.params 96
```

The default resolution is 14.9 million particles. `--gizmo` copies the cooling tables into the run
directory. The scripts in `cluster/` are written for the Genoa nodes of Flatiron's Rusty cluster: change
the `module` lines, `SYSTYPE` and `#SBATCH` directives for another machine.

The generator writes:

| File | Content |
|---|---|
| `mw_ics.hdf5` | ICs. Types: 0 gas, 1 dark matter, 2 stellar discs, 3 bulge and stellar halo. The header gives the ParticleID range of each component as `IDRange_<component>`. |
| `mw.params`, `Config.sh` | GIZMO parameter file and compile-time flags. |
| `mw_ics_rotcurve.txt`, `mw_ics_gas.txt` | Rotation curve by component; gas surface density, scale height and rotation speed. |
| `mw_ics_potential.ini` | The converged potential, readable by AGAMA. |

## Options

| Option | Default | Meaning |
|---|---|---|
| `-o` | `.` | Run directory. |
| `--setup` | `fire3` | `fire3`, or `isothermal`: gravity and isothermal hydro only, to test that the ICs stay in equilibrium. |
| `--gizmo` | none | GIZMO source directory to copy the cooling tables from. |
| `--mgas`, `--mstar`, `--mdm` | 1e4, 1e4, 1e5 | Particle masses [Msun]. Softenings scale from these as m^(1/3). |
| `--cs` | 10 | Isothermal sound speed of the gas disc [km/s]. |
| `--b0`, `--b-rho` | 10, 0.169 | Field strength [microgauss] at the reference density [Msun/pc^3]; `--b0 0` for no field. |
| `--gas-metallicity` | `lian23` | `lian23`, or `uniform` to use `InitMetallicity`. |
| `--metallicity` | 1.0 | `InitMetallicity` [solar], used for the stars in the ICs. |
| `--stellar-age` | 10 | `InitStellarAge` [Gyr]. |
| `--eps-ff` | 1.0 | `SfEffPerFreeFall`. |
| `--tmax`, `--dtsnap` | 1000, 10 | Run time and snapshot interval [Myr]. |
| `--maxmem` | 12000 | `MaxMemSize` [MB per MPI rank]. |
| `--gas-rmax` | 30 | Outer edge of the gas disc [kpc]. |
| `--threads` | 4 | Threads for the model iterations. |
| `--iterations`, `--tol`, `--seed` | 8, 0.01, 42 | Iteration limit, convergence threshold, random seed. |

## Analysis

Scripts are in `analysis/`. The Milky Way profiles in `analysis/data/` are from Kennicutt & Evans (2012).

| Script | Output |
|---|---|
| `check_equilibrium.py <run>/output` | Gas and stellar disc profiles at every snapshot and their drift from the first. |
| `fmol_series.py <run> [<run> ...]` | Molecular and atomic surface density and molecular fraction at the solar circle against time. |
| `plot_gas_profiles.py <run>` | Atomic and molecular surface density profiles at every snapshot, against the Milky Way. |
| `plot_profiles.py <snapshot>` | Atomic gas, molecular gas and star formation rate profiles, against the Milky Way. |
| `plot_sigma_sfr.py <run>` | Star formation rate surface density profile at every snapshot. |
| `plot_sfr.py <run> [<run> ...]` | Star formation histories; several runs are overlaid and labelled by the parameter they differ in. |
| `plot_ism.py <snapshot>` | Star formation rate profile and gas phase diagram. |
| `plot_phase_radial.py <run> [first] [last]` | Gas phase diagrams in radial annuli, stacked over a range of snapshots. |

## Parameter surveys

Clone a built run directory once per value of any parameter in `mw.params`, and run one per node:

```
cluster/make_grid.sh runs/mw SfEffPerFreeFall 1 0.3 0.1
sbatch --array=0-2 cluster/run_grid.sbatch runs/mw-SfEffPerFreeFall{1,0.3,0.1}
```

Submitting again resumes each run from its restart files.

## Limits

- The model is axisymmetric: a bar and spiral arms have to grow from the live disc.
- There is no hot halo gas, no initial turbulence and no gas inside about 1 kpc.
- GIZMO gives the stars in the ICs one metallicity and ages uniform between zero and twice
  `InitStellarAge`, whichever component they belong to.

## References

- Binney & Vasiliev 2023, MNRAS 520, 1832
- Vasiliev 2019, MNRAS 482, 1525 (AGAMA; the distribution-function code and `mw_bv23.ini` come from its
  Milky Way example)
- Tepper-García et al. 2024, arXiv:2406.00342
- Lian et al. 2023, arXiv:2306.14100
- Wibking & Krumholz 2023, MNRAS 526, 5065
- Kennicutt & Evans 2012, ARA&A 50, 531
