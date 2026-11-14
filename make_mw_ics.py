#!/usr/bin/env python
"""Milky Way initial conditions for GIZMO.

Collisionless components: the action-based DF model of Binney & Vasiliev (2023, MNRAS 520, 1832),
fit to Gaia DR2, built with AGAMA's self-consistent modelling machinery (DFs and parameters
taken from AGAMA's example_self_consistent_model_mw.py / SCM_MW.ini).

Gas: a live isothermal disc in vertical hydrostatic equilibrium in the total potential, iterated
together with the DF components, with rotation corrected for the radial pressure gradient
(Tepper-Garcia et al. 2024, arXiv:2406.00342, "Nexus", sec. 3).

Magnetic field: purely toroidal with |B| proportional to the gas density, as in Wibking & Krumholz
(2023, MNRAS 526, 5065), whose field and gas share one double-exponential profile with 10 microgauss
at 0.169 Msun/pc^3. Their relation is applied to the density of this disc rather than their profile.

Gas metallicity: the radial profile of stars younger than 4 Gyr from Lian et al. (2023, arXiv:2306.14100),
written to the ICs. GIZMO needs patches/gizmo_ic_metallicity.patch to read it.

Writes <outdir>/mw_ics.hdf5, mw.params and Config.sh for either a FIRE-3 run (default) or an
isothermal gravity + hydro run that tests the equilibrium of the ICs (--setup isothermal).
Internally AGAMA works in kpc, km/s, Msun; the snapshot and parameter file use pc, km/s, Msun.
"""
import argparse
import os
import shutil
import time
import warnings
from configparser import RawConfigParser

import agama
import h5py
import numpy as np
from scipy.integrate import cumulative_trapezoid, fixed_quad, trapezoid
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
from scipy.stats import qmc

HERE = os.path.dirname(os.path.realpath(__file__))


# --- DFs of Binney & Vasiliev (2023), as implemented in AGAMA's example script -------------------

def double_power_law_df(**p):
    J0, L0, slopeIn, slopeOut, rotFrac, Jphi0, alpha, beta, Fin, Fout, Jcutoff, cutoffStrength, Jcore = (
        float(p[k]) for k in ("J0", "L0", "slopeIn", "slopeOut", "rotFrac", "Jphi0", "alpha", "beta",
                              "Fin", "Fout", "Jcutoff", "cutoffStrength", "Jcore"))

    def df(J):
        modJphi = abs(J[:, 2])
        L = J[:, 1] + modJphi
        c = L / (L + J[:, 0])
        jta = ((1.5 * J[:, 0] + L) / L0) ** alpha
        xi = jta / (1 + jta)
        rat = (1 - xi) * Fin + xi * Fout
        a, b = 0.5 * (rat + 1), 0.5 * (rat - 1)
        cL = np.where(L > 0, J[:, 1] * (a + b * modJphi / L) + modJphi, 0)
        fac = np.exp(beta * np.sin(np.pi / 2 * c))
        hJ = J[:, 0] / fac + 0.5 * (1 + c * xi) * fac * cL
        result = norm / (2 * np.pi * J0) ** 3 * (1 + J0 / hJ) ** slopeIn * (1 + hJ / J0) ** (-slopeOut)
        if Jcutoff > 0:
            result *= np.exp(-((hJ / Jcutoff) ** cutoffStrength))
        if Jcore > 0:
            result *= (1 + Jcore / hJ * (Jcore / hJ - zeta)) ** (-0.5 * slopeIn)
        if rotFrac != 0:
            result *= 1 + rotFrac * np.tanh(J[:, 2] / Jphi0)
        return result

    if Jcore > 0:  # Cole & Binney core: zeta chosen so that the core conserves the total mass
        def root(zeta):
            def integrand(t):
                hJ = Jcore * t * t * (3 - 2 * t) / (1 - t) ** 2 / (1 + 2 * t)
                dhJdt = Jcore * 6 * t / (1 - t) ** 3 / (1 + 2 * t) ** 2
                return (hJ**2 * dhJdt * (1 + J0 / hJ) ** slopeIn * (1 + hJ / J0) ** (-slopeOut)
                        * ((1 + Jcore / hJ * (Jcore / hJ - zeta)) ** (-0.5 * slopeIn) - 1))
            return fixed_quad(integrand, 0.0, 1.0, n=20)[0]
        zeta = brentq(root, 0.0, 2.0)
    norm = 1.0
    norm = float(p["mass"]) / agama.DistributionFunction(df).totalMass()
    return df


def exponential_df(**p):
    Jr0, Jz0, Jphi0, pr, pz, addJden, addJvel = (
        float(p[k]) for k in ("Jr0", "Jz0", "Jphi0", "pr", "pz", "addJden", "addJvel"))

    def df(J):
        Jp = np.maximum(0, J[:, 2])
        Jvel, Jden = Jp + addJvel, Jp + addJden
        xr = (Jvel / Jphi0) ** pr / Jr0
        xz = (Jvel / Jphi0) ** pz / Jz0
        fr = xr * np.exp(-xr * J[:, 0])
        fz = xz * np.exp(-xz * J[:, 1])
        fp = norm / Jphi0**2 * abs(J[:, 2]) * np.exp(-Jden / Jphi0)
        return np.where(J[:, 2] > 0, fr * fz * fp, 0)

    norm = 1.0
    norm = float(p["mass"]) / agama.DistributionFunction(df).totalMass()
    return df


# (ini section, DF constructor, GIZMO particle type)
STELLAR_COMPONENTS = [
    ("bulge", "DF bulge", double_power_law_df, 3),
    ("disk_young", "DF young disk", exponential_df, 2),
    ("disk_middle", "DF middle disk", exponential_df, 2),
    ("disk_old", "DF old disk", exponential_df, 2),
    ("disk_thick", "DF highA disk", exponential_df, 2),
    ("stellar_halo", "DF stellar halo", double_power_law_df, 3),
]


# --- gas disc ------------------------------------------------------------------------------------

class GasDisc:
    """Isothermal (P = cs^2 rho) disc of prescribed surface density Sigma(R) in vertical hydrostatic
    equilibrium: rho(R,z) = rho0(R) exp(-[Phi(R,z) - Phi(R,0)] / cs^2)."""

    def __init__(self, surfaceDensity, scaleRadius, innerCutoffRadius, cs, Rmax, zmax=10.0):
        self.Sigma0, self.Rd, self.Rcut, self.cs = surfaceDensity, scaleRadius, innerCutoffRadius, cs
        self.R = np.geomspace(0.05, Rmax, 300)
        zs = 0.01  # kpc; z grid is linear below this and logarithmic above
        self.z = zs * np.sinh(np.linspace(0, np.arcsinh(zmax / zs), 300))
        self.Sigma = self.Sigma0 * np.exp(-self.R / self.Rd - self.Rcut / self.R)
        self.rho0 = self.profile = None

    def solve(self, pot):
        """Recompute the density in the potential `pot`; returns max |dln rho0| within the
        radii holding the central 99% of the mass (the previous solution is the reference)."""
        RR, zz = np.meshgrid(self.R, self.z, indexing="ij")
        Phi = pot.potential(np.column_stack((RR.ravel(), 0 * RR.ravel(), zz.ravel()))).reshape(RR.shape)
        lnprof = -(Phi - Phi[:, :1]) / self.cs**2
        self.profile = np.exp(lnprof)
        rho0 = self.Sigma / (2 * trapezoid(self.profile, self.z, axis=1))
        change = np.inf
        if self.rho0 is not None:
            m = cumulative_trapezoid(self.R * self.Sigma, self.R, initial=0)
            sel = (m > 0.005 * m[-1]) & (m < 0.995 * m[-1])
            change = np.max(np.abs(np.log(rho0 / self.rho0))[sel])
        self.rho0 = rho0
        lnrho = np.log(rho0)[:, None] + lnprof
        self._lnrho = RegularGridInterpolator((self.R, self.z), lnrho, bounds_error=False, fill_value=-np.inf)
        return change

    def density(self, xyz):
        R = np.maximum(np.hypot(xyz[:, 0], xyz[:, 1]), self.R[0])
        return np.exp(self._lnrho(np.column_stack((R, np.abs(xyz[:, 2])))))

    def agama_density(self):
        return agama.Density(self.density, symmetry="a")

    def scale_height(self):
        """Half-thickness Sigma / (2 rho0), the sech^2 scale height for a self-gravitating slab."""
        return self.Sigma / (2 * self.rho0)

    def vphi(self, pot):
        """Rotation speed on the R grid: v_phi^2 = v_c^2 + cs^2 dln(rho0)/dln(R).
        The disc is barotropic, so v_phi is independent of z."""
        xyz = np.column_stack((self.R, 0 * self.R, 0 * self.R))
        vc2 = -self.R * pot.force(xyz)[:, 0]
        v2 = vc2 + self.cs**2 * np.gradient(np.log(self.rho0), np.log(self.R))
        return np.sqrt(np.maximum(v2, 0)), np.sqrt(vc2)

    def sample(self, pot, mass, seed):
        """Equal-mass cells from a scrambled Halton sequence (lower shot noise than random draws)."""
        Mcum = cumulative_trapezoid(2 * np.pi * self.R * self.Sigma, self.R, initial=0)
        n = int(round(Mcum[-1] / mass))
        u = qmc.Halton(3, scramble=True, seed=seed).random(n)
        R = np.interp(u[:, 0] * Mcum[-1], Mcum, self.R)
        phi = 2 * np.pi * u[:, 1]
        # inverse vertical CDF, tabulated on (R, quantile)
        q = np.concatenate((np.linspace(0, 0.99, 200)[:-1], 1 - np.geomspace(1e-2, 1e-8, 80)))
        cdf = cumulative_trapezoid(self.profile, self.z, axis=1, initial=0)
        cdf /= cdf[:, -1:]
        zq = np.array([np.interp(q, c, self.z) for c in cdf])
        s = 2 * u[:, 2] - 1
        z = np.sign(s) * RegularGridInterpolator((self.R, q), zq)(
            np.column_stack((R, np.minimum(np.abs(s), q[-1]))))
        v = np.interp(R, self.R, self.vphi(pot)[0])
        pos = np.column_stack((R * np.cos(phi), R * np.sin(phi), z))
        vel = np.column_stack((-v * np.sin(phi), v * np.cos(phi), 0 * v))
        return pos, vel, np.full(n, Mcum[-1] / n)


# --- model construction --------------------------------------------------------------------------

def build_model(ini_file, cs, gas_rmax, n_iter, tol, threads):
    ini = RawConfigParser()
    ini.optionxform = str
    ini.read(ini_file)
    sec = lambda name: dict(ini.items(name))
    agama.setUnits(length=1, velocity=1, mass=1)  # kpc, km/s, Msun

    model = agama.SelfConsistentModel(**sec("SelfConsistentModel"))
    model.potential = agama.Potential(ini_file)  # initial guess

    t0 = time.time()
    dfs = {name: make(**sec(section)) for name, section, make, _ in STELLAR_COMPONENTS}
    dfs["dark_halo"] = double_power_law_df(**sec("DF dark halo"))
    print(f"DFs initialised in {time.time() - t0:.0f} s", flush=True)

    g = sec("Potential gas disk")
    gas = GasDisc(float(g["surfaceDensity"]), float(g["scaleRadius"]), float(g["innerCutoffRadius"]),
                  cs, gas_rmax)
    df_stars = agama.DistributionFunction(*[dfs[c[0]] for c in STELLAR_COMPONENTS])
    comp_stars = agama.Component(df=df_stars, disklike=True, **sec("SelfConsistentModel disk"))
    comp_halo = agama.Component(df=dfs["dark_halo"], disklike=False, **sec("SelfConsistentModel halo"))

    # Python DFs hold the GIL, so more than a few threads do not help the density integrals
    with agama.setNumThreads(threads):
        for it in range(1, n_iter + 1):
            t0 = time.time()
            phi0 = model.potential.potential(0, 0, 0)
            dgas = gas.solve(model.potential)
            model.components = [comp_stars, comp_halo, agama.Component(density=gas.agama_density(), disklike=True)]
            model.iterate()
            dphi = abs(model.potential.potential(0, 0, 0) / phi0 - 1)
            print(f"iteration {it}: dPhi(0)/Phi(0) = {dphi:.2e}, max |dln rho0_gas| = {dgas:.2e}, "
                  f"{time.time() - t0:.0f} s", flush=True)
            if it > 1 and dphi < tol and dgas < tol:
                break
    # residual: the gas disc is one potential update behind
    resid = GasDisc(gas.Sigma0, gas.Rd, gas.Rcut, cs, gas_rmax)
    resid.rho0 = gas.rho0
    print(f"gas hydrostatic residual in the final potential: max |dln rho0| = {resid.solve(model.potential):.2e}")
    return model, dfs, gas


def write_diagnostics(prefix, model, gas):
    R = np.geomspace(0.05, 50, 120)
    xyz = np.column_stack((R, 0 * R, 0 * R))
    vc = np.column_stack([np.sqrt(np.maximum(-R * p.force(xyz)[:, 0], 0)) for p in model.potential])
    np.savetxt(prefix + "_rotcurve.txt", np.column_stack((R, np.sqrt(np.sum(vc**2, axis=1)), vc)), fmt="%.6g",
               header="R[kpc] vcirc_total[km/s] vcirc of each potential component")
    vphi, vcg = gas.vphi(model.potential)
    np.savetxt(prefix + "_gas.txt", np.column_stack((gas.R, gas.Sigma * 1e-6, gas.rho0 * 1e-9,
                                                    gas.scale_height() * 1e3, vphi, vcg)), fmt="%.6g",
               header="R[kpc] Sigma[Msun/pc^2] rho0[Msun/pc^3] h=Sigma/(2rho0)[pc] vphi[km/s] vcirc[km/s]")
    model.potential.export(prefix + "_potential.ini")


ZSUN = 0.0142  # solar metal mass fraction in FIRE-3 (Asplund et al. 2009, proto-solar)


def lian23_young_feh(R):
    """[Fe/H] of stars younger than 4 Gyr at radius R [kpc], taken as the gas metallicity (Lian et al. 2023,
    arXiv:2306.14100): gradient -0.049 dex/kpc with the zero point read off their Fig. 1 (-0.07 at 8.2 kpc),
    held constant outside the measured range of 4.6 to 14.6 kpc."""
    return -0.131 - 0.049 * (np.clip(R, 4.6, 14.6) - 9.55)


def toroidal_field(gas, pos, b0, rho_ref):
    """B [Gauss] at the gas cells: B_phi = b0 [microgauss] * rho / rho_ref [Msun/pc^3]."""
    B = 1e-6 * b0 * gas.density(pos) * 1e-9 / rho_ref
    phi = np.arctan2(pos[:, 1], pos[:, 0])
    return np.column_stack((-B * np.sin(phi), B * np.cos(phi), 0 * B))


def write_gizmo(fname, parts, cs, gamma, component_ids, bfield=None, metallicity=None):
    """parts: {ptype: (pos[kpc], vel[km/s], mass[Msun], first_id)}; bfield: gas B [Gauss] or None;
    metallicity: gas metal mass fraction or None"""
    npart = np.zeros(6, dtype=np.int64)
    for t, p in parts.items():
        npart[t] = len(p[2])
    with h5py.File(fname, "w") as f:
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = npart.astype(np.uint32)
        h.attrs["NumPart_Total"] = npart.astype(np.uint32)
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6)
        h.attrs["Time"] = 0.0
        h.attrs["Redshift"] = 0.0
        h.attrs["BoxSize"] = 0.0
        h.attrs["NumFilesPerSnapshot"] = 1
        h.attrs["Flag_DoublePrecision"] = 0
        h.attrs["UnitLength_In_CGS"] = 3.085678e18
        h.attrs["UnitMass_In_CGS"] = 1.989e33
        h.attrs["UnitVelocity_In_CGS"] = 1e5
        if metallicity is not None:
            h.attrs["Flag_Metals"] = np.int32(1)  # Metallicity holds the total only; GIZMO assumes solar ratios
        for name, (first, n) in component_ids.items():  # ParticleID ranges of each model component
            h.attrs["IDRange_" + name] = np.array([first, first + n - 1], dtype=np.uint64)
        for t, (pos, vel, mass, first_id) in parts.items():
            g = f.create_group(f"PartType{t}")
            g["Coordinates"] = (pos * 1e3).astype(np.float32)  # kpc -> pc
            g["Velocities"] = vel.astype(np.float32)
            g["Masses"] = mass.astype(np.float32)
            g["ParticleIDs"] = np.arange(first_id, first_id + len(mass), dtype=np.uint32)
            if t == 0:
                g["InternalEnergy"] = np.full(len(mass), cs**2 / (gamma - 1), dtype=np.float32)
                if bfield is not None:
                    g["MagneticField"] = bfield.astype(np.float32)
                if metallicity is not None:
                    g["Metallicity"] = metallicity.astype(np.float32)


CONFIG = {
    "fire3": """FIRE_PHYSICS_DEFAULTS=3
SNE_BROADCAST_TIMESTEP                # without it, supernovae create energy in gas on long timesteps
INIT_STELLAR_METALS_AGES_DEFINED      # InitMetallicity where the ICs give none, ages of the IC stars from InitStellarAge
MULTIPLEDOMAINS=16
""",
    "isothermal": """# equilibrium test of the ICs: self-gravity + isothermal hydro, no cooling or star formation
HYDRO_MESHLESS_FINITE_MASS
EOS_GAMMA={gamma}
EOS_ENFORCE_ADIABAT={cs2}       # P = cs^2 rho^gamma, cs = {cs} km/s
ADAPTIVE_GRAVSOFT_FORGAS
MULTIPLEDOMAINS=16
""",
}

PARAMS = """% Milky Way ({setup}): m_gas = {mgas:g}, m_star = {mstar:g}, m_dm = {mdm:g} Msun
% Units: pc, Msun, km/s (time unit 0.978 Myr)
InitCondFile                       mw_ics
OutputDir                          output
TimeBegin                          0
TimeMax                            {tmax:.6g}      % {tmax_myr:g} Myr
TimeBetSnapshot                    {dtsnap:.6g}      % {dtsnap_myr:g} Myr
TimeOfFirstSnapshot                0
TimeBetStatistics                  1
MaxSizeTimestep                    1
MinSizeTimestep                    0
BoxSize                            1e6          % unused (open boundaries)
UnitLength_in_cm                   3.085678e18
UnitMass_in_g                      1.989e33
UnitVelocity_in_cm_per_s           1e5
UnitMagneticField_in_gauss         1
InitGasTemp                        0            % from ICs
MinGasTemp                         {mintemp}
{fire}DesNumNgb                          32
MaxNumNgbDeviation                 0.1
% softenings scale as m^(1/3): gas floor (adaptive above it), DM, IC stars (disc, bulge + halo), new stars
Softening_Type0                    {eps_gas:.3g}
Softening_Type1                    {eps_dm:.3g}
Softening_Type2                    {eps_star:.3g}
Softening_Type3                    {eps_star:.3g}
Softening_Type4                    {eps_new:.3g}
Softening_Type5                    {eps_new:.3g}
Softening_Type0_MaxPhysLimit       {eps_gas:.3g}
Softening_Type1_MaxPhysLimit       {eps_dm:.3g}
Softening_Type2_MaxPhysLimit       {eps_star:.3g}
Softening_Type3_MaxPhysLimit       {eps_star:.3g}
Softening_Type4_MaxPhysLimit       {eps_new:.3g}
Softening_Type5_MaxPhysLimit       {eps_new:.3g}
ErrTolIntAccuracy                  0.01
CourantFac                         0.2
MaxRMSDisplacementFac              0.1
ErrTolForceAcc                     0.001
ErrTolTheta                        0.7
TreeDomainUpdateFrequency          0.01
MaxMemSize                         {maxmem}         % MB per MPI rank
PartAllocFactor                    3.0
BufferSize                         200
TimeLimitCPU                       1e9
CpuTimeBetRestartFile              3600
ResubmitOn                         0
ResubmitCommand                    none
"""

FIRE_PARAMS = """InitMetallicity                    {Z:g}            % solar units: IC stars, and the gas if the ICs give no Metallicity
SfEffPerFreeFall                   {eps_ff:g}
CritPhysDensity                    {nth:g}            % star formation density threshold, n_H [cm^-3]
InitStellarAge                     {age:g}           % Gyr; IC stars get ages uniform in [0, 2x this]: {sfr:.2g} Msun/yr of stars forming
"""

# FIRE_MHD without MHD_B_SET_IN_PARAMS, which would replace the IC field by a uniform one
MHD_CONFIG = {
    "fire3": """MAGNETIC                              # field read from the ICs
CONDUCTION
CONDUCTION_SPITZER
VISCOSITY
VISCOSITY_BRAGINSKII
""",
    "isothermal": "MAGNETIC\n",
}

MYR = 0.9778  # code time unit in Myr


def write_setup(a, mstar_total):
    iso = a.setup == "isothermal"
    with open(os.path.join(a.outdir, "Config.sh"), "w") as f:
        f.write(CONFIG[a.setup].format(gamma=a.gamma, cs2=a.cs**2, cs=a.cs))
        if a.b0 > 0:
            f.write(MHD_CONFIG[a.setup])
    fire = "" if iso else FIRE_PARAMS.format(Z=a.metallicity, eps_ff=a.eps_ff, nth=a.sf_density, age=a.stellar_age,
                                             sfr=mstar_total / (2e9 * a.stellar_age))
    with open(os.path.join(a.outdir, "mw.params"), "w") as f:
        f.write(PARAMS.format(
            setup="isothermal equilibrium test" if iso else "FIRE-3", mgas=a.mgas, mstar=a.mstar, mdm=a.mdm,
            tmax=a.tmax / MYR, tmax_myr=a.tmax, dtsnap=a.dtsnap / MYR, dtsnap_myr=a.dtsnap,
            mintemp=0 if iso else 10, fire=fire, maxmem=a.maxmem,
            eps_gas=(a.mgas / 1e3) ** (1 / 3), eps_dm=50 * (a.mdm / 1e5) ** (1 / 3),
            eps_star=10 * (a.mstar / 1e4) ** (1 / 3), eps_new=15 * (a.mgas / 1e3) ** (1 / 3)))
    if not iso and a.gizmo:  # cooling tables
        shutil.copy(os.path.join(a.gizmo, "cooling", "TREECOOL"), a.outdir)
        if not os.path.exists(os.path.join(a.outdir, "spcool_tables")):
            tgz = os.path.join(a.outdir, "spcool_tables.tgz")
            from urllib.request import urlretrieve
            urlretrieve("https://users.flatironinstitute.org/~mgrudic/gizmo_tests/spcool_tables.tgz", tgz)
            shutil.unpack_archive(tgz, a.outdir)
            os.remove(tgz)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ini", default=os.path.join(HERE, "mw_bv23.ini"))
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--setup", choices=("fire3", "isothermal"), default="fire3")
    ap.add_argument("--gizmo", help="GIZMO source directory; if given, the FIRE cooling tables are put in outdir")
    ap.add_argument("--mgas", type=float, default=1e4, help="gas cell mass [Msun]")
    ap.add_argument("--mstar", type=float, default=1e4, help="star particle mass [Msun]")
    ap.add_argument("--mdm", type=float, default=1e5, help="dark matter particle mass [Msun]")
    ap.add_argument("--cs", type=float, default=10.0,
                    help="isothermal sound speed of the gas disc, sqrt(P/rho) [km/s]")
    ap.add_argument("--b0", type=float, default=10.0,
                    help="toroidal field strength [microgauss] at gas density --b-rho; 0 for no magnetic field")
    ap.add_argument("--b-rho", type=float, default=0.169, help="reference gas density for --b0 [Msun/pc^3]")
    ap.add_argument("--metallicity", type=float, default=1.0, help="InitMetallicity [solar], used for the IC stars")
    ap.add_argument("--gas-metallicity", choices=("lian23", "uniform"), default="lian23",
                    help="lian23: radial profile of young stars from Lian et al. (2023), written to the ICs "
                         "(needs a GIZMO that reads IC metallicities); uniform: InitMetallicity")
    ap.add_argument("--eps-ff", type=float, default=1.0, help="SfEffPerFreeFall, star formation efficiency per free-fall time")
    ap.add_argument("--sf-density", type=float, default=1.0,
                    help="CritPhysDensity, density threshold for star formation, n_H [cm^-3]; 1 is the FIRE-3 default")
    ap.add_argument("--stellar-age", type=float, default=10.0, help="InitStellarAge [Gyr]")
    ap.add_argument("--tmax", type=float, default=1000.0, help="run time [Myr]")
    ap.add_argument("--dtsnap", type=float, default=10.0, help="snapshot interval [Myr]")
    ap.add_argument("--maxmem", type=int, default=12000, help="MaxMemSize [MB per MPI rank]")
    ap.add_argument("--gas-rmax", type=float, default=30.0, help="outer edge of the gas disc [kpc]")
    ap.add_argument("--iterations", type=int, default=8, help="maximum number of self-consistency iterations")
    ap.add_argument("--tol", type=float, default=0.01, help="convergence threshold (relative change)")
    ap.add_argument("--threads", type=int, default=4, help="OpenMP threads for the model iterations")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    agama.setRandomSeed(a.seed)
    a.gamma = 1.001 if a.setup == "isothermal" else 5.0 / 3.0  # InternalEnergy = cs^2 / (gamma - 1)
    os.makedirs(a.outdir, exist_ok=True)
    a.output = os.path.join(a.outdir, "mw_ics.hdf5")
    prefix = os.path.join(a.outdir, "mw_ics")

    model, dfs, gas = build_model(a.ini, a.cs, a.gas_rmax, a.iterations, a.tol, a.threads)
    write_diagnostics(prefix, model, gas)

    t0 = time.time()
    parts, component_ids, next_id = {}, {}, 1

    def add(ptype, name, pos, vel, mass):
        nonlocal next_id
        component_ids[name] = (next_id, len(mass))
        next_id += len(mass)
        print(f"  {name:13s} type {ptype}: N = {len(mass):9d}, M = {mass.sum():.4g} Msun, "
              f"m = {mass[0]:.4g} Msun", flush=True)
        if ptype in parts:
            p = parts[ptype]
            parts[ptype] = (np.vstack((p[0], pos)), np.vstack((p[1], vel)), np.concatenate((p[2], mass)), p[3])
        else:
            parts[ptype] = (pos, vel, mass, component_ids[name][0])

    add(0, "gas", *gas.sample(model.potential, a.mgas, a.seed))
    components = [("dark_halo", 1, a.mdm)] + [(c[0], c[3], a.mstar) for c in STELLAR_COMPONENTS]
    for name, ptype, m in sorted(components, key=lambda c: c[1]):  # contiguous IDs within each type
        gm = agama.GalaxyModel(potential=model.potential, df=dfs[name], af=model.af)
        xv, mass = gm.sample(int(round(agama.DistributionFunction(dfs[name]).totalMass() / m)))
        add(ptype, name, xv[:, :3], xv[:, 3:], mass)
    print(f"sampled in {time.time() - t0:.0f} s")

    mtot = sum(p[2].sum() for p in parts.values())
    bfield = toroidal_field(gas, parts[0][0], a.b0, a.b_rho) if a.b0 > 0 else None
    if bfield is not None:
        print(f"toroidal B: max {1e6 * np.abs(bfield).max():.2f} microgauss")
    metallicity = None
    if a.setup == "fire3" and a.gas_metallicity == "lian23":
        metallicity = ZSUN * 10 ** lian23_young_feh(np.hypot(parts[0][0][:, 0], parts[0][0][:, 1]))
        print(f"gas metallicity: {metallicity.min() / ZSUN:.2f} to {metallicity.max() / ZSUN:.2f} solar, "
              f"mass-weighted mean {metallicity.mean() / ZSUN:.2f}")
    write_gizmo(a.output, parts, a.cs, a.gamma, component_ids, bfield, metallicity)
    write_setup(a, sum(parts[t][2].sum() for t in (2, 3)))
    print(f"wrote {a.output}: {next_id - 1} particles, {mtot:.4g} Msun; mw.params and Config.sh in {a.outdir}")


if __name__ == "__main__":
    main()
