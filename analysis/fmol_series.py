#!/usr/bin/env python
"""Time series of the molecular gas at the solar circle (7.5 < R < 9 kpc) for one or more runs:
Sigma_mol and Sigma_atomic (with helium, factor 1.36) and the molecular fraction of the neutral hydrogen.
Usage: fmol_series.py <run dir> [<run dir> ...]"""
import glob
import sys

import h5py
import numpy as np

MYR, RLO, RHI = 0.9778, 7.5, 9.0
AREA = np.pi * (RHI**2 - RLO**2) * 1e6  # pc^2

for run in sys.argv[1:]:
    print(run.rstrip("/"))
    print("   t[Myr]  Sigma_mol  Sigma_atomic  f_mol[%]   M_mol,total[Msun]")
    for snap in sorted(glob.glob(run + "/output/snapshot_*.hdf5")):
        try:
            with h5py.File(snap, "r") as f:
                g = f["PartType0"]
                if "MolecularMassFraction" not in g:
                    continue
                x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
                R = np.hypot(*(g["Coordinates"][:] - x0)[:, :2].T) / 1e3
                Z = g["Metallicity"][:]
                mH0 = g["Masses"][:] * (1 - Z[:, 0] - Z[:, 1]) * g["NeutralHydrogenAbundance"][:]
                mol = mH0 * g["MolecularMassFraction"][:]
                t = f["Header"].attrs["Time"] * MYR
        except OSError:  # snapshot still being written
            continue
        w = (R > RLO) & (R < RHI)
        print(f"   {t:6.0f}  {1.36 * mol[w].sum() / AREA:9.3f}  {1.36 * (mH0 - mol)[w].sum() / AREA:12.2f}"
              f"  {100 * mol[w].sum() / mH0[w].sum():8.2f}   {1.36 * mol.sum():.3g}")
