#!/usr/bin/env python
"""Star formation history of a run from sfr.txt. Usage: plot_sfr.py <run dir> [plot.png]"""
import glob
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt

MYR = 0.9778  # code time unit in Myr

run = sys.argv[1].rstrip("/")
plot = sys.argv[2] if len(sys.argv) > 2 else run + "/sfr.png"
d = np.loadtxt(run + "/output/sfr.txt")
t, sfr = d[:, 0] * MYR, d[:, 2]  # column 3: total SFR over gas cells [Msun/yr]

fig, ax = plt.subplots(figsize=(5, 3.2))
ax.plot(t, sfr, color="#2a6fb0", lw=2)
ax.set_xlabel("t [Myr]")
ax.set_ylabel(r"SFR [M$_\odot$ yr$^{-1}$]")
ax.set_xlim(0, t[-1])
ax.set_ylim(0)
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
print(f"latest: t = {t[-1]:.1f} Myr, SFR = {sfr[-1]:.2f} Msun/yr, peak {sfr.max():.2f} at {t[sfr.argmax()]:.1f} Myr")

# cross-check against the mass in new star particles in the latest snapshot
snap = sorted(glob.glob(run + "/output/snapshot_*.hdf5"))[-1]
with h5py.File(snap, "r") as f:
    if "PartType4" in f:
        tsnap = f["Header"].attrs["Time"]
        formed = f["PartType4/Masses"][:].sum()
        expect = np.trapezoid(sfr[d[:, 0] <= tsnap], t[d[:, 0] <= tsnap]) * 1e6
        print(f"{snap.split('/')[-1]}: {formed:.3g} Msun in new stars; time-integrated SFR to then {expect:.3g} Msun")
