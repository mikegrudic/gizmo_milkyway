#!/usr/bin/env python
"""Sigma_SFR(R) at every snapshot of a run, from stars formed in the WINDOW Myr before each.
Usage: plot_sigma_sfr.py <run dir> [plot.png]"""
import glob
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

MYR = 0.9778   # code time unit in Myr
WINDOW = 10.0  # Myr
edges = np.arange(0, 15.5, 1.0)
mid = 0.5 * (edges[1:] + edges[:-1])

run = sys.argv[1].rstrip("/")
plot = sys.argv[2] if len(sys.argv) > 2 else run + "/sigma_sfr_history.png"
times, profiles = [], []
for snap in sorted(glob.glob(run + "/output/snapshot_*.hdf5")):
    with h5py.File(snap, "r") as f:
        if "PartType4" not in f:
            continue
        time = f["Header"].attrs["Time"]
        x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
        R = np.hypot(*(f["PartType4/Coordinates"][:] - x0)[:, :2].T) / 1e3
        young = (time - f["PartType4/StellarFormationTime"][:]) * MYR < WINDOW
        m = f["PartType4/Masses"][:]
    sig = np.histogram(R[young], edges, weights=m[young])[0] / (np.pi * np.diff(edges**2)) / (WINDOW * 1e6)
    times.append(time * MYR)
    profiles.append(sig)
    print(f"t = {time * MYR:5.0f} Myr: SFR {m[young].sum() / (WINDOW * 1e6):.2f} Msun/yr, N = {young.sum():5d};",
          np.array2string(sig[:13] * 1e3, precision=1, floatmode="fixed"), "x 1e-3")

cmap = LinearSegmentedColormap.from_list("time", plt.get_cmap("Blues")(np.linspace(0.3, 1, 256)))
norm = Normalize(times[0], times[-1])
fig, ax = plt.subplots(figsize=(5.5, 3.4))
for t, sig in zip(times, profiles):
    ax.plot(mid, np.where(sig > 0, sig, np.nan), color=cmap(norm(t)), lw=1.5)
# observed Milky Way profile (Kennicutt & Evans 2012, Fig. 7), without their CMZ point
obs = np.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)), "data/kennicutt_evans_2012_mw_sigma_sfr.txt"))[1:]
y = 10.0 ** (obs[:, 1] - 3)  # Msun/yr/kpc^2
ax.errorbar(obs[:, 0], y, yerr=(y - y * 10 ** -obs[:, 2], y * 10 ** obs[:, 2] - y), fmt="o", ms=5, mfc="white",
            color="k", lw=1, capsize=2, zorder=5, label="Milky Way (Kennicutt \\& Evans 2012)")
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.set_yscale("log")
ax.set_xlim(0, edges[-1] + 0.5)
ax.set_xlabel("R [kpc]")
ax.set_ylabel(r"$\Sigma_{\rm SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]")
fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=ax, label="t [Myr]")
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
