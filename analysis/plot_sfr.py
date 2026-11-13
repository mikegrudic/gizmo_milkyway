#!/usr/bin/env python
"""Star formation histories of one or more runs from sfr.txt.
Usage: plot_sfr.py <run dir> [<run dir> ...] [-o plot.png]"""
import glob
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt

MYR = 0.9778   # code time unit in Myr
MW_SFR = 1.9   # Msun/yr, Chomiuk & Povich (2011), the normalisation used by Kennicutt & Evans (2012)

args = sys.argv[1:]
plot = None
if "-o" in args:
    plot = args[args.index("-o") + 1]
    del args[args.index("-o"):args.index("-o") + 2]
runs = [a.rstrip("/") for a in args]
plot = plot or (runs[0] + "/sfr.png" if len(runs) == 1 else os.path.dirname(runs[0]) + "/sfr_comparison.png")
names = [os.path.basename(r) for r in runs]
title = os.path.commonprefix(names) if len(runs) > 1 else ""
labels = [n[len(title):] for n in names]
# runs that differ in a single parameter are labelled by its value
params = [dict(line.split()[:2] for line in open(r + "/mw.params") if line.strip() and line[0] != "%")
          if os.path.exists(r + "/mw.params") else {} for r in runs]
differ = [k for k in params[0] if len({p.get(k) for p in params}) > 1]
if len(differ) == 1:
    title, labels = differ[0], [p[differ[0]] for p in params]

fig, ax = plt.subplots(figsize=(5.5, 3.4))
colors = plt.get_cmap("Blues")(np.linspace(1.0, 0.35, len(runs))) if len(runs) > 1 else ["#2a6fb0"]
tmax = 0
for run, name, label, color in zip(runs, names, labels, colors):
    d = np.loadtxt(run + "/output/sfr.txt")
    t, sfr = d[:, 0] * MYR, d[:, 2]  # column 3: total SFR over gas cells [Msun/yr]
    tmax = max(tmax, t[-1])
    ax.plot(t, sfr, color=color, lw=1.8, label=label)
    recent = t > t[-1] - 20
    line = f"{name}: t = {t[-1]:.0f} Myr, SFR now {sfr[-1]:.2f}, mean over the last 20 Myr {sfr[recent].mean():.2f}, peak {sfr.max():.2f} Msun/yr at {t[sfr.argmax()]:.0f} Myr"
    # cross-check against the mass in new star particles in the latest snapshot
    snap = sorted(glob.glob(run + "/output/snapshot_*.hdf5"))[-1]
    with h5py.File(snap, "r") as f:
        if "PartType4" in f:
            tsnap = f["Header"].attrs["Time"]
            formed = f["PartType4/Masses"][:].sum()
            expect = np.trapezoid(sfr[d[:, 0] <= tsnap], t[d[:, 0] <= tsnap]) * 1e6
            line += f"; new stars at {tsnap * MYR:.0f} Myr: {formed:.3g} Msun, integrated SFR {expect:.3g} Msun"
    print(line)
ax.axhline(MW_SFR, color="0.5", lw=1, ls="--")
ax.text(0.01 * tmax, MW_SFR * 1.02, "Milky Way", color="0.35", fontsize=8, va="bottom")
ax.set_xlabel("t [Myr]")
ax.set_ylabel(r"SFR [M$_\odot$ yr$^{-1}$]")
ax.set_xlim(0, tmax)
ax.set_ylim(0)
if len(runs) > 1:
    ax.legend(frameon=False, fontsize=8, title=title.replace("_", r"\_"), title_fontsize=8, loc="upper center")
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
