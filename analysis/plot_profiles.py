#!/usr/bin/env python
"""Radial profiles of atomic gas, molecular gas and star formation rate surface density in one
snapshot, on one axis with the Milky Way data, after Kennicutt & Evans (2012) Fig. 7.
Gas includes helium as a factor 1.36 on the hydrogen mass; Sigma_SFR uses stars younger than WINDOW.
Usage: plot_profiles.py <snapshot.hdf5> [plot.png]"""
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

MYR = 0.9778   # code time unit in Myr
WINDOW = 10.0  # Myr
HERE = os.path.dirname(os.path.realpath(__file__))
edges = np.arange(0, 16.5, 1.0)
mid = 0.5 * (edges[1:] + edges[:-1])
area = np.pi * np.diff(edges**2) * 1e6  # pc^2

snap = sys.argv[1]
plot = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(snap)[0].replace("/output/snapshot", "/snapshot") + "_profiles.png"
with h5py.File(snap, "r") as f:
    time = f["Header"].attrs["Time"]
    x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
    g = f["PartType0"]
    R = np.hypot(*(g["Coordinates"][:] - x0)[:, :2].T) / 1e3
    Z = g["Metallicity"][:]
    mH0 = g["Masses"][:] * (1 - Z[:, 0] - Z[:, 1]) * g["NeutralHydrogenAbundance"][:]  # neutral hydrogen
    fmol = g["MolecularMassFraction"][:]  # of the neutral hydrogen
    Rs = np.hypot(*(f["PartType4/Coordinates"][:] - x0)[:, :2].T) / 1e3
    young = (time - f["PartType4/StellarFormationTime"][:]) * MYR < WINDOW
    ms = f["PartType4/Masses"][:]
sim = {"atomic": 1.36 * np.histogram(R, edges, weights=mH0 * (1 - fmol))[0] / area,
       "mol": 1.36 * np.histogram(R, edges, weights=mH0 * fmol)[0] / area,
       "sfr": np.histogram(Rs[young], edges, weights=ms[young])[0] / area / (WINDOW * 1e-3)}  # Msun/pc^2/Gyr

series = [("atomic", "#2a78d6", "s", "atomic gas"), ("mol", "#eb6834", "o", "molecular gas"), ("sfr", "#1baf7a", "*", "SFR")]
fig, ax = plt.subplots(figsize=(6, 4.2))
for key, color, marker, _ in series:
    ax.plot(mid, np.where(sim[key] > 0, sim[key], np.nan), color=color, lw=2)
    obs = np.loadtxt(os.path.join(HERE, f"data/kennicutt_evans_2012_mw_sigma_{key}.txt"))
    y = 10.0 ** obs[:, 1]
    lo, hi = (obs[:, 2], obs[:, 3]) if obs.shape[1] > 3 else (obs[:, 2], obs[:, 2])
    ax.errorbar(obs[:, 0], y, yerr=(y - y * 10**-lo, y * 10**hi - y), fmt=marker, ms=9 if marker == "*" else 5.5,
                mfc="white", color=color, lw=1, capsize=2)
ax.set_yscale("log")
ax.set_xlim(0, edges[-1])
ax.set_ylim(1e-3, 1e3)
ax.set_xlabel(r"$R_{\rm gal}$ [kpc]")
ax.set_ylabel(r"$\Sigma_{\rm gas}$ [M$_\odot$ pc$^{-2}$] or $\Sigma_{\rm SFR}$ [M$_\odot$ pc$^{-2}$ Gyr$^{-1}$]")
handles = [Line2D([], [], color=c, marker=m, mfc="white", ms=8 if m == "*" else 5.5, lw=2, label=lab) for _, c, m, lab in series]
handles += [Line2D([], [], color="0.4", lw=2, label=f"simulation, t = {time * MYR:.0f} Myr"),
            Line2D([], [], color="0.4", marker="o", mfc="white", ls="", label="Milky Way (Kennicutt \\& Evans 2012)")]
ax.legend(handles=handles, frameon=False, fontsize=8, ncol=2, loc="upper right", columnspacing=1.5)
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
print(f"t = {time * MYR:.0f} Myr; SFR over the last {WINDOW:.0f} Myr: {ms[young].sum() / (WINDOW * 1e6):.2f} Msun/yr ({young.sum()} particles)")
print("R [kpc]:", edges[:-1].astype(int).tolist())
for key in sim:
    print(f"  {key:7s}", np.array2string(sim[key], precision=3, floatmode="maxprec", suppress_small=True))
