#!/usr/bin/env python
"""Sigma_SFR(R) and the gas phase diagram of one snapshot.
Usage: plot_ism.py <snapshot.hdf5> [output prefix]"""
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

MYR = 0.9778            # code time unit in Myr
RHO_CGS = 6.770e-23     # Msun/pc^3 in g/cm^3
MP = 1.6726e-24
WINDOW = 10.0           # Myr; stars younger than this count towards Sigma_SFR

snap = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(snap)[0].replace("/output/snapshot", "/snapshot")
with h5py.File(snap, "r") as f:
    time = f["Header"].attrs["Time"]
    x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
    star_R = np.hypot(*(f["PartType4/Coordinates"][:] - x0)[:, :2].T) / 1e3
    star_age = (time - f["PartType4/StellarFormationTime"][:]) * MYR
    star_m = f["PartType4/Masses"][:]
    g = f["PartType0"]
    m, T, Z = g["Masses"][:], g["Temperature"][:], g["Metallicity"][:]
    nH = g["Density"][:] * RHO_CGS * (1 - Z[:, 0] - Z[:, 1]) / MP
label = f"t = {time * MYR:.0f} Myr"


# Sigma_SFR from stars formed in the last WINDOW Myr
edges = np.arange(0, 15.5, 1.0)
young = star_age < WINDOW
sig = np.histogram(star_R[young], edges, weights=star_m[young])[0] / (np.pi * np.diff(edges**2)) / (WINDOW * 1e6)
fig, ax = plt.subplots(figsize=(5, 3.2))
ax.stairs(np.where(sig > 0, sig, np.nan), edges, color="#2a6fb0", lw=2)
ax.set_yscale("log")
ax.set_xlim(0, edges[-1])
ax.set_xlabel("R [kpc]")
ax.set_ylabel(r"$\Sigma_{\rm SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]")
ax.set_title(f"{label}, stars younger than {WINDOW:.0f} Myr", fontsize=10)
fig.savefig(prefix + "_sigma_sfr.png", dpi=200, bbox_inches="tight")
print("plot:", prefix + "_sigma_sfr.png")
print(f"  SFR over the last {WINDOW:.0f} Myr: {star_m[young].sum() / (WINDOW * 1e6):.2f} Msun/yr from {young.sum()} particles")
print("  R[kpc]:", edges[:-1].astype(int).tolist())
print("  Sigma_SFR:", np.array2string(sig, precision=4))

# mass-weighted phase diagram
xe, ye = np.linspace(-6, 3, 181), np.linspace(0.5, 8, 151)
H = np.histogram2d(np.log10(nH), np.log10(T), bins=(xe, ye), weights=m)[0] / (np.diff(xe)[0] * np.diff(ye)[0])
fig, ax = plt.subplots(figsize=(5, 3.8))
im = ax.pcolormesh(xe, ye, np.where(H > 0, H, np.nan).T, norm=LogNorm(H.max() * 1e-5, H.max()), cmap="Blues", rasterized=True)
ax.set_xlabel(r"$\log_{10}\, n_{\rm H}$ [cm$^{-3}$]")
ax.set_ylabel(r"$\log_{10}\, T$ [K]")
ax.set_title(label, fontsize=10)
fig.colorbar(im, ax=ax, label=r"d$M$/d$\log n_{\rm H}$ d$\log T$ [M$_\odot$ dex$^{-2}$]")
fig.savefig(prefix + "_phase.png", dpi=200, bbox_inches="tight")
print("plot:", prefix + "_phase.png")
lt = np.log10(T)
for name, s in (("T < 100 K", lt < 2), ("100 - 8000 K", (lt >= 2) & (lt < 3.9)), ("8000 - 1e5 K", (lt >= 3.9) & (lt < 5)), ("> 1e5 K", lt >= 5)):
    print(f"  {name:13s} {m[s].sum() / m.sum():.3f} of gas mass")
print(f"  n_H: median {np.median(nH):.3g}, max {nH.max():.3g} cm^-3; T: min {T.min():.3g}, max {T.max():.3g} K")
