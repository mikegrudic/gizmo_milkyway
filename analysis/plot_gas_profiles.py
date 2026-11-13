#!/usr/bin/env python
"""Atomic and molecular gas surface density profiles at every snapshot of a run, with the Milky Way
profiles of Kennicutt & Evans (2012). Both include helium as a factor 1.36 on the hydrogen mass.
Usage: plot_gas_profiles.py <run dir> [plot.png]"""
import glob
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

MYR = 0.9778  # code time unit in Myr
HERE = os.path.dirname(os.path.realpath(__file__))
edges = np.arange(0, 16.5, 1.0)
mid = 0.5 * (edges[1:] + edges[:-1])
area = np.pi * np.diff(edges**2) * 1e6  # pc^2

run = sys.argv[1].rstrip("/")
plot = sys.argv[2] if len(sys.argv) > 2 else run + "/gas_profiles.png"
times, atomic, mol = [], [], []
for snap in sorted(glob.glob(run + "/output/snapshot_*.hdf5")):
    with h5py.File(snap, "r") as f:
        g = f["PartType0"]
        if "MolecularMassFraction" not in g:
            continue
        x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
        R = np.hypot(*(g["Coordinates"][:] - x0)[:, :2].T) / 1e3
        Z = g["Metallicity"][:]
        mH0 = g["Masses"][:] * (1 - Z[:, 0] - Z[:, 1]) * g["NeutralHydrogenAbundance"][:]  # neutral hydrogen
        fmol = g["MolecularMassFraction"][:]  # of the neutral hydrogen
        times.append(f["Header"].attrs["Time"] * MYR)
    atomic.append(1.36 * np.histogram(R, edges, weights=mH0 * (1 - fmol))[0] / area)
    mol.append(1.36 * np.histogram(R, edges, weights=mH0 * fmol)[0] / area)
    print(f"t = {times[-1]:4.0f} Myr: M_atomic = {atomic[-1] @ area:.3g}, M_mol = {mol[-1] @ area:.3g} Msun (with He)")

cmap = LinearSegmentedColormap.from_list("time", plt.get_cmap("Blues")(np.linspace(0.3, 1, 256)))
norm = Normalize(times[0], times[-1])
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), sharex=True, gridspec_kw=dict(wspace=0.28))
for ax, prof, name in zip(axes, (atomic, mol), ("atomic", "mol")):
    for t, p in zip(times, prof):
        ax.plot(mid, np.where(p > 0, p, np.nan), color=cmap(norm(t)), lw=1.5)
    obs = np.loadtxt(os.path.join(HERE, f"data/kennicutt_evans_2012_mw_sigma_{name}.txt"))
    y = 10.0 ** obs[:, 1]
    ax.errorbar(obs[:, 0], y, yerr=(y - y * 10 ** -obs[:, 2], y * 10 ** obs[:, 3] - y), fmt="o", ms=5, mfc="white",
                color="k", lw=1, capsize=2, zorder=5, label="Milky Way (Kennicutt \\& Evans 2012)")
    ax.set_yscale("log")
    ax.set_xlim(0, edges[-1])
    ax.set_xlabel("R [kpc]")
axes[0].set_ylabel(r"$\Sigma_{\rm atomic}$ [M$_\odot$ pc$^{-2}$]")
axes[1].set_ylabel(r"$\Sigma_{\rm mol}$ [M$_\odot$ pc$^{-2}$]")
axes[0].set_ylim(0.3, 40)
axes[1].set_ylim(1e-3, 60)
axes[0].legend(frameon=False, fontsize=8, loc="lower right")
fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=axes, label="t [Myr]", pad=0.02)
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
print("latest profiles [Msun/pc^2]  R:", edges[:-1].astype(int).tolist())
print("  atomic:", np.array2string(atomic[-1], precision=2, floatmode="fixed"))
print("  mol   :", np.array2string(mol[-1], precision=3, floatmode="fixed"))
