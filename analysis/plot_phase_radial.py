#!/usr/bin/env python
"""Mass-weighted gas phase diagrams (temperature against hydrogen number density) in radial annuli,
stacked over a range of snapshots. Each panel is normalised to the gas mass in its annulus.
Usage: plot_phase_radial.py <run dir> [first snapshot] [last snapshot] [plot.png]"""
import glob
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

MYR = 0.9778         # code time unit in Myr
RHO_CGS = 6.770e-23  # Msun/pc^3 in g/cm^3
MP = 1.6726e-24
ANNULI = [(2, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 15)]  # kpc
xe, ye = np.linspace(-5, 3, 161), np.linspace(0.8, 7.5, 135)

run = sys.argv[1].rstrip("/")
snaps = sorted(glob.glob(run + "/output/snapshot_*.hdf5"))
first = int(sys.argv[2]) if len(sys.argv) > 2 else len(snaps) - 1
last = int(sys.argv[3]) if len(sys.argv) > 3 else len(snaps) - 1
plot = sys.argv[4] if len(sys.argv) > 4 else f"{run}/phase_radial_{first:03d}_{last:03d}.png"

H = np.zeros((len(ANNULI), len(xe) - 1, len(ye) - 1))
mass, cold, mol, neutral = (np.zeros(len(ANNULI)) for _ in range(4))
times = []
for snap in snaps[first:last + 1]:
    with h5py.File(snap, "r") as f:
        g = f["PartType0"]
        x0 = np.median(f["PartType2/Coordinates"][:], axis=0)  # disc centre
        R = np.hypot(*(g["Coordinates"][:] - x0)[:, :2].T) / 1e3
        m, T, Z = g["Masses"][:], g["Temperature"][:], g["Metallicity"][:]
        nH = g["Density"][:] * RHO_CGS * (1 - Z[:, 0] - Z[:, 1]) / MP
        mH0 = m * (1 - Z[:, 0] - Z[:, 1]) * g["NeutralHydrogenAbundance"][:]
        mmol = mH0 * g["MolecularMassFraction"][:]
        times.append(f["Header"].attrs["Time"] * MYR)
    for i, (lo, hi) in enumerate(ANNULI):
        w = (R >= lo) & (R < hi)
        H[i] += np.histogram2d(np.log10(nH[w]), np.log10(T[w]), bins=(xe, ye), weights=m[w])[0]
        mass[i] += m[w].sum()
        cold[i] += m[w][T[w] < 100].sum()
        mol[i] += mmol[w].sum()
        neutral[i] += mH0[w].sum()

H /= mass[:, None, None] * np.diff(xe)[0] * np.diff(ye)[0]  # mass fraction per dex^2
fig, axes = plt.subplots(2, 3, figsize=(9, 5.2), sharex=True, sharey=True, gridspec_kw=dict(wspace=0.04, hspace=0.04))
norm = LogNorm(1e-4, H.max())
for i, (ax, (lo, hi)) in enumerate(zip(axes.flat, ANNULI)):
    im = ax.pcolormesh(xe, ye, np.where(H[i] > 0, H[i], np.nan).T, norm=norm, cmap="Blues", rasterized=True)
    ax.text(0.96, 0.95, f"{lo}--{hi} kpc", transform=ax.transAxes, ha="right", va="top")
    ax.text(0.96, 0.85, f"$T<100$ K: {100 * cold[i] / mass[i]:.1f}\\% of mass\nmolecular: {100 * mol[i] / neutral[i]:.2f}\\% of H\\,I+H$_2$",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
for ax in axes[1]:
    ax.set_xlabel(r"$\log_{10}\, n_{\rm H}$ [cm$^{-3}$]")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\log_{10}\, T$ [K]")
when = f"t = {times[0]:.0f} Myr" if len(times) == 1 else f"t = {times[0]:.0f}--{times[-1]:.0f} Myr, {len(times)} snapshots"
fig.suptitle(when, fontsize=10, y=0.93)
fig.colorbar(im, ax=axes, pad=0.015, label=r"mass fraction per dex$^2$")
fig.savefig(plot, dpi=200, bbox_inches="tight")
print("plot:", plot)
print(f"{when}\n  R [kpc]   T<100 K   100-8000 K   molecular/(HI+H2)")
for i, (lo, hi) in enumerate(ANNULI):
    print(f"  {lo:2d}-{hi:<2d}    {100 * cold[i] / mass[i]:6.2f}%   {100 * H[i][:, (ye[:-1] >= 2) & (ye[:-1] < 3.9)].sum() * np.diff(xe)[0] * np.diff(ye)[0]:6.1f}%      {100 * mol[i] / neutral[i]:6.2f}%")
