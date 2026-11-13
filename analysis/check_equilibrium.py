#!/usr/bin/env python
"""Equilibrium check of a Milky Way run: radial profiles of the gas and stellar discs in every
snapshot, compared with the first one. Usage: check_equilibrium.py <output dir> [plot.png]"""
import glob
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt

RBINS = np.linspace(0, 20, 41)  # kpc
RMID = 0.5 * (RBINS[1:] + RBINS[:-1])
AREA = np.pi * np.diff(RBINS**2) * 1e6  # pc^2


def profiles(f, ptype):
    g = f[f"PartType{ptype}"]
    x0, v0 = (np.median(f["PartType2"][k][:], axis=0) for k in ("Coordinates", "Velocities"))  # disc centre
    x, v, m = (g["Coordinates"][:] - x0) * 1e-3, g["Velocities"][:] - v0, g["Masses"][:]
    R = np.hypot(x[:, 0], x[:, 1])
    vphi = (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]) / R
    vR = (x[:, 0] * v[:, 0] + x[:, 1] * v[:, 1]) / R
    idx = np.digitize(R, RBINS) - 1
    ok = idx < len(RMID)
    idx, w = idx[ok], m[ok]
    msum = np.bincount(idx, w, len(RMID))
    mean = lambda q: np.bincount(idx, w * q[ok], len(RMID)) / np.maximum(msum, 1e-30)
    vphi_mean = mean(vphi)
    return dict(Sigma=msum / AREA, zrms=1e3 * np.sqrt(mean(x[:, 2] ** 2)), vphi=vphi_mean,
                sigR=np.sqrt(mean(vR**2) - mean(vR) ** 2), sigz=np.sqrt(mean(v[:, 2] ** 2)),
                vR=mean(vR))


def main():
    outdir = sys.argv[1]
    plot = sys.argv[2] if len(sys.argv) > 2 else outdir.rstrip("/") + "/equilibrium.png"
    snaps = sorted(glob.glob(outdir + "/snapshot_*.hdf5"))
    times, gas, stars = [], [], []
    for s in snaps:
        with h5py.File(s, "r") as f:
            times.append(f["Header"].attrs["Time"] * 0.9778)  # Myr
            gas.append(profiles(f, 0))
            stars.append(profiles(f, 2))
    panels = [(gas, "Sigma", r"gas $\Sigma$ [M$_\odot$ pc$^{-2}$]", True),
              (gas, "zrms", r"gas $z_{\rm rms}$ [pc]", False),
              (gas, "vphi", r"gas $\langle v_\phi\rangle$ [km/s]", False),
              (gas, "sigR", r"gas $\sigma_R$ [km/s]", False),
              (stars, "Sigma", r"stellar disc $\Sigma$ [M$_\odot$ pc$^{-2}$]", True),
              (stars, "zrms", r"stellar disc $z_{\rm rms}$ [pc]", False),
              (stars, "sigR", r"stellar disc $\sigma_R$ [km/s]", False),
              (stars, "sigz", r"stellar disc $\sigma_z$ [km/s]", False)]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), sharex=True)
    cmap = plt.get_cmap("viridis")
    for ax, (data, key, label, log) in zip(axes.flat, panels):
        for i, d in enumerate(data):
            ax.plot(RMID, d[key], color=cmap(i / max(len(data) - 1, 1)), lw=1)
        ax.set_ylabel(label)
        if log:
            ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("R [kpc]")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(times[0], times[-1]))
    fig.colorbar(sm, ax=axes, label="t [Myr]", pad=0.01, fraction=0.02)
    fig.savefig(plot, dpi=150, bbox_inches="tight")
    print("plot:", plot)

    # drift of the final snapshot relative to the first, over 2 < R < 16 kpc
    sel = (RMID > 2) & (RMID < 16)
    print(f"t = {times[-1]:.0f} Myr; median (max) |fractional change| over 2 < R < 16 kpc:")
    for data, name in ((gas, "gas"), (stars, "stellar disc")):
        for key in ("Sigma", "zrms", "vphi", "sigR" if name != "gas" else None, "sigz" if name != "gas" else None):
            if key is None:
                continue
            d = np.abs(data[-1][key][sel] / data[0][key][sel] - 1)
            print(f"  {name:12s} {key:6s} {np.median(d):.3f} ({d.max():.3f})")
    print(f"  gas sigma_R (final): median {np.median(gas[-1]['sigR'][sel]):.2f} km/s; "
          f"|<v_R>| max {np.abs(gas[-1]['vR'][sel]).max():.2f} km/s")


if __name__ == "__main__":
    main()
