# -*- coding: utf-8 -*-
"""
Extract the fast component of the diffuse intensity change in SnSe
as a function of radius.
"""
from math import sqrt
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import scipy.stats
import skued
from crystals import Crystal
from iris import DiffractionDataset
from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis
from dissutils.snse import overnight4
from skimage.filters import gaussian

DATADIR = Path("data") / "snse"
DATADIR.mkdir(exist_ok=True)

HEADER = """# Fast component of the diffuse intensity change integrated in a ring around zone-center.
# k [1/A], Amplitude [a.u], Error in amplitude [a.u.]
"""

CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")

INDICES_DIFFUSE = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
    (0, -3, 5),
    (0, -5, 3),
]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, 0, *args, **kwargs)


def polaron(q, A, rp):
    return A * q * rp ** 3 * np.exp(-(q ** 2 * rp ** 2) / 4) / (1 + (q * rp) ** 2) ** 2


if __name__ == "__main__":
    radii = range(10, 70, 1)
    timeseries = dict()
    with DiffractionDataset(overnight4.path, mode="r") as dset:
        timedelays = dset.time_points

        for r in radii:
            inner_radius = r - 5
            outer_radius = r + 5
            timeseries[r] = np.zeros_like(timedelays)
            for indices in INDICES_DIFFUSE:
                q2 = np.linalg.norm(CRYSTAL.scattering_vector(indices)) ** 2

                yi, xi = overnight4.miller_to_arrindex(*indices)
                diffuse_selection = skued.RingSelection(
                    shape=dset.resolution,
                    center=(xi, yi),
                    inner_radius=inner_radius,
                    outer_radius=outer_radius,
                )
                timeseries[r] += dset.time_series_selection(diffuse_selection) / q2

    # Normalize all time-series to pre-time-zero
    for k, ts in timeseries.items():
        timeseries[k] /= np.mean(ts[timedelays < 0])

    kx, _ = overnight4.kgrid()
    dk = kx[0, 1] - kx[0, 0]

    amplitudes = list()
    amplitudes_err = list()
    radii_ = list()
    for (r, ts) in timeseries.items():

        params, pcov = opt.curve_fit(
            biexponential,
            timedelays[timedelays < 15],
            ts[timedelays < 15],
            p0=(-0.01, 0.01, 0.4, 3.7, 1),
        )

        amp = abs(params[0])
        err = np.sqrt(np.diag(pcov))[0]

        if err < amp:
            radii_.append(r)
            amplitudes.append(amp)
            amplitudes_err.append(err)

    ks = dk * np.asarray(radii_)
    amplitudes = np.asarray(amplitudes)
    amplitudes_err = np.asarray(amplitudes_err)

    data = np.empty(shape=(len(ks), 3))
    data[:, 0] = ks
    data[:, 1] = amplitudes
    data[:, 2] = amplitudes_err

    np.savetxt(DATADIR / "fast-diffuse.csv", data, delimiter=",", header=HEADER)
