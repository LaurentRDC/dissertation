# -*- coding: utf-8 -*-
"""
Extract the fast component of the diffuse intensity change in SnSe
as a function of radius.
"""
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import skued
import itertools as it
from crystals import Crystal
from iris import DiffractionDataset
from dissutils.snse import overnight4

DATADIR = Path("data") / "snse"
DATADIR.mkdir(exist_ok=True)

HEADER = """# Fast component of the diffuse intensity change integrated in a ring around zone-center.
# k [1/A], Amplitude [a.u], Error in amplitude [a.u.]
"""

CRYSTAL = Crystal.from_cif(DATADIR / "snse_pnma.cif")

INDICES_DIFFUSE_FAST = [
    (0, -1, 3),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -1, 7),
    (0, -3, 5),
    (0, -5, 3),
]

INDICES_DIFFUSE_SLOW = [
    (0, -1, 3),
    (0, 2, 0),
    (0, -2, 0),
    (0, 0, 4),
    (0, -2, 4),
    (0, -1, 5),
    (0, -3, 5),
]

# FWHM of the gaussian IRF
IRF = 130  # femtoseconds


@skued.with_irf(IRF / 1e3)
def biexponential(time, *args, **kwargs):
    return skued.biexponential(time, 0, *args, **kwargs)


# Fast time-scale -------------------------------------------------------------

radii = range(10, 50, 1)
timeseries = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

    for r in radii:
        inner_radius = r - 11
        outer_radius = r + 11
        timeseries[r] = np.zeros_like(timedelays)
        for indices in INDICES_DIFFUSE_FAST:
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
        method="trf",
        bounds=([-np.inf, 0, 0, 2.5, -np.inf], [0, np.inf, 0.7, np.inf, np.inf]),
        max_nfev=1000,
        ftol=1e-10,
    )
    amp = abs(params[0])
    err = np.sqrt(np.diag(pcov))[0]

    radii_.append(r)
    amplitudes.append(amp)
    amplitudes_err.append(err)

ks = dk * np.asarray(radii_)
amplitudes = np.asarray(amplitudes)
amplitudes_err = np.asarray(amplitudes_err)

_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters

data = np.empty(shape=(2 * len(ks), 3))
data[:, 0] = np.concatenate(
    [ks, np.linspace(ks.max(), np.hypot(bstar / 2, cstar / 2), num=len(ks))]
)
data[:, 1] = np.concatenate([amplitudes, np.zeros_like(ks)])
data[:, 2] = np.concatenate(
    [amplitudes_err, np.full_like(amplitudes_err, fill_value=amplitudes_err.mean() / 2)]
)

np.savetxt(DATADIR / "fast-diffuse-profile.csv", data, delimiter=",", header=HEADER)

# Slow time-scale -------------------------------------------------------------

with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points
    eq = dset.diff_eq()
    t1 = np.argmin(np.abs(timedelays - 4))
    t2 = np.argmin(np.abs(timedelays - 6))
    IMAGE = np.mean(dset.diffraction_group["intensity"][:, :, t1:t2], axis=2)
IMAGE -= eq
IMAGE /= eq


def diffuse_5ps(h, k, l):
    extent = np.linspace(start=-1 / 2, stop=1 / 2, num=64, endpoint=True)
    bz = np.zeros(shape=(len(extent), len(extent)), dtype=float)
    err = np.zeros_like(bz)
    xoffset = yoffset = -5  # Reflections are not perfectly aligned

    for y, z in it.product(extent, repeat=2):
        yi, xi = overnight4.miller_to_arrindex(h, k + y, l + z)
        bz[np.argmin(np.abs(extent - y)), np.argmin(np.abs(extent - z))] = np.mean(
            IMAGE[
                xi - 15 + xoffset : xi + 15 + xoffset,
                yi - 15 + yoffset : yi + 15 + yoffset,
            ]
        )
        err[np.argmin(np.abs(extent - y)), np.argmin(np.abs(extent - z))] = np.std(
            IMAGE[
                xi - 15 + xoffset : xi + 15 + xoffset,
                yi - 15 + yoffset : yi + 15 + yoffset,
            ]
        )

    return bz, err


result = np.zeros(shape=(64, 64))
err = np.zeros_like(result)

for (h, k, l) in INDICES_DIFFUSE_SLOW:
    r, e = diffuse_5ps(h, k, l)
    result += r
    err += e ** 2

result /= len(INDICES_DIFFUSE_SLOW)
err = np.sqrt(err) / len(INDICES_DIFFUSE_SLOW)

_, bstar, cstar, *_ = CRYSTAL.reciprocal.lattice_parameters
kx, ky = np.meshgrid(
    np.linspace(-bstar / 2, bstar / 2, num=result.shape[0]),
    np.linspace(-cstar / 2, cstar / 2, num=result.shape[1]),
)
kk = np.hypot(kx, ky)
_, profile = skued.azimuthal_average(result, center=np.asarray(result.shape) // 2)
_, profile_err = skued.azimuthal_average(err, center=np.asarray(err.shape) // 2)
_, kr = skued.azimuthal_average(kk, center=np.asarray(result.shape) // 2)

data = np.empty(shape=(len(kr), 3))
data[:, 0] = kr
data[:, 1] = profile
data[:, 2] = profile_err

np.savetxt(DATADIR / "slow-diffuse-profile.csv", data, delimiter=",", header=HEADER)
