# -*- coding: utf-8 -*-
"""
Decompose diffraction patterns into mode populations and store them.
"""
import itertools as it
from collections import defaultdict
from pathlib import Path

import h5py as h5
import npstreams as ns
import numpy as np
import scipy.interpolate as interpolate
import scipy.optimize as opt
from iris import DiffractionDataset
from mkoneph import (
    IN_PLANE_MODES,
    debye_waller_factors,
    mapply,
    one_phonon_structure_factor,
    prepare_modes,
)
from skimage.filters import gaussian
from skimage.transform import rotate
from skued import detector_scattvectors, nfold
from tqdm import tqdm

from dissutils import GRAPHITE_ANGLE, GRAPHITE_CAMERA_LENGTH

INPUT = Path("data") / "graphite"
OUTPUT = INPUT / "populations"
OUTPUT.mkdir(exist_ok=True)

with DiffractionDataset(INPUT / "graphite_time_corrected_iris5.hdf5", mode="r") as dset:
    c, r = np.asarray(dset.center).astype(int)

DECIMATION = 15
GAMMA_RADIUS = 0.45

EXCLUDES = defaultdict(lambda: {"LO2"})
EXCLUDES[0.5] = {"LO2", "LO3"}
EXCLUDES[0.588] = {"LO2", "LO3"}
EXCLUDES[1.5] = {"LO2", "LO3"}
EXCLUDES[1.588] = {"LO2", "LO3"}

QX, QY, _ = detector_scattvectors(
    keV=90,
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6,
    center=(c, r),
)
QQ = np.sqrt(QX ** 2 + QY ** 2)


def gaussian_table(points, values, sigma=5):
    """
    Perform gaussian smoothing on a table of points

    Parameters
    ----------
    points : ndarray, shape (N, 3)
        Every row corresponds to points
    values : ndarray, shape (N,1)
        Value at every
    sigma : float
        Standard deviation for Gaussian kernel
    """
    makx = max([abs(points[:, 0].min()), abs(points[:, 0].max())])
    maxy = max([abs(points[:, 1].min()), abs(points[:, 1].max())])
    kx, ky = np.meshgrid(np.linspace(-makx, makx, 512), np.linspace(-maxy, maxy, 512))

    interpolated = np.squeeze(
        interpolate.griddata(
            points=points[:, 0:2],
            values=values,
            xi=(kx, ky),
            method="linear",
            fill_value=0.0,
        )
    )

    smoothed = gaussian(interpolated, sigma=sigma)

    return interpolate.RegularGridInterpolator(
        points=(kx[0, :], ky[:, 0]),
        values=smoothed,
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    ).__call__(points[:, 0:2])


def extract_scattering(time, q_points):
    """
    Interpolate the scattering pattern onto `q_points`.
    Negative scattering is no doubt due to Debye-Waller factor's
    influence on elastic scattering; we remove it.

    Parameters
    ----------
    time : float
        Time-delay [ps].
    qx, qy : ndarray, shape (X,Y)
        Meshgrid on which the image should be interpolated.

    Returns
    -------
    interpolated : ndarray, shape (X, Y)
        Image interpolated at `qx` and `qy`.
    """
    # Decomposition is very sensitive; results were
    # optimize before Lorentz factor correction.
    # Therefore, we keep it that way.
    with DiffractionDataset(INPUT / "graphite_time_corrected_iris5.hdf5") as dset:
        image = nfold(
            dset.diff_data(time) - dset.diff_eq(),
            mod=6,
            center=(c, r),
            mask=dset.valid_mask,
            fill_value=0,
        )
    image[QQ < 1.5] = 0
    image[:] = rotate(image, angle=GRAPHITE_ANGLE, center=(c, r), mode="reflect")

    gaussian(image, sigma=10, output=image)

    interpolated = interpolate.RegularGridInterpolator(
        points=(QX[0, :], QY[:, 0]),
        values=image.T,
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    ).__call__(q_points[:, 0:2])

    # Negative scattering is no doubt due to Debye-Waller factor's influence on
    # elastic scattering. We remove it.
    interpolated[interpolated < 0] = 0.0
    return interpolated


def optimize(factors, intensities, **kwargs):
    """
    Optimize for the distribution of populations.
    """
    assert factors.shape[0:2] == intensities.shape

    nk, _, nmodes = factors.shape

    final_shape = (nk, nmodes)
    solution = np.empty(shape=final_shape, dtype=float)
    errors = np.empty_like(solution)

    for k_index in range(nk):
        system = factors[k_index]
        sol, *_ = opt.nnls(A=system, b=intensities[k_index], maxiter=100)
        solution[k_index, :] = sol

    # Estimation of the errors can be done in one step because
    # numpy.linalg.pinv can act on stacks of matrices
    #
    # See Eq 5.98
    #   E. L. Robinson, `Data Analysis for Scientists and Engineers` (2016)
    #   Princeton University Press
    errors = np.sqrt(np.abs(np.diagonal(np.linalg.pinv(factors), axis1=1, axis2=2)))

    return solution, errors


def population(kx, ky, time, exclude=tuple()):
    """
    Determine the population breakdown for all modes,
    inside a single Brillouin zone. The results are interpolated
    on the grid `kx` and `ky`

    Parameters
    ----------
    kx, ky : ndarray, shape (N,M)
        meshgrid-like arrays [Inverse angstroms]
    time : float

    Returns
    -------
    pop : dict[str, ndarray]
        Population changes for each in-plane mode
    """
    exclude = set(exclude)
    mode_names = sorted(set(IN_PLANE_MODES) - exclude)

    refls = [-3, -2, -1, 0, 1, 2, 3]
    reflections = tuple(filter(lambda t: t != (0, 0, 0), it.product(refls, refls, [0])))
    modes = prepare_modes(reflections=reflections, decimate=DECIMATION)
    modes = {
        name: mode.filter_gamma(GAMMA_RADIUS) for name, mode in modes.items()
    }  # Too close to the bragg peak, intensity is misleading

    Ms = debye_waller_factors(modes)

    nbz = len(reflections)
    nk = int(modes["LA"].q_points.shape[0] / nbz)
    nmodes = len(IN_PLANE_MODES) - len(exclude)

    extracted_intensities = extract_scattering(time=time, q_points=modes["LA"].q_points)

    intensities_list = np.split(extracted_intensities, indices_or_sections=nbz, axis=0)

    # There is one system of equation per reduced wavevector k
    # Therefore, axes 1 and 2 are along modes and reflections
    # while axis 0 is along wavevectors
    system = np.zeros(shape=(nk, nbz, nmodes), dtype=float)
    intensities = np.zeros(shape=(nk, nbz), dtype=float)

    for mode_index, mode_name in enumerate(mode_names):

        mode = modes[mode_name]
        F1j = np.abs(one_phonon_structure_factor(mode, dw_factors=Ms)) ** 2
        weights = F1j / mode.frequencies.reshape((-1, 1))

        # Perform smoothing
        # this is very expensive...
        weights = gaussian_table(mode.q_points, weights, sigma=5)

        # Split quantities row-wise such that every quantity
        # is related to a single reflection
        weights_list = np.split(weights, indices_or_sections=nbz, axis=0)

        for refl_index, _ in enumerate(reflections):
            system[:, refl_index, mode_index] = np.squeeze(weights_list[refl_index])
            intensities[:, refl_index] = np.squeeze(intensities_list[refl_index])

    pop, _ = optimize(system, intensities)

    result = dict()
    for index, mode_name in enumerate(mode_names):
        result[mode_name] = pop[:, index]

    return result


if __name__ == "__main__":

    with DiffractionDataset(
        INPUT / "graphite_time_corrected_iris5.hdf5", mode="r"
    ) as dset:
        times = dset.time_points

    # Roundabout way of getting k-points
    modes = prepare_modes(reflections=tuple([(0, 0, 0)]), decimate=DECIMATION)
    modes = {
        name: mode.filter_gamma(GAMMA_RADIUS) for name, mode in modes.items()
    }  # Too close to the bragg peak, intensity is misleading

    kpoints = modes[
        "LA"
    ].q_points  # This only works if the modes were calculated using (0,0,0) reflection only

    kx, ky, _ = np.hsplit(kpoints, 3)

    with h5.File(OUTPUT / "population_timeseries.hdf5", mode="w") as f:
        f.attrs["times"] = times

        f.create_dataset("kx", data=kx)
        f.create_dataset("ky", data=ky)

        for mode_name, mode in modes.items():
            f.create_dataset(
                mode_name,
                shape=(len(kx), len(times)),
                dtype=float,
                chunks=True,
                shuffle=True,
                compression="gzip",
                compression_opts=9,
            )

            f.create_dataset(
                f"{mode_name}_frequencies",
                data=mode.frequencies,
            )

        for time_index, time in enumerate(tqdm(times)):
            pop = population(kx=kx, ky=ky, time=time, exclude=EXCLUDES[round(time, 3)])
            for mode_name in tqdm(IN_PLANE_MODES):
                f[mode_name][:, time_index] = np.squeeze(
                    pop.get(mode_name, np.zeros_like(kx))
                )
