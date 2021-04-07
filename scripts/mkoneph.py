# -*- coding: utf-8 -*-
"""
Calculate the in-plane one-phonon structure factors for Graphite.

The calculation of those is required for multiple plots. 
"""
import multiprocessing as mp
from functools import lru_cache, partial
import json
from math import isclose
from os import cpu_count
from pathlib import Path

import npstreams as ns
import numpy as np
from crystals import Crystal
from crystals.affine import change_of_basis
from crystals.parsers import PWSCFParser
from scipy.constants import physical_constants
from scipy.interpolate import griddata
from skimage.filters import gaussian
from skued import affe, detector_scattvectors
from tqdm import tqdm

INPUT = Path(__file__).parent.parent / "data" / "graphite"
OUTPUT = INPUT / "oneph"
OUTPUT.mkdir(exist_ok=True)

# Mode ordering of graphite according to the file
# Gra-C_XDM_mode_grid_new2.json
MODE_ORDERING = {
    "LA": 0,
    "TA": 1,
    "ZA": 2,
    "LO1": 3,
    "LO2": 4,
    "LO3": 5,
    "TO1": 6,
    "TO2": 7,
    "TO3": 8,
    "ZO1": 9,
    "ZO2": 10,
    "ZO3": 11,
}
MODES = sorted(MODE_ORDERING.keys())
IN_PLANE_MODES = sorted(set(MODE_ORDERING.keys()) - {"ZA", "ZO1", "ZO2", "ZO3"})

NCORES = cpu_count() - 1
EPS = np.finfo(float).eps


def ncells(crystal):
    """Calculate the number of unit cells N present in a crystal of graphite,
    500um x 500um x 50nm"""
    side_length = 5000000  # 500 um in angstroms
    depth = 500  # 50 nm in angstroms
    return int(side_length * side_length * depth / crystal.volume)


def get_eigenvector_frequency(data):
    eigenvector = []

    for key in data:
        eigenvectors = np.zeros((3, 1), dtype=np.complex128)
        if key == "freq":
            frequency = float(data[key])
        else:
            eigenvectors[0] = np.array(float(data[key][0]) + float(data[key][1]) * 1j)
            eigenvectors[1] = np.array(float(data[key][2]) + float(data[key][3]) * 1j)
            eigenvectors[2] = np.array(float(data[key][4]) + float(data[key][5]) * 1j)
            eigenvector.append(eigenvectors)

    return eigenvector, frequency


def extract_info_ordered(fname, crystal, **kwargs):
    """
    Extract the information from JSON file

    Returns
    -------
    q_points : ndarray, shape (N, 3)
    frequencies: ndarray, shape (N, nmodes)
    polarizations: ndarray, shape (N, nmodes, natoms, 3)
    """
    # Reciprocal lattice vectors in units of 2pi/alat
    # so that the output q-points are indeed in "fractional coordinates"
    p = PWSCFParser(crystal.source)
    transf = np.linalg.inv(np.array(p.reciprocal_vectors_alat()))

    with open(fname, mode="r") as f:
        modes = json.load(f)

    modes_op, modes_ac = modes["Optical"], modes["Acoustic"]

    q_points = []
    eig_vector = []
    freq = []

    for key in modes_op.keys():

        for key2 in modes_op[key].keys():

            q_pt = np.array(
                [
                    float(modes_op[key][key2]["q_point"][0]),
                    float(modes_op[key][key2]["q_point"][1]),
                    float(modes_op[key][key2]["q_point"][2]),
                ]
            )
            q_points.append(np.array(q_pt).dot(transf))

            e_vec = []
            frequencies = []

            for i in range(3):

                eigenvector, frequency = get_eigenvector_frequency(
                    modes_ac[key][key2][str(i)]
                )
                e_vec.append(eigenvector)
                frequencies.append(frequency)

            for i in range(3, 12):

                eigenvector, frequency = get_eigenvector_frequency(
                    modes_op[key][key2][str(i)]
                )
                e_vec.append(eigenvector)
                frequencies.append(np.array(frequency))

            eig_vector.append(e_vec)
            freq.append(np.array(frequencies))

    q_points = np.asarray(q_points)
    eig_vector = np.squeeze(np.asarray(eig_vector))
    speed_of_light_cm = (
        physical_constants["speed of light in vacuum"][0] * 100
    )  # in cm/s
    freq = np.asarray(freq) * speed_of_light_cm  # frequencies in Hz

    # Certain acoustic modes may have slightly negative frequencies
    # We shift frequencies up so that the minimum is always 0
    freq -= freq.min() + EPS

    return q_points, freq, eig_vector


def rowdot(arr, brr):
    """ Row-wise dot product """
    # This is much, much faster than np.inner for some reason.
    return np.einsum("ij,ij->i", arr, brr).reshape((-1, 1))


def coth(*args, **kwargs):
    """ Hyperbolic cotangent """
    return np.cosh(*args, **kwargs) / np.sinh(*args, **kwargs)


def mapply(matrix, table):
    """ Apply a matrix transformation to a table where every row is considered a vector. """
    return np.transpose(matrix @ table.T)


def unique_by_row_idx(arr):
    """ Return the indices of the unique rows in arr """
    # np.unique can be rather slow when checking for uniqueness in rows
    # The optimization below results in ~3X faster performance
    #
    # A discussion of why is located here:
    #   https://github.com/numpy/numpy/issues/11136
    arr = np.ascontiguousarray(arr)
    arr_row_view = arr.view("|S%d" % (arr.itemsize * arr.shape[1]))
    _, unique_row_indices = np.unique(arr_row_view, return_index=True)
    return unique_row_indices


def unique_by_rows(*arrs):
    """ Filter arrays by the unique rows of the first array """
    unique_row_indices = unique_by_row_idx(arrs[0])
    return tuple(arr[unique_row_indices] for arr in arrs)


def roughly_unique_by_rows(*arrs, decimals, axis=0):
    """ Apply uniqueness rules on an array, based on lower precision. """
    rough = np.copy(arrs[0])
    np.around(rough, decimals=decimals, out=rough)
    return unique_by_rows(rough, *arrs[1:])


def tile_over_rows(*arrs):
    """ Tile arrays over rows time until all arrays have the same number of rows as the first array"""
    nrows = arrs[0].shape[0]

    arrays = [arrs[0]]
    for array in arrs[1:]:
        missing_reps = int(nrows / array.shape[0])
        reps_tuple = [1] * array.ndim
        reps_tuple[0] = missing_reps
        newarr = np.tile(array, reps=tuple(reps_tuple))
        arrays.append(newarr)

    return tuple(arrays)


def is_in_plane(transformation):
    """ Determine if a symmetry transformation is in the a-b plane """
    translation = transformation[0:3, -1]

    if not isclose(translation[2], 0, abs_tol=1e-5):
        return False

    return True


def apply_symops(kpoints, polarizations, crystal, symprec=1e-1):
    """
    Apply symmetry operations to polarizations vectors and q-points

    kpoints : ndarray, shape (N,3)
        Scattering vector within one Brillouin zone
    polarizations : ndarray, shape (N, natoms, 3)
        Complex polarization vectors. Every row is associated with the corresponding row
        in `kpoints`.
    crystal: crystals.Crystal
        Crystal object with the appropriate symmetry.
    symprec : float, optional
        Symmetry-determination precision.
    """
    # Change of basis matrices allow to express
    # transformations in other bases
    to_reciprocal = change_of_basis(
        np.array(crystal.lattice_vectors), np.array(crystal.reciprocal_vectors)
    )
    from_reciprocal = np.linalg.inv(to_reciprocal)
    reciprocal = lambda m: to_reciprocal @ m @ from_reciprocal

    to_frac = change_of_basis(np.eye(3), np.array(crystal.lattice_vectors))
    from_frac = np.linalg.inv(to_frac)
    cartesian = lambda m: from_frac @ m @ to_frac

    transformed_k, transformed_p = list(), list()

    transformations = crystal.symmetry_operations(symprec=symprec)
    in_plane_transformations = list(filter(is_in_plane, transformations))
    for transf in in_plane_transformations:
        rotation = transf[0:3, 0:3]
        transformed_k.append(mapply(reciprocal(rotation), kpoints))

        # Transforming polarizations is a little more complex
        # because of the extra dimension.
        newpols = np.copy(polarizations)
        for atm_index in range(len(crystal)):
            newpols[:, atm_index, :] = mapply(
                cartesian(rotation), polarizations[:, atm_index, :]
            )

        transformed_p.append(newpols)

    return np.vstack(transformed_k), np.vstack(transformed_p)


class Mode:
    """
    Parameters
    ----------
    name : str
        Mode name, e.g. "LA".
    q_points : ndarray, shape (N, 3)
        Table of reciprocal space vectors where the mode is defined.
    frequencies : ndarray, shape (N, 1)
        Mode frequencies at every point in ``k_points`` [Hz]
    polarizations : ndarray, shape (N, 3, natoms), dtype complex
        Complex mode polarization PER ATOM.
    crystal : crystals.Crystal instance
    hkls : ndarray, shape (N, 3), dtype int, optional
        Nearest Bragg-peak associated with each row in k_points.
        Default is the (000) reflection only.
    """

    def __init__(self, name, q_points, frequencies, polarizations, crystal, hkls=None):
        if hkls is None:
            hkls = np.zeros_like(q_points)

        self.name = name
        self.q_points = q_points
        self.frequencies = frequencies
        self.polarizations = polarizations
        self.crystal = crystal
        self.hkls = hkls

    def save(self, fname):
        """ Save all mode information """
        np.savez(
            fname,
            q_points=self.q_points,
            frequencies=self.frequencies,
            polarizations=self.polarizations,
            hkls=self.hkls,
        )

    def k_points(self):
        """ Determine the unique k-points in this mode. """
        from_miller = change_of_basis(
            np.array(self.crystal.reciprocal_vectors), np.eye(3)
        )
        bragg = mapply(from_miller, self.hkls)
        return self.q_points - bragg

    def filter_gamma(self, radius):
        """ Filter information so that k-points near Gamma are removed. """
        not_near_gamma = np.greater(np.linalg.norm(self.k_points(), axis=1), radius)

        return Mode(
            name=self.name,
            q_points=self.q_points[not_near_gamma],
            frequencies=self.frequencies[not_near_gamma],
            polarizations=self.polarizations[not_near_gamma],
            crystal=self.crystal,
            hkls=self.hkls[not_near_gamma],
        )


def remove_by_rows(remove, *arrs):
    """
    Filter arrays where rows in `remove` are True are removed.
    """
    for arr in arrs:
        assert remove.shape[0] == arr.shape[0]

    return tuple(arr[np.logical_not(remove)] for arr in arrs)


def symmetrize(mode):
    """
    Extend mode information by symmetrization.

    The input is assumed to be a Mode representing information in
    a single Brillouin zone (000), unsymmetrized.
    """
    assert np.allclose(mode.hkls, np.zeros_like(mode.q_points))

    k_points_frac = mode.q_points  # called k_points because only one brillouin zone
    polarizations = mode.polarizations
    frequencies = mode.frequencies
    cryst = mode.crystal

    # Conversion matrices
    to_fractional = change_of_basis(np.eye(3), np.array(cryst.reciprocal_vectors))
    from_fractional = np.linalg.inv(to_fractional)

    # Apply symmetry information and tile arrays to same shape
    k_points_frac, polarizations = apply_symops(
        k_points_frac, polarizations, crystal=cryst
    )
    k_points_frac, polarizations, frequencies = tile_over_rows(
        k_points_frac, polarizations, frequencies
    )

    # Change of basis to inverse angstroms
    # Called k_points because still inside a single Brillouin zone.
    k_points = mapply(from_fractional, k_points_frac)

    return Mode(
        name=mode.name,
        q_points=k_points,
        frequencies=frequencies,
        polarizations=polarizations,
        crystal=cryst,
        hkls=mode.hkls,
    )


def extend_bragg(mode, reflections):
    """
    Expand mode information so that it covers the entire detector range.

    The input is assumed to be a Mode representing information in a single Brillouin zone (000).
    Symmetry operations are applied, as well as some filtering.
    """
    q_points = mode.q_points
    polarizations = mode.polarizations
    frequencies = mode.frequencies

    over_reflections, hkls = list(), list()
    astar, bstar, cstar = mode.crystal.reciprocal_vectors
    for (h, k, l) in reflections:
        H = h * astar + k * bstar + l * cstar
        over_reflections.append(q_points + H[None, :])

        # Quick way to copy a chunk of h, k, l row-wise
        hkls.append(np.zeros_like(q_points) + np.array([h, k, l], dtype=int)[None, :])

    # Important to distinguish
    q_points = np.vstack(over_reflections)
    hkls = np.vstack(hkls)
    q_points, polarizations, frequencies, hkls = tile_over_rows(
        q_points, polarizations, frequencies, hkls
    )

    return Mode(
        name=mode.name,
        q_points=q_points,
        frequencies=frequencies,
        polarizations=polarizations,
        crystal=mode.crystal,
        hkls=hkls,
    )


def decimate(mode, decimals=2):
    """
    Decimate the information contained in modes based on similarity
    i.e. q-points that are roughly the same will be purged.
    """
    # Filter qs to a slightly lower precision
    # Because of QE rounding error + symmetry operations,
    # lots of duplicated points...
    q_points, polarizations, frequencies, hkls = roughly_unique_by_rows(
        mode.q_points,
        mode.polarizations,
        mode.frequencies,
        mode.hkls,
        decimals=decimals,
    )

    return Mode(
        name=mode.name,
        q_points=q_points,
        polarizations=polarizations,
        frequencies=frequencies,
        hkls=hkls,
        crystal=mode.crystal,
    )


def prepare_modes(reflections, decimate=3):
    """
    Prepare all modes for further calculations. Caching included.

    Parameters
    ----------
    reflections : iterable of 3-tuples
        Miller indices of the reflections to consider.
    decimate : int, optional
        Decimation number, i.e., only one in every `decimate` q-points will be kept. Increase this factor to
        lower the number of q-points considered in the calculations.
    temperatures : dict[str, float] or None, optional
        Mode temperatures [K], e.g. {"LA": 100}. Default value is room temperature.

    Returns
    -------
    modes : dict[str, Mode]
    """
    return _prepare_modes(reflections=tuple(reflections), decimate=int(decimate))


@lru_cache(maxsize=16)
def _prepare_modes(reflections, decimate=3):
    cryst = Crystal.from_pwscf(INPUT / "graphite.out")

    k_points, frequencies, polarizations = extract_info_ordered(
        INPUT / "Gra-C_XDM_mode_grid_new2.json", crystal=cryst
    )

    # NOTE
    # Optimization trick is to skip over a few k-points
    # This doesn't seem to mangle the output of this program
    # and it saves time.
    k_points = k_points[0::decimate, :]
    frequencies = frequencies[0::decimate, :]
    polarizations = polarizations[0::decimate, :]

    # Int -> Str for in-plane modes
    names = {v: k for k, v in MODE_ORDERING.items()}
    modes = [
        Mode(
            name=names[mode_index],
            q_points=k_points,
            frequencies=frequencies[:, mode_index],
            polarizations=polarizations[:, mode_index, :, :],
            crystal=cryst,
        )
        for mode_index in range(frequencies.shape[1])
    ]

    with mp.Pool(NCORES) as pool:
        symmetrized_modes = pool.map(symmetrize, modes)
        extended_modes = pool.map(
            partial(extend_bragg, reflections=reflections), symmetrized_modes
        )

    return {mode.name: mode for mode in extended_modes}


def phonon_amplitude(frequencies, temperature):
    """
    Phonon amplitude within the Debye-Waller factor

    Parameters
    ----------
    frequencies : ndarray, shape (N,)
        Frequencies of a single phonon mode [Hz].
    temperature : float
        Mode temperature [K].

    Returns
    -------
    amplitude : ndarray, shape (N,)
        Amplitude SQUARED

    References
    ----------
    Xu and Chiang (2005) Eq. 23
    """
    # Factor of 1/N is calculated in the parent caller
    hbar = physical_constants["Planck constant over 2 pi in eV s"][0]
    kB = physical_constants["Boltzmann constant in eV/K"][0]
    return (hbar / frequencies) * coth(hbar * frequencies / (2 * kB * temperature))


def _debye_waller_factor(modes, temperatures, atm_index):
    """ Calculate a Debye-Waller factor for one atom. """
    # This calculation assumes that the Debye-Waller factor is isotropic
    # i.e. it only depends on the magnitude of q and polarizations ek
    # The anisotropic factor is much more computationally expensive
    n = ncells(modes["LA"].crystal)

    # The sum happens for all q's, but the polarization vectors of a single
    # Brillouin zone. This is annoying to keep track of. Since polarization
    # vectors are the same across different Brillouin zones, we sum over all
    # zones. This leads to "double" counting, hence a correction factor based
    # on the number of Brilluoin zones `nzones`.
    hkls, *_ = unique_by_rows(modes["LA"].hkls)
    nzones = hkls.shape[0]

    def accumulator(mode):
        # This is the sum over k
        # Really, it is a sum over all q, and it will be corrected
        # in the parent function.
        temp = temperatures[mode.name]
        return (
            np.sum(
                (1 / n)
                * phonon_amplitude(mode.frequencies, temp)
                * np.linalg.norm(mode.polarizations[:, atm_index, :], axis=1) ** 2
            )
            / nzones
        )

    # This is the sum over modes
    return ns.sum(accumulator(m) for m in modes.values())


def debye_waller_factors(modes, temperatures=None):
    """
    Compute the debye-waller factor based on all mode information.
    These modes are assumed to have been expanded, i.e. represent mode information
    over the entire detector range.

    For performance reasons, we consider the isotropic Debye-Waller effect.

    Parameters
    ----------
    modes : dict[str, Mode]
    temperatures : dict[str, float] or None, optional
        Mode temperatures [K]. Defaults to room temperature.

    Returns
    -------
    factors : iterable of narrays
        One factor for every atom in the unit cell.

    References
    ----------
    Xu and Chiang (2005) Eq. 19 (anisotropic) and Eq. 20 (isotropic)
    """
    amu_to_kg = 1.6605e-27  # atomic mass units to Kg

    if temperatures is None:
        temperatures = {m: 300 for m in modes.keys()}

    # We loop through atoms in order that they are visible in the PWSCF file
    # That's the tag properties on Atom objects
    atoms = sorted(modes["LA"].crystal, key=lambda a: a.tag)
    q2 = np.linalg.norm(modes["LA"].q_points, axis=1) ** 2
    prefactor = lambda atm: q2 / (12 * atm.mass * amu_to_kg)

    # Parallelizing this calculation is actually slower
    # The correction factor `nzones` accounts for the "double" counting
    # of multiple BZ in the sum
    return [
        prefactor(atom) * _debye_waller_factor(modes, temperatures, index)
        for index, atom in enumerate(atoms)
    ]


def one_phonon_structure_factor(mode, dw_factors):
    """
    Compute the one-phonon structure factor associated with a mode.

    Parameters
    ----------
    mode : Mode
        Mode defined at `N` q-points.
    dw_factors : iterable of ndarray, shapes (N,)
        Debye-Waller factors, at every q-point of `mode`, for each atom in the unit cell.

    Returns
    -------
    oneph: ndarray, shape (N,)
        One-phonon structure factor for `mode`.
    """
    qpoints, polarizations, hkls, crystal = (
        mode.q_points,
        mode.polarizations,
        mode.hkls,
        mode.crystal,
    )

    assert dw_factors[0].shape == (qpoints.shape[0],)
    assert qpoints.shape == hkls.shape
    assert polarizations[:, 0, :].shape == qpoints.shape

    q_norm = np.linalg.norm(qpoints, axis=1, keepdims=True)

    # We loop through atoms in order that they are visible in the PWSCF file
    atoms = sorted(crystal, key=lambda a: a.tag)
    accumulator = np.zeros(shape=(qpoints.shape[0], 1), dtype=complex)
    for atm_index, atm in enumerate(atoms):
        # Accumulator is built in pieces
        # because all of these arrays are pretty big
        arg = np.ones_like(
            accumulator, dtype=complex
        )  # because polarization are complex vectors

        arg *= np.exp(-1 * dw_factors[atm_index].reshape(-1, 1))
        arg *= affe(atm, q_norm) / np.sqrt(atm.mass)
        arg *= rowdot(qpoints, polarizations[:, atm_index, :])
        accumulator += arg

    return np.nan_to_num(accumulator)


def render(mode_str, reflections, smoothing_sigma=15):
    """
    Render the one-phonon structure factor map as visible on
    the Siwick research group detector, for a specific mode.

    Parameters
    ----------
    mode_str : str
        Mode name, e.g. "LA"
    reflections : iterable of 3-tuple
        Reflections to use in the render.
    smoothing_sigma : int, optional
        Size in pixel of the smoothing gaussian kernel.

    Returns
    -------
    image : ndarray, shape (2048, 2048)
        One-phonon structure factor amplitude squared
    qx, qy : ndarray, shape (2048, 2048)
        Q-point mesh on which the one-phonon structure factor was calculated.
    f : ndarray, shape (2048, 2048)
        Mode frequencies [meV]
    """
    reflections = tuple(reflections)  # need to be hashable for caching to work
    modes = prepare_modes(reflections)
    Ms = debye_waller_factors(modes)

    mode = modes[mode_str]

    F1j = np.abs(one_phonon_structure_factor(mode, dw_factors=Ms)) ** 2

    # Create a grid of wavevectors that are visible on the detector
    # Also determine what reflections are visible on the detector
    qx, qy, _ = detector_scattvectors(
        keV=90,
        camera_length=0.25,
        shape=(2048, 2048),
        pixel_size=14e-6,
        center=(1024, 1024),
    )

    # Interpolation sucks
    # Here is an idea for further performance enhancements:
    #    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids?noredirect=1
    #
    # Because the number of reciprocal space points is so large, we can do with only nearest interpolation
    interpolated_oneph = griddata(
        points=mode.q_points[:, 0:2],
        values=F1j,
        xi=(qx, qy),
        method="nearest",
        fill_value=0.0,
    )

    # Frequencies are required for the oneph-majority figure
    interpolated_f = griddata(
        points=mode.q_points[:, 0:2],
        values=mode.frequencies.reshape((-1, 1)),
        xi=(qx, qy),
        method="nearest",
        fill_value=0.0,
    )

    image = np.squeeze(interpolated_oneph)
    return gaussian(image, sigma=smoothing_sigma), qx, qy, np.squeeze(interpolated_f)


def calculate(mode_str, reflections):
    """
    Plot the one-phonon structure factor map as visible on
    the Siwick research group detector, for a specific mode,
    on a Matplotlib `Axes` object.

    Parameters
    ----------
    mode_str : str
        Mode name, e.g. "LA"
    reflections : iterable of 3-tuple
        Reflections to use in the render.
    """
    # Calculate the locations of all reflections
    # that we  can plot as scatter.
    # Note that this is different than the hkls arrays store in the Mode class
    modes = prepare_modes(
        reflections
    )  # it's ok to compute this here since results are cached.
    cryst = modes["LA"].crystal
    astar, bstar, cstar = cryst.reciprocal_vectors
    bragg_peaks = np.vstack(
        [h * astar + k * bstar + l * cstar for (h, k, l) in reflections]
    )

    image, qx, qy, f = render(mode_str, reflections=reflections, smoothing_sigma=15)
    np.save(OUTPUT / f"{mode_str}_oneph.npy", image)
    np.save(OUTPUT / f"{mode_str}_freq.npy", f)
    np.save(OUTPUT / f"qx.npy", qx)
    np.save(OUTPUT / f"qy.npy", qy)
    np.save(OUTPUT / "bragg_peaks.npy", bragg_peaks)


if __name__ == "__main__":
    in_plane_refls = filter(
        lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
    )
    in_plane_refls = tuple(in_plane_refls)

    for mode in tqdm(IN_PLANE_MODES):
        calculate(mode, in_plane_refls)
