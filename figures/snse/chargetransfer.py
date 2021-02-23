"""
Calculate the expected effects of Se -> Sn charge transfer on the Bragg reflections
"""

from enum import Enum
from math import exp
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal, Element, is_element
import skued
from iris import DiffractionDataset
from plotutils.snse_datasets import overnight4

from skued.simulation.form_factors import aspherical_ff
from skued import spectrum_cmap

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")


def pnma_valid(h, k, l):
    """
    True for Miller indices associated with allowed reflections in Pnma space group.
    False otherwise.
    """
    even = lambda n: n % 2 == 0

    return (h == 0) and (
        even(k + l) or (even(k) and (l == 0)) or ((k == 0) and even(l))
    )


class Occupation(Enum):
    Total = "total"
    P1 = "p1"
    P0 = "p0"


def affe_(atom, occ, h, k, l):
    """
    Aspherical atomic form factors for electron scattering.
    Only atoms lighter than Xe (and including) are supported (Z <= 54).

    Parameters
    ----------
    atom : skued.Atom, int, or str
        Atomic number, atomic symbol, or Atom instance.
        Atomic symbols are expected to be properly capitalized, e.g. ``As`` or ``W``.
    occ : Occupation
        Electronic occupation
    h, k, l : int
        Miller indices

    Returns
    -------
    eff : float
        Atomic form factor for electron scattering.

    Raises
    ------
    ValueError : scattering information is not available, for example if the atomic number is larger than 54

    References
    ----------
    .. [#] Jin-Cheng Zheng, Lijun Wu and Yimei Zhu. "Aspherical electron scattering factors and their
           parameterizations for elements from H to Xe" (2009). J. Appl. Cryst. vol 42, pp. 1043 - 1053.
    """
    if isinstance(atom, (int, str)):
        atom = Element(atom)
    element = atom.element

    params_a = aspherical_ff[element][occ.value]["a"]
    params_b = aspherical_ff[element][occ.value]["b"]

    s = np.linalg.norm(CRYSTAL.scattering_vector(np.asarray([h, k, l]))) / (4 * np.pi)
    s2 = s ** 2

    return sum(a * exp(-b * s2) for (a, b) in zip(params_a, params_b))


def intensity_charge_transfer(h, k, l):
    """
    Difference between the diffracted intensity value at equilibrium vs. with
    charge transfer between Sn and Se

    Parameters
    ----------
    h, k, l : int
        Miller indices.

    Returns
    -------
    diff : float
        Fraction of equilibrium intensity
    """
    # Distribute input
    # This works whether G is a list of 3 numbers, a ndarray shape(3,) or
    # a list of meshgrid arrays.
    Gx, Gy, Gz = CRYSTAL.scattering_vector((h, k, l))
    nG = np.sqrt(Gx ** 2 + Gy ** 2 + Gz ** 2)

    se_atoms = filter(is_element("Se"), CRYSTAL)
    sn_atoms = filter(is_element("Sn"), CRYSTAL)

    SFsin, SFcos = 0, 0
    for index, atom in enumerate(se_atoms):
        atomff = skued.affe(atom, nG)
        # Selenium has 4 p electrons, and we transfer only one
        if index == 2:
            atomff -= affe_(atom, Occupation.P1, h, k, l) / 4

        x, y, z = atom.coords_cartesian
        arg = x * Gx + y * Gy + z * Gz

        SFsin += atomff * np.sin(arg)
        SFcos += atomff * np.cos(arg)

    for index, atom in enumerate(sn_atoms):
        atomff = skued.affe(atom, nG)
        # Tin has 2 p electrons, and we transfer only one
        if index == 2:
            atomff += affe_(atom, Occupation.P1, h, k, l) / 2

        x, y, z = atom.coords_cartesian
        arg = x * Gx + y * Gy + z * Gz

        SFsin += atomff * np.sin(arg)
        SFcos += atomff * np.cos(arg)

    eq_intensity = np.abs(skued.structure_factor(CRYSTAL, h, k, l)) ** 2
    transfer_intensity = SFcos ** 2 + SFsin ** 2
    return transfer_intensity - eq_intensity


if __name__ == "__main__":

    # Calculate the normalization and uncertainty
    rj, cj = overnight4.miller_to_arrindex(0, 2, 0)

    with DiffractionDataset(overnight4.path, mode="r") as dset:
        diff_eq = dset.diff_data(0)
        timedelays = dset.time_points
        ts020 = dset.time_series((rj - 15, rj + 15, cj - 15, cj + 15))

        rj_, cj_ = overnight4.miller_to_arrindex(0, 4, 4)
        uncertainty = np.std(dset.time_series((rj_ - 15, rj_ + 15, cj_ - 15, cj_ + 15)))

    I020 = np.mean(diff_eq[rj - 15 : rj + 15, cj - 15 : cj + 15])
    normalization = I020 / np.linalg.norm(skued.structure_factor(CRYSTAL, 0, 2, 0)) ** 2

    figure, ax = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)

    norms = []
    diffs = []
    for (h, k, l) in filter(
        lambda refl: pnma_valid(*refl), CRYSTAL.bounded_reflections(10)
    ):
        if (h, k, l) == (0, 0, 0):
            continue

        norm = np.linalg.norm(CRYSTAL.scattering_vector((h, k, l)))
        diff = abs(normalization * intensity_charge_transfer(h, k, l))
        norms.append(norm)
        diffs.append(diff)

    norms = np.asarray(norms)
    diffs = np.asarray(diffs)

    ax.scatter(
        norms,
        diffs / uncertainty,
        c=spectrum_cmap(norms / norms.max()),
        marker="o",
    )
    ax.axhline(y=1, linestyle="dashed", color="k", linewidth=1)
    ax.set_xlabel(r"$|\mathbf{q}|$ [$\AA^{-1}$]")
    ax.set_ylabel("$\Delta I_{(hkl)}$ [$\sigma_e$]")
