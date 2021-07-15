"""
Fit SnSe's electronic structure from:
Thermoelectric Figure-of-Merit of Fully Dense Single-Crystalline SnSe
(DOI:10.1021/acsomega.8b03323).
"""

import numpy as np
import scipy.optimize as opt
from skued import gaussian
import skimage.filters as filters
from pathlib import Path

DATADIR = Path("data") / "snse" / "estructure"


def loadcsv(path):
    """Load, sort and symmetrize by |k| = 0"""
    k, en = np.loadtxt(path, delimiter=",", comments="#", unpack=True)
    # Sort in ascending |k|
    idx = np.argsort(k)
    k = np.take(k, idx)
    en = np.take(en, idx)
    # Symmetrize about |k|=0
    return np.concatenate((-k[::-1], k)), np.concatenate((en[::-1], en))


valence_g2y_k, valence_g2y = loadcsv(DATADIR / "valence_gamma_y.csv")
valence_g2z_k, valence_g2z = loadcsv(DATADIR / "valence_gamma_z.csv")
conduct_g2y_k, conduct_g2y = loadcsv(DATADIR / "conduction_gamma_y.csv")
conduct_g2z_k, conduct_g2z = loadcsv(DATADIR / "conduction_gamma_z.csv")


def valence(ky, kz, args):

    (ay, cy, sy, oy, az, cz, sz, oz, ag, sg, og) = args
    e_valence = np.zeros_like(ky, dtype=float)

    for a, cy_, cz_, s, o in [
        (ay, cy, 0, sy, oy),
        (ay, -cy, 0, sy, oy),
        (az, 0, cz, sz, oz),
        (az, 0, -cz, sz, oz),
    ]:
        e_valence += a * gaussian(coordinates=[ky, kz], center=(cy_, cz_), std=s) + o

    # Gamma is special
    # It does not "add" on top of the other valleys, but rather is slotted in
    center = (ky / 0.334) ** 2 + (kz / 0.368) ** 2 < 1
    e_valence[center] = (
        ag * (gaussian(coordinates=[ky, kz], center=(0, 0), std=sg) + og)[center]
    )

    return e_valence


def conduction(ky, kz, args):

    (ay, cy, sy, oy, az, cz, sz, oz, agy, agz, og) = args
    e_conduction = np.zeros_like(ky, dtype=float)

    for a, cy_, cz_, s, o in [
        (ay, cy, 0, sy, oy),
        (ay, -cy, 0, sy, oy),
        (az, 0, cz, sz, oz),
        (az, 0, -cz, sz, oz),
    ]:
        e_conduction += a * gaussian(coordinates=[ky, kz], center=(cy_, cz_), std=s) + o

    # Gamma is special; it looks more like a parabola
    # It does not "add" on top of the other valleys, but rather is slotted in
    center = (ky / 0.477) ** 2 + (kz / 0.563) ** 2 < 1
    e_conduction[center] = ((agy * ky ** 2 + agz * kz ** 2) + og)[center]

    return e_conduction


def objective_valence(args):
    """
    Objective function to fit the valence band.
    """
    val_y = valence(ky=valence_g2y_k, kz=np.zeros_like(valence_g2y_k), args=args)
    val_z = valence(ky=np.zeros_like(valence_g2z_k), kz=valence_g2z_k, args=args)

    return np.sum((val_y - valence_g2y) ** 2) + np.sum((val_z - valence_g2z) ** 2)


def objective_conduct(args):
    """
    Objective function to fit the conduction band.
    """
    con_y = valence(ky=conduct_g2y_k, kz=np.zeros_like(conduct_g2y_k), args=args)
    con_z = valence(ky=np.zeros_like(conduct_g2z_k), kz=conduct_g2z_k, args=args)

    return np.sum((con_y - conduct_g2y) ** 4) + np.sum((con_z - conduct_g2z) ** 4)


if __name__ == "__main__":
    positive = (0.001, None)
    negative = (None, -0.001)
    anything = (None, None)
    valence_results = opt.minimize(
        objective_valence, x0=(1, 0.6, 0.3, -0.5, 1, 3 / 4, 0.2, -0.5, 1, 0.4, -1)
    )

    conduct_results = opt.minimize(
        objective_conduct,
        x0=(-0.03, 0.67, 0.1, 1, -0.04, 0.75, 0.15, 0.33, 1.2, 1.1, 0.8),
        bounds=[
            # Along Y
            negative,
            (0.5, 1),
            positive,
            anything,
            # Along Z
            (-0.0401, -0.0399),
            (0.749, 0.751),
            (0.149, 0.151),
            (0.3299, 0.3301),
            # Near Gamma
            (1, 1.5),
            (1, 1.5),
            (0.75, 0.85),
        ],
        method="SLSQP",
    )

    ky, kz = np.meshgrid(np.linspace(-1, 1, num=92), np.linspace(-1, 1, num=92))
    conduction_band = conduction(ky, kz, args=conduct_results.x)
    valence_band = valence(ky, kz, valence_results.x)

    np.save(DATADIR / "ky.npy", ky)
    np.save(DATADIR / "kz.npy", kz)
    # We filter a little because of the hard boundary near Gamma
    np.save(
        DATADIR / "inplane-conduction.npy", filters.gaussian(conduction_band, sigma=1)
    )
    np.save(DATADIR / "inplane-valence.npy", filters.gaussian(valence_band, sigma=1))
