"""
Important constants for datasets
"""
import numpy as np
from pathlib import Path
from crystals import Crystal
import skued

CAMERA_LENGTH = 0.2939  # m


class DatasetInfo:
    """
    This class aggregates information specific to each dataset,
    most importantly calibration of lattice vectors.
    """

    def __init__(self, name, path, reciprocal_basis):
        self.name = name
        self.path = Path(path)
        assert self.path.exists()

        # Determine the center twice, from two peaks each
        center010 = (0.5) * (reciprocal_basis[(0, 1, 0)] + reciprocal_basis[(0, -1, 0)])
        center001 = (0.5) * (reciprocal_basis[(0, 0, 1)] + reciprocal_basis[(0, 0, -1)])
        self.center = tuple(np.rint(0.5 * center010 + 0.5 * center001))
        self.bstar = np.array(reciprocal_basis[(0, 1, 0)]) - np.array(self.center)
        self.cstar = np.array(reciprocal_basis[(0, 0, 1)]) - np.array(self.center)

    def miller_to_arrindex(self, h, k, l):
        """ Determine the array indices for Miller indices (h, k, l). """
        return (k * self.bstar + l * self.cstar + np.array(self.center)).astype(int)

    def kgrid(self):
        """ Calculate the kx, ky meshgrid of the dataset. """
        wavelength = skued.electron_wavelength(keV=90)
        # Grid of detector dimention in meters
        pixel_width = 14e-6  # Pixel width of a Gatan Ultrascan 895
        cx, cy = self.center
        extent_x = np.arange(0, 2048) - cx
        extent_y = np.arange(0, 2048) - cy
        xx, yy = np.meshgrid(pixel_width * extent_x, pixel_width * extent_y)

        r, phi = np.sqrt(xx ** 2 + yy ** 2), np.arctan2(yy, xx)
        angle = np.arctan(r / CAMERA_LENGTH)

        # Scattering vector norm (inverse Angs)
        nG = 4 * np.pi * np.sin(angle / 2) / wavelength
        extent_kx, extent_ky = (nG * np.cos(phi))[0, :], (nG * np.sin(phi))[:, 0]
        return np.meshgrid(extent_kx, extent_ky)


class DatasetInfo200(DatasetInfo):
    """
    Proxy class to handle datasets defined with {002} reflections rather than {001}
    """

    def __init__(self, name, path, reciprocal_basis):
        self.name = name
        self.path = Path(path)
        assert self.path.exists()

        # Determine the center twice, from two peaks each
        center020 = (0.5) * (reciprocal_basis[(0, 2, 0)] + reciprocal_basis[(0, -2, 0)])
        center002 = (0.5) * (reciprocal_basis[(0, 0, 2)] + reciprocal_basis[(0, 0, -2)])
        self.center = tuple(np.rint(0.5 * center020 + 0.5 * center002))
        self.bstar = (np.array(reciprocal_basis[(0, 2, 0)]) - np.array(self.center)) / 2
        self.cstar = (np.array(reciprocal_basis[(0, 0, 2)]) - np.array(self.center)) / 2


# -----------------------------------------------------------------------------
# Sample 1 ID 1 Overnight 4 2mj/cm2
# -----------------------------------------------------------------------------

overnight4 = DatasetInfo(
    name="Overnight 4",
    path=Path("data") / "snse" / "overnight4.hdf5",
    reciprocal_basis={
        # Miller indices : array indices
        (0, 0, 1): np.array((641, 1166)),
        (0, -1, 0): np.array((672, 1407)),
        (0, 1, 0): np.array((405, 1174)),
        (0, 0, -1): np.array((435, 1413)),
    },
)

# -----------------------------------------------------------------------------
# Exfoliated 2 Static
# -----------------------------------------------------------------------------

# For static datasets, the sample was much thinner, and {001} reflections are not visible
# Therefore, we use a proxy class to handle {002} reflections instead.
static = DatasetInfo200(
    name="Static",
    path=Path("data") / "snse" / "static.hdf5",
    reciprocal_basis={
        # Miller indices : array indices
        (0, 0, 2): np.array((1324, 860)),
        (0, -2, 0): np.array((1134, 1353)),
        (0, 2, 0): np.array((827, 698)),
        (0, 0, -2): np.array((638, 1195)),
    },
)
