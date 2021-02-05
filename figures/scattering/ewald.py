import itertools as it

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal, Atom
from matplotlib.ticker import FixedFormatter, FixedLocator
from plotutils import FIGURE_WIDTH, ImageGrid, named_arrow
from skued import lorentzian, indices_to_text, electron_wavelength

EWALD_RADIUS = 2*np.pi / electron_wavelength(keV=100)
EWALD_RADIUS_XRAY = 2*np.pi/0.95 # 13 keV x-rays

ELECTRONS_COLOR = "k"
XRAY_COLOR = "indigo"

# Abstract simple cubic crystal with 3Angs sides
CRYSTAL = Crystal(unitcell=[Atom("C", (0, 0, 0))], lattice_vectors=5 * np.eye(3))

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 1.5))
(ax,) = ImageGrid(fig, rect=111, nrows_ncols=(1, 1), cbar_location="top")

ky, kz = np.meshgrid(np.linspace(-6, 6, 256), np.linspace(-4, 7, 256))
im = np.zeros_like(ky)

for k, l in it.product([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5], repeat=2):
    _, qy, qz = CRYSTAL.scattering_vector((0,k,l))
    im += lorentzian(coordinates=[ky, kz], center=[qy, qz], fwhm=0.15)

m = ax.imshow(
    im,
    vmin=0,
    cmap="CMRmap_r",
    extent=[ky.min(), ky.max(), kz.max(), kz.min()],
)
cbar = ax.cax.colorbar(
    m, ticks=FixedLocator(locs=[0, im.max()]), format=FixedFormatter(["0", "1"])
)
ax.cax.set_xlabel("$|\hat{V}(\mathbf{q})|$ [a.u.]")

# Ewald spheres
ax.add_patch(
    mpatches.Circle(
        xy=(0, EWALD_RADIUS),
        radius=EWALD_RADIUS,
        fc="none",
        ec=ELECTRONS_COLOR,
    )
)
ax.add_patch(
    mpatches.Circle(
        xy=(0, EWALD_RADIUS_XRAY),
        radius=EWALD_RADIUS_XRAY,
        fc="none",
        ec=XRAY_COLOR,
        linestyle="dashed",
    )
)

for (h, k, l) in [(0, 0, 0)]:
    ax.annotate(
        xy=(k, l),
        ha="center",
        va="bottom",
        text=indices_to_text(h, k, l),
        xytext=(k, l + 0.2),
    )

# Lattice vectors
_, _y, _z = CRYSTAL.scattering_vector((0,-4,-2))
arrow_kwds = dict(
    x=_y, y=_z, length_includes_head=True, width=0.001, head_width=0.1, fc="k"
)

named_arrow(
    ax,
    dx=np.linalg.norm(CRYSTAL.reciprocal_vectors[1]),
    dy=0,
    text=r"$\mathbf{b}_2$",
    toffset=(0, -0.1),
    tkwds=dict(va="top", ha="center"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=0,
    dy=np.linalg.norm(CRYSTAL.reciprocal_vectors[2]),
    text=r"$\mathbf{b}_3$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    **arrow_kwds
)


electron_handle = mlines.Line2D(
    [], [], color=ELECTRONS_COLOR, marker=None, linestyle="solid", label="Electrons (100 keV)"
)
xray_handle = mlines.Line2D(
    [], [], color=XRAY_COLOR, marker=None, linestyle="dashed", label="X-rays (13 keV)"
)

fig.legend(
    handles=[electron_handle, xray_handle],
    loc="center",
    ncol=2,
    bbox_to_anchor=(0.5, 0.05),
)

ax.set_xlim([ky.min(), ky.max()])
ax.set_ylim([-2.7, 4.8])
ax.axis("off")
