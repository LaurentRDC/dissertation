import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter, FixedLocator
import numpy as np
from crystals import Crystal
from skued import electrostatic, affe
from plotutils import FIGURE_WIDTH, ImageGrid, named_arrow, tag_axis

from scipy.constants import hbar, m_e

CRYSTAL = Crystal.from_database("Pu-epsilon")
a, *_ = CRYSTAL.lattice_parameters

extent = np.linspace(start=-3.1 * a, stop=3.1 * a, num=256)
extent_z = np.linspace(start=-1, stop=1, num=32)
xx, yy, zz = np.meshgrid(extent, extent, extent_z)

V = electrostatic(CRYSTAL, x=xx, y=yy, z=zz).mean(axis=2)
np.clip(V, a_min=0, a_max=2000)

# Compute the Fourier transform
# We scale the fourier transform so  that its maximum is 1 (which is always as q=0),
# then scale according to skued's affe
f = np.abs(np.fft.ifftshift(np.fft.ifft2(V, s=4 * np.array(V.shape))))
f *= affe("Pu", 0) / f.max()

l = f.shape[0]
f = f[3 * l // 8 : 5 * l // 8, 3 * l // 8 : 5 * l // 8]
dx = abs(xx[0, 0, 0] - xx[1, 1, 0])
frequencies = np.fft.ifftshift(np.fft.fftfreq(xx.shape[0], d=dx))

fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 1.5))
(ax_r, ax_f) = ImageGrid(
    fig, 111, nrows_ncols=(1, 2), cbar_mode="each", cbar_location="top"
)

for im, ax, label in zip([V, f], [ax_r, ax_f], ["a)", "b)"]):
    # The value of 'extent' needs to be the same for both
    # plots so that the scaling of images is equal
    m = ax.imshow(
        im, vmin=0, cmap="CMRmap_r", extent=[xx.min(), xx.max(), yy.max(), yy.min()]
    )
    cbar = ax.cax.colorbar(m)
    ax.yaxis.set_visible(False)
    tag_axis(ax, label)

# ImageGrid axes share the yaxis
# So we must use the x-axis to show position/frequency scales
ax_r.set_xlabel("Position [$\AA$]")

factor = abs(xx.max() / frequencies.max())
ax_f.xaxis.tick_bottom()
ax_f.xaxis.set_label_position("bottom")
ax_f.xaxis.set_major_locator(
    FixedLocator([-factor * 5, -factor * 2.5, 0, factor * 2.5, factor * 5])
)
ax_f.xaxis.set_major_formatter(FixedFormatter(["-5", "-2.5", "0", "2.5", "5"]))
ax_f.set_xlabel("Spatial frequency [$\AA^{-1}$]")

ax_r.cax.set_xlabel("$V(\mathbf{x})$ [V]")
ax_f.cax.set_xlabel("$f(\mathbf{q})$ [$\AA$]")

# lattice vectors

arrow_kwds = dict(
    x=0, y=0, length_includes_head=True, width=0.001, head_width=0.3, fc="k"
)
named_arrow(
    ax_r,
    dx=CRYSTAL.lattice_parameters[0],
    dy=0,
    text=r"$\mathbf{a}_1$",
    toffset=(0, 0.1),
    tkwds=dict(va="top", ha="center"),
    **arrow_kwds
)
named_arrow(
    ax_r,
    dx=0,
    dy=-CRYSTAL.lattice_parameters[0],  # vertical
    text=r"$\mathbf{a}_2$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    **arrow_kwds
)
a1, a2, _ = CRYSTAL.lattice_vectors
b1, b2, _ = CRYSTAL.reciprocal_vectors
named_arrow(
    ax_f,
    dx=CRYSTAL.lattice_parameters[0] / np.linalg.norm(b1),
    dy=0,
    text=r"$\mathbf{b}_1$",
    toffset=(0, 0.1),
    tkwds=dict(va="top", ha="center"),
    **arrow_kwds
)
named_arrow(
    ax_f,
    dx=0,
    dy=-CRYSTAL.lattice_parameters[0] / np.linalg.norm(b1),  # vertical
    text=r"$\mathbf{b}_2$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    **arrow_kwds
)

plt.tight_layout()
