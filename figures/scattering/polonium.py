import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from crystals import Crystal, Atom
from skued import electrostatic, affe
from plotutils import LARGE_FIGURE_WIDTH, CBAR_SIZE, named_arrow, tag_axis

from scipy.constants import hbar, m_e

Pu = Crystal.from_database("Pu-epsilon")

CRYSTAL = Crystal(
    unitcell=[Atom("Pu", [0, 0, 0])], lattice_vectors=[[5, 0, 0], [0, 3, 0], [0, 0, 3]]
)
a, b, *_ = CRYSTAL.lattice_parameters

extent = np.linspace(start=-6 * a, stop=6 * a, num=512)
extent_z = np.linspace(start=-1, stop=1, num=32)
xx, yy, zz = np.meshgrid(extent, extent, extent_z)

V = electrostatic(CRYSTAL, x=xx, y=yy, z=zz).mean(axis=2)
np.clip(V, a_min=0, a_max=2000)

# Compute the Fourier transform
# We scale the fourier transform so  that its maximum is 1 (which is always at q=0),
# then scale according to skued's affe
f = np.abs(np.fft.fftshift(np.fft.fft2(V, s=6 * np.array(V.shape))))
f *= affe("Pu", 0) / f.max()
kx = 2 * np.pi * np.fft.fftfreq(f.shape[0], d=abs(xx[1, 1, 1] - xx[0, 0, 0]))
ky = 2 * np.pi * np.fft.fftfreq(f.shape[1], d=abs(yy[1, 1, 1] - yy[0, 0, 0]))

fig, (ax_r, ax_f) = plt.subplots(
    1, 2, figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 1.5)
)
cbar_ax_r = make_axes_locatable(ax_r).append_axes("top", size=CBAR_SIZE, pad=0.05)
cbar_ax_f = make_axes_locatable(ax_f).append_axes("top", size=CBAR_SIZE, pad=0.05)

m_r = ax_r.imshow(
    V,
    vmin=0,
    cmap="CMRmap_r",
    extent=[xx.min(), xx.max(), yy.max(), yy.min()],
    origin="lower",
)
ax_r.set_xlim([-3 * a, 3 * a])
ax_r.set_ylim([-3 * a, 3 * a])
ax_r.yaxis.set_visible(False)
ax_r.set_xlabel("Position [$\AA$]")

m_f = ax_f.imshow(
    f,
    vmin=0,
    cmap="CMRmap_r",
    extent=[kx.min(), kx.max(), ky.min(), ky.max()],
    origin="lower",
)
ax_f.set_xlim([-a, a])
ax_f.set_ylim([-a, a])
ax_f.yaxis.set_visible(False)
ax_f.set_xlabel("Spatial frequency [$\AA^{-1}$]")

for mappable, cbar_ax, label in zip(
    [m_r, m_f],
    [cbar_ax_r, cbar_ax_f],
    ["$V(\mathbf{x})$ [V]", "$f(\mathbf{q})$ [$\AA$]"],
):
    plt.colorbar(mappable=mappable, cax=cbar_ax, orientation="horizontal")
    cbar_ax.set_xlabel(label)
    cbar_ax.xaxis.set_label_position("top")
    cbar_ax.xaxis.tick_top()

# lattice vectors
arrow_kwds = dict(x=0, y=0, length_includes_head=True, width=0.001, fc="k")
named_arrow(
    ax_r,
    dx=CRYSTAL.lattice_parameters[0],
    dy=0,
    text=r"$\mathbf{a}_1$",
    toffset=(0, -0.1),
    tkwds=dict(va="top", ha="center"),
    head_width=0.3,
    **arrow_kwds
)
named_arrow(
    ax_r,
    dx=0,
    dy=CRYSTAL.lattice_parameters[1],  # vertical
    text=r"$\mathbf{a}_2$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    head_width=0.3,
    **arrow_kwds
)

named_arrow(
    ax_f,
    dx=CRYSTAL.reciprocal.lattice_parameters[0],
    dy=0,
    text=r"$\mathbf{b}_1$",
    toffset=(0, -0.1),
    tkwds=dict(va="top", ha="center"),
    head_width=0.1,
    **arrow_kwds
)
named_arrow(
    ax_f,
    dx=0,
    dy=CRYSTAL.reciprocal.lattice_parameters[1],
    text=r"$\mathbf{b}_2$",
    toffset=(-0.1, 0),
    tkwds=dict(va="center", ha="right"),
    head_width=0.1,
    **arrow_kwds
)

tag_axis(ax_r, "a)")
tag_axis(ax_f, "b)")

plt.tight_layout()
