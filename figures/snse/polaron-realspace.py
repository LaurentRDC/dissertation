from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib
from crystals import Crystal
from dissutils import named_arrow, tag_axis, LARGE_FIGURE_WIDTH

NY = NZ = 5  # Supercell factor in b and c directions
CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")
SNCOLOR, SECOLOR = plt.get_cmap("inferno")([0.2, 0.8])
PATCH_PARAMS = {"Sn": dict(radius=0.4, fc=SNCOLOR), "Se": dict(radius=0.3, fc=SECOLOR)}
COLORMAP = "Purples"  # Important that lowest values map to ~(1,1,1,1) (white)


def signum(x):
    return (x > 0) - (x < 0)


def gaussian2d(x, y, sx, sy):
    """2D Gaussian centered at (0,0)"""
    return np.exp(-(0.5 * (x / sx) ** 2 + 0.5 * (y / sy) ** 2))


def radial(x, y):
    direction = np.array([x, y])
    return direction / np.linalg.norm(direction)


def linear(x, y):
    return np.array([0.1 * signum(float(x)), signum(float(y))])


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 0.4 * LARGE_FIGURE_WIDTH))
gs = gridspec.GridSpec(
    2, 3, figure=figure, height_ratios=(1, 0.1), width_ratios=(1, 1, 0.05)
)

ax_1d = figure.add_subplot(gs[0, 0])
ax_3d = figure.add_subplot(gs[0, 1], sharex=ax_1d)
ax_elements = figure.add_subplot(gs[1, 0:2])
ax_disp = figure.add_subplot(gs[0, -1])

cell = CRYSTAL.supercell(1, NY, NZ)
a, b, c, *_ = cell.lattice_parameters
yy, zz = np.meshgrid(
    np.linspace(-NY / 2, NY / 2, num=256), np.linspace(-NZ / 2, NZ / 2, num=256)
)

for ax, (sy, sz), dirf, amp, label in zip(
    [ax_1d, ax_3d], [(1, 5), (20, 20)], [linear, radial], [2, 1], ["a)", "b)"]
):

    ax.axis("off")
    ax.set_aspect(1)
    im = gaussian2d(yy, zz, sy, sz)
    vmin = max([im[:, 0].max(), im[:, -1].max(), im[0, :].max(), im[-1, :].max()])
    vmax = im.max()
    m = ax.imshow(
        im,
        vmax=vmax + (vmax - vmin) / 3,
        vmin=vmin,
        cmap=COLORMAP,
        extent=[-NY * b / 2, NY * b / 2, -NZ * c / 2, NZ * c / 2],
    )

    line_params = dict(color="k", linewidth=1)
    for i in range(-NY // 2, (NY + 1) // 2):
        ax.axvline(x=i * b, **line_params)
    for j in range(-NZ // 2, (NZ + 1) // 2):
        ax.axhline(y=j * c, **line_params)

    for atm in cell:
        x, y, z = atm.coords_cartesian
        # Center supercell about (0,0)
        y -= NY * b / 2
        z -= NZ * c / 2

        uy, uz = amp * gaussian2d(y, z, sy, sz) * dirf(y, z)

        ax.add_patch(
            Circle(
                xy=(y + uy, z + uz),
                **PATCH_PARAMS[atm.element],
                clip_on=False,
                zorder=x
            )
        )
    ax.set_xlim([-NY * b / 2, NY * b / 2])
    ax.set_ylim([-NZ * c / 2, NZ * c / 2])

    tag_axis(
        ax,
        label,
        x=0,
        y=1.03,
        verticalalignment="bottom",
        horizontalalignment="right",
    )

# Crystal axes
arrow_kwds = dict(
    x=-(NY + 1) * b / 2,
    y=-(NZ + 1) * c / 2,
    clip_on=False,
    width=0.001,
    head_width=0.5,
    head_length=0.5,
    fc="k",
)
named_arrow(
    ax_1d,
    dx=0,
    dy=c,
    toffset=(-0.5, 0),
    text="$\mathbf{c}$",
    tkwds=dict(
        ha="right",
        va="center",
    ),
    **arrow_kwds
)
named_arrow(
    ax_1d,
    dx=b,
    dy=0,
    toffset=(0, -0.5),
    text="$\mathbf{b}$",
    tkwds=dict(
        ha="center",
        va="top",
    ),
    **arrow_kwds
)

# Show element names ----------------------------------------------------------
ax_elements.axis("off")
ax_elements.set_aspect(1)
ax_elements.set_xlim([-1, 1])
ax_elements.set_ylim([-1, 1])
for xy, elem, ha, toffset in zip(
    [(-1, 0), (1, 0)], ["Sn", "Se"], ["right", "left"], [-0.75, 0.75]
):
    ax_elements.add_patch(Circle(xy=xy, **PATCH_PARAMS[elem], clip_on=False))
    ax_elements.text(x=xy[0] + toffset, y=0, s=elem, ha=ha, va="center")

# Displacement colorbar -------------------------------------------------------
plt.colorbar(
    mappable=matplotlib.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=COLORMAP),
    cax=ax_disp,
    orientation="vertical",
)
ax_disp.set_ylabel("$|\mathbf{u}(\mathbf{r})|$ [a.u.]")
ax_disp.yaxis.set_major_locator(ticker.FixedLocator([0, 1]))

plt.subplots_adjust(
    top=0.889, bottom=0.072, left=0.072, right=0.844, hspace=0.1, wspace=0.131
)
