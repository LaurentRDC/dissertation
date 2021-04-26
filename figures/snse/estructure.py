import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.colors as cl
import numpy as np
from dissutils import LARGE_FIGURE_WIDTH, tag_axis
import skued


def colormap():
    """ Modify the colormap `base` so the last color is completely transparent """
    purples = plt.get_cmap("Purples")
    yellows = plt.get_cmap("YlOrBr_r")

    colors = [purples(i) for i in np.linspace(0, 1, 128)] + [
        yellows(i) for i in np.linspace(0, 1, 128)
    ]
    return cl.LinearSegmentedColormap.from_list(name="", colors=colors)


def electron_dispersion(ky, kz):
    elec_energy = np.zeros_like(ky)

    for (cy, cz), height, width in zip(
        [(2 / 3, 0), (-2 / 3, 0), (0, 0), (0, 2 / 3), (0, -2 / 3)],
        [1.0, 1.0, 0.8, 0.7, 0.7],
        [0.4, 0.4, 0.6, 0.4, 0.4],
    ):
        gauss = skued.gaussian(coordinates=[ky, kz], center=[cy, cz], fwhm=width)
        gauss /= gauss.max()
        elec_energy += height * gauss

    elec_energy /= elec_energy.max()
    elec_energy *= -1
    elec_energy -= elec_energy.min()
    return elec_energy


def hole_dispersion(ky, kz):
    hole_energy = np.zeros_like(ky)

    for (cy, cz), height, width in zip(
        [
            (-2 / 3, 0),
            (2 / 3, 0),
            (0, 0),
            (0, -3 / 4),
            (0, -3 / 4),
            (0, 3 / 4),
            (0, 3 / 4),
        ],
        # Pudding mold band represented by the difference of a tall wide gaussian
        # and a narrow, small gaussian
        [0.7, 0.7, 0.6, 1.0, -0.2, 1.0, -0.2],
        [0.4, 0.4, 0.7, 0.4, 0.08, 0.4, 0.08],
    ):
        gauss = skued.gaussian(coordinates=[ky, kz], center=[cy, cz], fwhm=width)
        gauss /= np.abs(gauss).max()
        hole_energy += height * gauss
    hole_energy /= hole_energy.max()

    hole_energy -= hole_energy.max()
    return hole_energy


ky, kz = np.meshgrid(*[np.linspace(-1, 1, num=256)] * 2)
elec_energy = electron_dispersion(ky, kz) + 0.1
hole_energy = hole_dispersion(ky, kz) - 0.1


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH / 2))
# The small bottom and top rows are a way to control padding of the 2D plot.
# Otherwise, it's hard to have the 2D and 3D plots the same height
gs = gridspec.GridSpec(
    nrows=3,
    ncols=2,
    figure=figure,
    width_ratios=[1.5, 1],
    height_ratios=[0.05, 1, 0.05],
)
ax_3d = figure.add_subplot(gs[:, 0], projection="3d", elev=10, azim=-60)
ax_2d = figure.add_subplot(gs[1, 1])

# 3D plot ---------------------------------------------------------------------
surface_kwds = dict(
    antialiased=True, rcount=128, ccount=128, cmap=colormap(), vmin=-1, vmax=1
)
ax_3d.plot_surface(ky, kz, hole_energy, **surface_kwds)
ax_3d.plot_surface(ky, kz, elec_energy, **surface_kwds)

ax_3d.set_box_aspect((1, 1, 1.2))

locator = ticker.FixedLocator([-1, -0.5, 0, 0.5, 1])
formatter = ticker.FixedFormatter(["-1", "-½", "0", "½", "1"])
for axis in [ax_3d.xaxis, ax_3d.yaxis]:
    axis.set_major_locator(locator)
    axis.set_major_formatter(formatter)

ax_3d.zaxis.set_major_locator(ticker.FixedLocator([-1, -0.5, 0, 0.5, 1]))
ax_3d.zaxis.set_major_formatter(ticker.FixedFormatter(5 * [""]))

ax_3d.set_xlabel(r"$\mathbf{k} / |\mathbf{b}_2|$")
ax_3d.set_ylabel(r"$\mathbf{k} / |\mathbf{b}_3|$")

ax_3d.set_xlim3d([-1, 1])
ax_3d.set_ylim3d([-1, 1])

ax_3d.set_zlabel("$E(\mathbf{k})$ [a.u.]")

# 2D plot ---------------------------------------------------------------------

# Higher density for line plot
ky, kz = np.meshgrid(*[np.linspace(-1, 1, num=1024)] * 2)
elec_energy = electron_dispersion(ky, kz) + 0.2
hole_energy = hole_dispersion(ky, kz) - 0.2

n, m = ky.shape
assert n == m  # assumption for what comes next
line = np.linspace(-1, 1, num=ky.shape[0])
elec_energy_2d = np.zeros_like(line)
hole_energy_2d = np.zeros_like(line)

# Path is Z - Gamma - Y
elec_energy_2d[0 : n // 2] = elec_energy[0 : n // 2, n // 2]
elec_energy_2d[n // 2 : :] = elec_energy[n // 2, n // 2 : :]
hole_energy_2d[0 : n // 2] = hole_energy[0 : n // 2, n // 2]
hole_energy_2d[n // 2 : :] = hole_energy[n // 2, n // 2 : :]

ax_2d.axvline(x=0, color="k", linestyle="--", linewidth=0.5)
ax_2d.scatter(
    line,
    hole_energy_2d,
    s=3,
    c=plt.get_cmap("Purples")(hole_energy_2d - hole_energy_2d.min()),
)
ax_2d.scatter(line, elec_energy_2d, s=3, c=plt.get_cmap("YlOrBr_r")(elec_energy_2d))

ax_2d.set_ylim(-1, 1)
ax_2d.set_xlim([line.min(), line.max()])

ax_2d.xaxis.set_major_locator(ticker.FixedLocator([-1, 0, 1]))
ax_2d.xaxis.set_major_formatter(ticker.FixedFormatter(["$Z$", "$\Gamma$", "$Y$"]))

ax_2d.axhline(y=-0.1, linestyle="dashed", linewidth=0.5, color="k")
ax_2d.text(
    x=0.97, y=-0.1, s="$E_f$", va="bottom", ha="right", transform=ax_2d.transData
)

ax_2d.set_yticks([])

plt.subplots_adjust(
    top=1.0, bottom=0.094, left=0.022, right=0.969, hspace=0.505, wspace=0.026
)
