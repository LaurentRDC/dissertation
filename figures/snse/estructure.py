import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.colors as cl
import numpy as np
from dissutils import MEDIUM_FIGURE_WIDTH
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
        [0.7, 0.7, 0.6, 1.0, -0.4, 1.0, -0.4],
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

figure, ax_3d = plt.subplots(
    1,
    1,
    figsize=(MEDIUM_FIGURE_WIDTH, MEDIUM_FIGURE_WIDTH),
    subplot_kw=dict(projection="3d", elev=10, azim=-60),
)

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

plt.tight_layout()
