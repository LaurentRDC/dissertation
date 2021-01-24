"""
Visualization of the electronic structure of graphene

Based on Castro-Neto et al, The electronic properties of graphene, Rev Mod Phys (2009)
"""
import numpy as np
from math import sqrt
from functools import reduce

# Documentation for mplot3d:
#   https://matplotlib.org/api/toolkits/mplot3d.html
from mpl_toolkits.mplot3d import axes3d
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

from plotutils import FIGURE_WIDTH

# Hopping parameter
t = 2.7  # eV
tprime = -0.2 * t


def E(kx, ky):
    a = 1.42  # Angstroms. carbon - carbon distance
    fk = 2 * np.cos(np.sqrt(3) * ky * a) + 4 * np.cos(np.sqrt(3) * ky * a / 2) * np.cos(
        3 * kx * a / 2
    )

    E_plus = t * np.sqrt(3 + fk) - tprime * fk
    E_minus = -t * np.sqrt(3 + fk) - tprime * fk
    return E_plus, E_minus


fig, ax1 = plt.subplots(
    1,
    1,
    figsize=(FIGURE_WIDTH, FIGURE_WIDTH / 1.5),
    # gridspec_kw=dict(width_ratios=[3, 1]),
    subplot_kw=dict(projection="3d", elev=10, azim=-45),
)

extent = np.linspace(-1.7, 1.7, num=256)
kx, ky = np.meshgrid(extent, extent)

kx = kx.flatten()
ky = ky.flatten()
Eplus, Eminus = E(kx, ky)

trim = mtri.Triangulation(kx, ky)  # inside masked
tri = mtri.Triangulation(kx, ky)  # inside not masked
rtri = np.sqrt(
    kx[tri.triangles].mean(axis=1) ** 2 + ky[tri.triangles].mean(axis=1) ** 2
)
imask = rtri < 1.1
omask = rtri > 1.9
trim.set_mask(np.logical_or(imask, omask))
tri.set_mask(omask)


surface_kwargs = dict(
    cmap="plasma",
    vmin=Eminus.min(),
    vmax=6,
    alpha=0.9,
    antialiased=True,
)

ax1.plot_trisurf(trim, Eplus - Eplus.min(), **surface_kwargs)
ax1.plot_trisurf(tri, Eminus - Eplus.min(), **surface_kwargs)

# Draw Brillouin zone
# mplot3d does not respect z-order
# Therefore, the only location to place the BZ is below the dispersion curve
K = np.array([0, -1.7])
R60deg = np.array([[1 / 2, -sqrt(3) / 2], [sqrt(3) / 2, 1 / 2]])
vertices = np.empty((6, 3), dtype=np.float)
for i in range(0, 6):
    point = reduce(lambda m1, m2: m1 @ m2, i * [R60deg], R60deg) @ K
    vertices[i, :] = np.array([point[0], point[1], Eplus.min()])

for i in range(0, 6):
    x1 = vertices[i, 0]
    x2 = vertices[i - 1, 0]

    y1 = vertices[i, 1]
    y2 = vertices[i - 1, 1]
    ax1.plot3D([x1, x2], [y1, y2], zs=Eminus.min() - 1 / 2, color="gray", linewidth=2)
    ax1.plot3D(
        [x1, x1],
        [y1, y1],
        zs=[Eminus.min() - 1 / 2, (Eminus - Eplus.min()).max()],
        color="gray",
        linestyle="dashed",
        linewidth=2,
    )

ax1.plot3D([1.5, 1.5], [-2, 2], zs=Eminus.min() - 1 / 2, color="k", linestyle=":")

ax1.set_xticks([-2, -1, 0, 1, 2])
ax1.set_yticks([-2, -1, 0, 1, 2])
ax1.set_zlim([Eminus.min(), 6])
ax1.set_xlabel(r"$\mathbf{k} \cdot \mathbf{b}_1$ [$\AA^{-1}$]")
ax1.set_ylabel(r"$\mathbf{k} \cdot \mathbf{b}_2$ [$\AA^{-1}$]")
ax1.set_zlabel(r"$E(\mathbf{k})$ [eV]")

plt.subplots_adjust(bottom=0.01, top=0.99)
