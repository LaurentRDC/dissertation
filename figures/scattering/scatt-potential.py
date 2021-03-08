import matplotlib.pyplot as plt
import numpy as np
from crystals import Crystal, Element, Atom
from skued import electrostatic, affe
from warnings import simplefilter
from plotutils import FIGURE_WIDTH, tag_axis, discrete_colors

# For division by zero warnings
simplefilter("ignore", category=RuntimeWarning)

ELEMENTS = [
    "C",
    "Cu",
    "Sn",
    "Au",
]
COLORS = discrete_colors(len(ELEMENTS))


def potential(element, r):
    """ Radial electrostatic potential for an element [V] """
    return electrostatic(
        Crystal(
            unitcell=[Atom(element, coords=(0, 0, 0))],
            lattice_vectors=2 * r.max() * np.eye(3),
        ),
        x=r,
        y=np.zeros_like(r),
        z=np.zeros_like(r),
    )


fig, (ax, ax_r) = plt.subplots(1, 2, figsize=(FIGURE_WIDTH, 3))

r = np.linspace(0.0, 0.5, num=128)
q = np.linspace(0, 12, num=128) / 4 * np.pi

handles = list()
for e, c in zip(ELEMENTS, COLORS):
    (l,) = ax.plot(r, potential(e, r), color=c, label=e)
    ax_r.plot(q, affe(e, q), color=c)
    handles.append(l)

fig.legend(
    handles=handles,
    loc="center",
    ncol=len(ELEMENTS),
    bbox_to_anchor=(0.5, 0.925),
    edgecolor="none",
)

ax.set_xlabel("Radius [$\AA$]")
ax.set_ylabel("Potential [V]")
ax.set_xlim([r.min(), r.max()])
ax.set_ylim([-10, 10000])

ax_r.set_xlabel("Scattering vector $\mathbf{q}$ [$\AA^{-1}$]")
ax_r.set_ylabel("Scattering amplitude")
ax_r.set_xlim([q.min(), q.max()])
ax_r.yaxis.tick_right()
ax_r.yaxis.set_label_position("right")

for ax_, tag in zip([ax, ax_r], ["a) $V_a(\mathbf{x})$", "b) $f_e(\mathbf{q})$"]):
    tag_axis(ax_, text=tag, x=0.95, horizontalalignment="right")

plt.tight_layout()
plt.subplots_adjust(top=0.85)
