import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from functools import partial
from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis
from pathlib import Path

DATADIR = Path("data") / "snse"


def polaron(q, A, rp, C):
    return A * q * rp ** 3 * np.exp(-(q ** 2 * rp ** 2) / 4) + C


rf, pf, ef = np.loadtxt(
    DATADIR / "fast-diffuse-profile.csv", delimiter=",", unpack=True
)
pf -= pf.min()
ef /= pf.max()
pf /= pf.max()

rs, ps, es = np.loadtxt(
    DATADIR / "slow-diffuse-profile.csv", delimiter=",", unpack=True
)
ps -= ps.min()
es /= ps.max()
ps /= ps.max()


paramsf, _ = opt.curve_fit(
    polaron,
    rf[rf < 0.32],
    pf[rf < 0.32],
    sigma=ef[rf < 0.32],
    p0=(8e-3, 12, 0),
    bounds=([0, 0, -1], [np.inf, np.inf, 1]),
)

paramss, _ = opt.curve_fit(
    polaron,
    rs,
    ps,
    sigma=es,
    bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]),
)

figure, (ax_fast, ax_slow) = plt.subplots(1, 2, figsize=(LARGE_FIGURE_WIDTH, 2))

for ax, (r, p, e), params, color in zip(
    [ax_fast, ax_slow],
    [(rf, pf, ef), (rs, ps, es)],
    [paramsf, paramss],
    discrete_colors(2),
):
    ax.plot(
        r,
        p,
        marker="o",
        color=color,
        markersize=2,
        linestyle="None",
    )

    ax.fill_between(
        r,
        y1=p + e,
        y2=p - e,
        facecolor="gray",
        edgecolor="k",
        linestyle="dashed",
        linewidth=1,
        alpha=0.2,
    )

    r_ = np.linspace(r.min(), r.max(), 1024)
    ax.plot(r_, polaron(r_, *params), color=color, linestyle="-")

    ax.set_ylabel("$\Delta I / I_0$ [a.u.]")
    ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

    ax.set_xlim([r.min(), r.max()])

ax_slow.yaxis.tick_right()
ax_slow.yaxis.set_label_position("right")

tag_axis(ax_fast, "a)", x=0.9, y=0.9, horizontalalignment="left")
tag_axis(ax_slow, "b)", y=0.9)

plt.tight_layout()
