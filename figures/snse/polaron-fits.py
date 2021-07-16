from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter

from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors, tag_axis

DATADIR = Path("data") / "snse"


def polaron(q, A, rp):
    return A * q * rp ** 2 * np.exp(-(q ** 2 * rp ** 2) / 2)


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
    p0=(8e-3, 12),
    bounds=([0, 0], [np.inf, np.inf]),
)

paramss, _ = opt.curve_fit(
    polaron,
    rs,
    ps,
    sigma=es,
    bounds=([-np.inf, 0], [np.inf, np.inf]),
)

figure, ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 3.5))

for (r, p, e), params, marker, color, label in zip(
    [(rf, pf, ef), (rs, ps, es)],
    [paramsf, paramss],
    ["^", "o"],
    discrete_colors(2),
    ["$\\tau=1$ ps", "$\\tau=5$ ps"],
):
    p_smooth = gaussian_filter(p, sigma=3)
    ax.plot(
        r,
        p_smooth,
        marker=marker,
        markerfacecolor="none",
        color=color,
        markersize=6,
        linestyle="None",
        label=label,
    )

    ax.fill_between(
        r,
        y1=p_smooth + e,
        y2=p_smooth - e,
        facecolor=color,
        edgecolor="none",
        linewidth=1,
        alpha=0.2,
    )

    r_ = np.linspace(rf.min(), rs.max(), 1024)
    ax.plot(r_, polaron(r_, *params), color=color, linestyle="-")

    ax.set_ylabel("$\Delta I / I_0$ [a.u.]")
    ax.set_xlabel("$|\mathbf{k}|$ [1/$\AA$]")

ax.set_xlim([rf.min(), rs.max()])

ax.legend(
    loc="center", ncol=2, bbox_to_anchor=(0.5, 1.05), edgecolor="none", facecolor="none"
)

plt.tight_layout()
