import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scipy.optimize as opt
import scipy.stats
import scipy.signal as signal
import scipy.interpolate as interpolate
from pathlib import Path
from crystals import Crystal
from plotutils import FIGURE_WIDTH, tag_axis, discrete_colors
from plotutils.snse_datasets import overnight4
import skued
from iris import DiffractionDataset

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [(0, 0, 2), (0, 2, 0), (0, 1, -3), (0, 0, 4)]


figure, axes = plt.subplots(
    4,
    2,
    sharex=True,
    figsize=(FIGURE_WIDTH, FIGURE_WIDTH),
    gridspec_kw=dict(hspace=0.05, wspace=0.025),
)

# Fit to extract width
def peak(f, A, c, fwhm, o):
    return A * skued.gaussian(f, center=c, fwhm=fwhm) + o


fwhms = dict()
with DiffractionDataset(overnight4.path, mode="r") as dset:
    timedelays = dset.time_points

for refl in INDICES_DIFFUSE:

    yi_, xi_ = overnight4.miller_to_arrindex(*refl)
    width = 40
    kx, _ = overnight4.kgrid()
    dk = kx[0, 1] - kx[0, 0]
    k = dk * (np.arange(start=yi_ - width, stop=yi_ + width) - yi_)  # inverse angstroms
    with DiffractionDataset(overnight4.path, mode="r") as dset:
        timedelays = dset.time_points
        cube = dset.diffraction_group["intensity"][
            xi_ - width : xi_ + width,
            yi_ - width : yi_ + width,
            :,
        ]
        cube = cube[:, :, np.less_equal(timedelays, 5)]
        timedelays = timedelays[np.less_equal(timedelays, 5)]

    ws = list()
    centers = list()
    errs_fwhm = list()
    errs_center = list()

    for index, time in enumerate(timedelays):
        if time > 5:
            continue

        linecut = cube[width, 0 : 2 * width, index]
        linecut -= linecut.min()

        try:
            params, pcov = opt.curve_fit(
                peak,
                k,
                linecut,
                p0=(3e3, -1.6e-2, 4, 7),
            )
        except RuntimeError:
            peak_fwhm = 0
        else:
            _, err_center, err_fwhm, _ = np.sqrt(np.diag(pcov))
            peak_cut = peak(k, *params)
            peak_cut -= peak_cut.min()

            interp = interpolate.UnivariateSpline(k, peak_cut - peak_cut.max() / 2)
            r1, r2 = interp.roots()

        ws.append(abs(r1 - r2))
        centers.append((r1 + r2) / 2)
        errs_fwhm.append(err_fwhm)
        errs_center.append(err_center)

    fwhms[refl] = (ws, centers, errs_fwhm, errs_center)

for index, ((refl, (ws, centers, err_fwhm, err_center)), color) in enumerate(
    zip(fwhms.items(), discrete_colors(len(INDICES_DIFFUSE)))
):

    ax_width = axes[index, 0]
    ax_center = axes[index, 1]

    ws = np.asarray(ws)
    ax_width.errorbar(
        timedelays,
        ws - np.mean(ws[timedelays < 0]),
        yerr=err_fwhm,
        linestyle="none",
        markersize=2,
        marker="o",
        color=color,
    )
    ax_width.axhline(y=0, linestyle="dashed", color="k", linewidth=0.5, zorder=np.inf)
    ax_width.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5, zorder=np.inf)

    centers = np.asarray(centers)
    centers -= np.mean(centers[timedelays < 0])
    centers /= np.mean(ws)

    ax_center.errorbar(
        timedelays,
        centers,
        yerr=err_center,
        linestyle="none",
        markersize=2,
        marker="o",
        color=color,
    )
    ax_center.axhline(y=0, linestyle="dashed", color="k", linewidth=0.5, zorder=np.inf)
    ax_center.axvline(x=0, linestyle="dashed", color="k", linewidth=0.5, zorder=np.inf)

    if index != len(INDICES_DIFFUSE) - 1:
        ax_width.xaxis.set_visible(False)
        ax_center.xaxis.set_visible(False)

    ax_width.set_ylabel("$\Delta \sigma$ [$\AA^{-1}$]")
    ax_width.set_ylim([-0.015, 0.015])

    ax_center.set_ylabel(r"$\Delta x_c / \bar{\sigma}$")
    ax_center.yaxis.set_label_position("right")
    ax_center.yaxis.tick_right()
    ax_center.set_ylim([-0.035, 0.035])
    ax_center.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    tag_axis(
        ax_width,
        skued.indices_to_text(*refl),
        y=0.9,
    )

axes[0, 0].set_xlim([-1.5, 5])


axes[-1, 0].set_xlabel("Time-delay [ps]")
axes[-1, 1].set_xlabel("Time-delay [ps]")

plt.subplots_adjust(
    top=0.975, bottom=0.085, left=0.125, right=0.9, hspace=0.2, wspace=0.2
)
