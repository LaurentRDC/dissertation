
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import scipy.stats
import scipy.signal as signal
import scipy.interpolate as interpolate
from pathlib import Path
from crystals import Crystal
from plotutils import FIGURE_WIDTH
from plotutils.snse_datasets import overnight4
import skued
from iris import DiffractionDataset

CRYSTAL = Crystal.from_cif(Path("data") / "snse" / "snse_pnma.cif")

# Determine the peak position and the time-series integration bounding box
INDICES_DIFFUSE = [(0, 0, 2), (0, 2, 0), (0, 1, -3), (0, 0, 4)]


figure, ax_widths = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, 0.5 * FIGURE_WIDTH))

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
        cube = cube[:,:,np.less_equal(timedelays, 5)]
        timedelays = timedelays[np.less_equal(timedelays, 5)]

    ws = list()
    errs = list()

    for index, time in enumerate(timedelays):
        if time > 5:
            ws.append(0)
            errs.append(0)
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
            _, _, err_fwhm, _ = np.sqrt(np.diag(pcov))
            peak_cut = peak(k, *params)
            peak_cut -= peak_cut.min()

            interp = interpolate.UnivariateSpline(k, peak_cut - peak_cut.max() / 2)
            r1, r2 = interp.roots()
            peak_fwhm = abs(r1 - r2)

        ws.append(peak_fwhm)
        errs.append(err_fwhm)

    fwhms[refl] = (ws, errs)

for (refl, (ws, errs)), marker, color in zip(fwhms.items(), [".", "D", "^", "*"], skued.spectrum_colors(len(INDICES_DIFFUSE))):
    ax_widths.errorbar(
        timedelays,
        ws,
        yerr=errs,
        linestyle="none",
        marker=marker,
        color=color,
        label=skued.indices_to_text(*refl),
    )

ax_widths.axvline(x=0, linestyle="--", linewidth=0.5, color="k")
ax_widths.legend(
    loc="center",
    ncol=len(INDICES_DIFFUSE),
    bbox_to_anchor=(0.5, 1.1),
    bbox_transform=ax_widths.transAxes,
    framealpha=1,
    edgecolor="none",
)
ax_widths.set_xlim([-1, 5])
ax_widths.set_xlabel("time-delay [ps]")
ax_widths.set_ylabel("FWHM [$\AA^{-1}$]")

plt.tight_layout()
