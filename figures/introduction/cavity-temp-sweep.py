import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LinearLocator, StrMethodFormatter

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis

DATADIR = Path("data") / "introduction" / "VNA measurements"
MAXPOINTS = 2048  # Maximum number of points for each amplitude trace


def kelvins(C):
    return C + 273


def load_spectrum(fname):
    """
    Read a file from this directory

    Parameters
    ----------
    fname : Path
        File name

    Returns
    -------
    freq : ndarray, shape (N,)
        Frequency [GHz]
    amplitude : ndarray, shape (N,)
        Amplitude [dB]
    phase : ndarray, shape(N,)
        Phase [Radians]
    """
    freq, real1, imag1, *_ = np.loadtxt(
        fname=fname, dtype=float, delimiter=",", skiprows=8, unpack=True
    )

    freq *= 1e-9

    # Conversion from ratio of voltages to dB
    # https://en.wikipedia.org/wiki/Decibel
    to_db = lambda d: 20 * np.log10(d)

    amplitude = to_db(np.abs(real1 + 1j * imag1))
    phase = np.arctan2(imag1, real1)

    return freq, amplitude, phase


fig, (ax_amp, ax_center) = plt.subplots(1, 2, figsize=(LARGE_FIGURE_WIDTH, 2.75))

temps = list(range(11, 32 + 1))
colors = discrete_colors(len(temps))

# Build a key-value store of the bath and cavity temperatures
# based on filenames
with open(DATADIR / "temperatures.csv", mode="r") as f:
    reader = csv.reader(f)
    for i in range(3):
        next(reader, None)  # skip rows of headers
    temperatures = {Path(row[0]): (float(row[1]), float(row[2])) for row in reader}

# We accumulate the resonance across temperatures
# in order to determine the variation in kHz / deg C
# First column is the cavity temperature
# Second column is the resonance frequency
sweep = np.empty(shape=(len(temps), 2), dtype=float)

# Not all traces have the same range
# We accumulate the min and max to restrict
# the plot range later
f_min, f_max = 0, 4

# Reversed list so that the curves at lower temperatures appear
# above the ones with lower temperatures (i.e. plotted later)
for row, (temp, color) in reversed(list(enumerate(zip(temps, colors)))):

    fname = DATADIR / f"{temp}.csv"
    bath_temp, cav_temp = temperatures[fname.relative_to(DATADIR)]
    freq, amp, phase = load_spectrum(fname)
    center = freq[np.argmax(amp)]

    sweep[row, 0] = cav_temp
    sweep[row, 1] = center

    f_min = max([f_min, center - 0.004])
    f_max = min([f_max, center + 0.004])

    downsampling = len(freq) // MAXPOINTS  # Traces have between 16k and 32k points!
    ax_amp.plot(
        freq[::downsampling],
        amp[::downsampling],
        color=color,
    )
    ax_center.errorbar(
        kelvins(cav_temp), center, marker=".", markersize=10, xerr=0.1, color=color
    )

cavity_temps = sweep[:, 0]

# Show the center of the lowest and highest temperatures
# We do this last so that the vertical lines appear above the
# traces.
for temp, zorder in zip([min(temps), max(temps)], [np.inf, 0]):
    freq, amp, _ = load_spectrum(DATADIR / f"{temp}.csv")
    center = freq[np.argmax(amp)]
    ax_amp.axvline(
        x=center, color="k", linestyle="dashed", linewidth=0.5, zorder=zorder
    )

ax_amp.text(
    x=0.7,
    y=0.7,
    ha="left",
    va="bottom",
    s=f"{kelvins(min(cavity_temps)):.1f} K",
    color=colors[0],
    transform=ax_amp.transAxes,
)

ax_amp.text(
    x=0.3,
    y=0.7,
    ha="right",
    va="bottom",
    s=f"{kelvins(max(cavity_temps)):.1f} K",
    color=colors[-1],
    transform=ax_amp.transAxes,
)

# Relationship between resonance and cavity temperature should be linear
(slope, intercept), cov = np.polyfit(
    x=kelvins(sweep[:, 0]), y=sweep[:, 1], deg=1, full=False, cov=True
)
slope_err = np.sqrt(cov[0, 0])
Ts = np.linspace(284, 305, num=32)
ax_center.plot(
    Ts, np.polyval([slope, intercept], Ts), "-k", linewidth=1
)

ax_amp.set_xlabel("Frequency [GHz]")
ax_amp.set_ylabel("Amplitude [dB]")
ax_amp.set_xlim([f_min, f_max])
ax_amp.set_ylim([-57, -34])

ax_center.set_xlim([Ts.min(), Ts.max()])
ax_center.set_xticks([284, 288, 292, 296, 300, 304])
ax_center.yaxis.set_major_locator(LinearLocator(numticks=5))
ax_center.set_xlabel("Cavity temperature [K]")
ax_center.set_ylabel("Frequency [GHz]")
ax_center.yaxis.tick_right()
ax_center.yaxis.set_label_position("right")
ax_center.yaxis.set_major_formatter(StrMethodFormatter("{x:.4f}"))

tag_axis(ax_amp, "a)", y=0.93)
tag_axis(ax_center, "b)", x=0.95, y=0.93, horizontalalignment="right")
plt.tight_layout()
