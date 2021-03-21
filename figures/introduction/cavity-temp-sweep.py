import csv
from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import numpy as np

import skued
from plotutils import FIGURE_WIDTH, discrete_colors, tag_axis

DATADIR = Path("data") / "introduction" / "VNA measurements"


def csv_to_kvstore(fname):
    """Parse CSV file into key-value store where keys are
    always filepaths, and values are tuples of
    (bath temperature, cavity temperature)"""
    with open(fname, mode="r") as f:
        reader = csv.reader(f)
        for i in range(3):
            next(reader, None)  # skip rows of headers
        return {Path(row[0]): (float(row[1]), float(row[2])) for row in reader}


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


fig, (ax_amp, ax_center) = plt.subplots(1, 2, figsize=(FIGURE_WIDTH, 3))

temps = list(range(11, 32 + 1))
colors = discrete_colors(len(temps))

# Build a key-value store of the bath and cavity temperatures
# based on filenames
temperatures = csv_to_kvstore(DATADIR / "temperatures.csv")

# We accumulate the resonance across temperatures
# in order to determine the variation in kHz / deg C
# First column is the cavity temperature
# Second column is the resonance frequency
sweep = np.empty(shape=(len(temps), 2), dtype=float)

# Not all traces have the same range
# We accumulate the min and max to restrict
# the plot range later
f_min, f_max = 0, 4
for row, (temp, color) in enumerate(zip(temps, colors)):

    fname = DATADIR / f"{temp}.csv"
    bath_temp, cav_temp = temperatures[fname.relative_to(DATADIR)]
    freq, amp, phase = load_spectrum(fname)
    center = freq[np.argmax(amp)]
    # print(f'temp: {temp}, freq: {center}, maxamp: {np.max(amp)}')

    sweep[row, 0] = cav_temp
    sweep[row, 1] = center

    f_min = max([f_min, center - 0.004])
    f_max = min([f_max, center + 0.004])

    offset = row * (amp.max() - amp.min()) / 5
    ax_amp.plot(
        freq, amp + offset, marker=".", markersize=2, linestyle="none", color=color
    )
    ax_center.errorbar(
        cav_temp, center, marker=".", markersize=10, xerr=0.1, color=color
    )

# Show the center of the lowest and highest temperatures
# We do this last so that the vertical lines appear above the
# traces.
for temp in {min(temps), max(temps)}:
    freq, amp, _ = load_spectrum(DATADIR / f"{temp}.csv")
    center = freq[np.argmax(amp)]
    ax_amp.axvline(x=center, color="k", linestyle="solid", linewidth=1)

ax_amp.text(
    x=0.98,
    y=0.98,
    ha="right",
    va="top",
    s=f"{max(temps)} $^{{\circ}}$C",
    color=colors[-1],
    transform=ax_amp.transAxes,
)

ax_amp.text(
    x=0.98,
    y=0.02,
    ha="right",
    va="bottom",
    s=f"{min(temps)} $^{{\circ}}$C",
    color=colors[0],
    transform=ax_amp.transAxes,
)

# Relationship between resonance and cavity temperature should be linear
(slope, intercept), cov = np.polyfit(
    x=sweep[:, 0], y=sweep[:, 1], deg=1, full=False, cov=True
)
slope_err = np.sqrt(cov[0, 0])
ax_center.plot(
    sweep[:, 0], np.polyval([slope, intercept], sweep[:, 0]), "-k", linewidth=1
)

# For presentation purposes, we report it in kHz/C
slope *= 1e6  # kHz / C
slope_err *= 1e6  # kHz / C

ax_center.text(
    x=0.25,
    y=0.8,
    s=f"(${slope:.1f} \pm {slope_err:.1f} $)" + " kHz / $^{\circ}$C",
    horizontalalignment="left",
    verticalalignment="bottom",
    transform=ax_center.transAxes,
)
ax_center.yaxis.set_major_locator(LinearLocator(numticks=5))

# We hide the y-axis amplitude tick labels because of offsets that makes
# them meaningless.
ax_amp.yaxis.set_ticks([])
ax_amp.set_xlabel("Frequency [GHz]")
ax_amp.set_ylabel("Amplitude [a.u.]")
ax_amp.set_xlim([f_min, f_max])

ax_center.set_xlabel("Cavity temperature [$^{\circ}$C]")
ax_center.set_ylabel("Frequency [GHz]")
ax_center.yaxis.tick_right()
ax_center.yaxis.set_label_position("right")


tag_axis(ax_amp, "a)")
tag_axis(ax_center, "b)", x=0.95, horizontalalignment="right")
plt.tight_layout()
