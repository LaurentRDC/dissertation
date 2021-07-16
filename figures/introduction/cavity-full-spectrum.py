from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from dissutils import LARGE_FIGURE_WIDTH, discrete_colors, tag_axis

DATADIR = Path("data") / "introduction" / "VNA measurements"
COLOR = discrete_colors(1)[0]

# Approximate location of modes which can be used for compression
COMPRESSION_MODES = [3, 8.6, 9, 9.8, 10.2]

freq_full, amp_full, *_ = np.loadtxt(
    fname=DATADIR / "00.csv", dtype=float, delimiter=",", skiprows=8, unpack=True
)
freq_full *= 1e-9  # GHz

# Measurements near resonance are reported as the 2-port matrix (S-matrix?)
freq, real1, imag1, *_ = np.loadtxt(
    fname=DATADIR / "13.csv", dtype=float, delimiter=",", skiprows=8, unpack=True
)
freq *= 1e-9

# Conversion from ratio of voltages to dB https://en.wikipedia.org/wiki/Decibel
to_db = lambda d: 20 * np.log10(d)
amp = to_db(np.abs(real1 + 1j * imag1))

phase = np.arctan2(imag1, real1)
phase += np.pi
phase[phase > np.pi] -= 2 * np.pi


figure = plt.figure(figsize=(LARGE_FIGURE_WIDTH, 4))
gs = gridspec.GridSpec(2, 2, figure=figure)

ax_amp = figure.add_subplot(gs[0, 0])
ax_phase = figure.add_subplot(gs[0, 1], sharex=ax_amp)
ax_amp_full = figure.add_subplot(gs[1, :])

f0 = freq[np.argmax(amp)]
for ax in [ax_amp, ax_phase, ax_amp_full]:
    ax.axvline(x=f0, linestyle="dashed", color="k", linewidth=1)

ax_amp.plot(freq, amp, "-", color=COLOR)
ax_amp.set_xlim([freq.min(), freq.max()])
ax_amp.set_ylabel("Amplitude [dB]")


ax_phase.plot(freq, phase, "-", color=COLOR)
ax_phase.yaxis.tick_right()
ax_phase.yaxis.set_label_position("right")
ax_phase.set_ylabel("Phase [rad]")

ax_amp_full.plot(freq_full, amp_full, "-", color=COLOR)
ax_amp_full.set_ylabel("Amplitude [dB]")
ax_amp_full.set_xlabel("Frequency [GHz]")
ax_amp_full.set_xlim([freq_full.min(), freq_full.max()])

minamp, maxamp = ax_amp_full.get_ylim()
for fc, patch_width in zip(
    COMPRESSION_MODES, [0.3] + [0.1] * (len(COMPRESSION_MODES) - 1)
):
    # Find exact mode location
    masked_amp = np.array(amp_full, copy=True)
    masked_amp[
        np.logical_or(freq_full < fc - 0.1, freq_full > fc + 0.1)
    ] = amp_full.min()
    fc = freq_full[np.argmax(masked_amp)]

    ax_amp_full.add_patch(
        Rectangle(
            xy=(fc - patch_width / 2, minamp),
            width=patch_width,
            height=maxamp - minamp,
            transform=ax_amp_full.transData,
            fc="grey",
            alpha=0.4,
            ec="none",
        )
    )


tag_axis(
    ax_amp,
    "a)",
    y=0.9,
)
tag_axis(ax_phase, "b)", y=0.9, x=0.95, horizontalalignment="right")
tag_axis(ax_amp_full, "c)", x=0.025, y=0.9)

plt.tight_layout()
