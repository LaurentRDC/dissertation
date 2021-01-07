from skued import spectrum_colors
import numpy as np
import matplotlib.pyplot as plt
from plotutils import FIGURE_WIDTH, tag_axis

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(FIGURE_WIDTH, 4), sharex=True)

x = np.linspace(-0.5, 0.5, 512)
pulse = np.zeros_like(x)
mode = lambda x, n: np.cos(np.pi * (2 * n + 1) * x)

colors = reversed(list(spectrum_colors(4)))
for n, c in enumerate(colors, start=1):
    ax1.plot(x, mode(x, n), color=c)
    pulse += mode(x, n)

for n in range(1, 45):
    pulse += mode(x, n)

ax2.plot(x, pulse / pulse.max(), "-k")

for ax in ax1, ax2:
    ax.xaxis.set_ticks([])

ax1.yaxis.set_ticks([0])
ax2.yaxis.set_ticks([0])

ax2.set_xlabel("Cavity position (a.u.)")
ax1.set_ylabel("E-field amplitude (a.u.)")
ax2.set_ylabel("E-field amplitude (a.u.)")

tag_axis(ax1, "a)", y=0.8)
tag_axis(ax2, "b)", y=0.8)
