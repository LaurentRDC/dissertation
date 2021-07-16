from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from dissutils import MEDIUM_FIGURE_WIDTH, discrete_colors
from dissutils.snse import SAMPLE_THICKNESS, photocarrier_density

figure, carrier_ax = plt.subplots(1, 1, figsize=(MEDIUM_FIGURE_WIDTH, 2.5))

fluences = np.linspace(0, 14, 1024)
densities = (
    photocarrier_density(fluences, thickness=SAMPLE_THICKNESS) / 1e20
)  # 10^20 / cm^3

carrier_ax.fill_between(
    fluences,
    photocarrier_density(fluences, thickness=SAMPLE_THICKNESS + 5e-9) / 1e20,
    photocarrier_density(fluences, thickness=SAMPLE_THICKNESS - 5e-9) / 1e20,
    facecolor="gray",
    edgecolor="k",
    linestyle="dashed",
    linewidth=1,
    alpha=0.2,
)

carrier_ax.plot(fluences, densities, "-", color=discrete_colors(1)[0])

carrier_ax.set_xlabel("Fluence [mJ/cm$^2$]")
carrier_ax.set_ylabel("$N_{\gamma}$ [$10^{20}$ cm$^{-3}$]")

carrier_ax.grid("on")

carrier_ax.set_xlim([fluences.min(), fluences.max()])
carrier_ax.set_ylim([densities.min(), densities.max()])

plt.tight_layout()
