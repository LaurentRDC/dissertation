import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from plotutils import (
    ImageGrid,
    tag_axis,
    draw_hexagon_field,
    draw_hexagon,
    set_height_auto,
)

# CONSTANTS -------------------------------------------------------------------

FIGURE_WIDTH = 6 + 3 / 4  # inches

# Diffraction patterns are rotated by 8 degrees clockwise from aligned
GRAPHITE_ANGLE = 8  # degrees
GRAPHITE_CAMERA_LENGTH = 0.25  # centi-meters
_peak1 = np.array((754, 905))
_peak2 = np.array((1265, 1318))
GRAPHITE_CENTER = np.array(0.5 * (_peak1 + _peak2), dtype=np.int)

# STYLE -----------------------------------------------------------------------

plt.rcParams["figure.figsize"] = (FIGURE_WIDTH, FIGURE_WIDTH)
plt.rcParams["font.size"] = 10
plt.rcParams["savefig.pad_inches"] = 0.0

# -----------------------------------------------------------------------------
