# STYLE -----------------------------------------------------------------------
import matplotlib.pyplot as plt

from dissutils import FONTSIZE, LARGE_FIGURE_WIDTH

plt.rcParams["figure.figsize"] = (LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Source Sans 3"
plt.rcParams["font.monospace"] = "Monoid"
plt.rcParams["font.size"] = FONTSIZE
plt.rcParams["savefig.pad_inches"] = 0.0

# WARNINGS --------------------------------------------------------------------
import warnings

import numpy

warnings.filterwarnings("ignore", category=UserWarning, module="iris*")
warnings.filterwarnings("ignore", category=RuntimeWarning)
numpy.seterr(divide="ignore")

# PLOT ------------------------------------------------------------------------
