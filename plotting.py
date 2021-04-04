# STYLE -----------------------------------------------------------------------
import matplotlib.pyplot as plt
from dissutils import LARGE_FIGURE_WIDTH, FONTSIZE

plt.rcParams["figure.figsize"] = (LARGE_FIGURE_WIDTH, LARGE_FIGURE_WIDTH)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Source Sans 3"
plt.rcParams["font.size"] = FONTSIZE
plt.rcParams["savefig.pad_inches"] = 0.0
plt.rcParams["mpl_toolkits.legacy_colorbar"] = False

# WARNINGS --------------------------------------------------------------------
import warnings
import numpy

warnings.filterwarnings("ignore", category=UserWarning, module="iris*")
warnings.filterwarnings("ignore", category=RuntimeWarning)
numpy.seterr(divide="ignore")
