# STYLE -----------------------------------------------------------------------
import matplotlib.pyplot as plt
from plotutils import FIGURE_WIDTH, FONTSIZE

plt.rcParams["figure.figsize"] = (FIGURE_WIDTH, FIGURE_WIDTH)
plt.rcParams["font.size"] = FONTSIZE
plt.rcParams["savefig.pad_inches"] = 0.0
plt.rcParams["mpl_toolkits.legacy_colorbar"] = False

# WARNINGS --------------------------------------------------------------------
import warnings
import numpy

warnings.filterwarnings("ignore", category=UserWarning, module="iris*")
warnings.filterwarnings("ignore", category=RuntimeWarning)
numpy.seterr(divide="ignore")
