# STYLE -----------------------------------------------------------------------
import matplotlib.pyplot as plt
from plotutils import FIGURE_WIDTH

plt.rcParams["figure.figsize"] = (FIGURE_WIDTH, FIGURE_WIDTH)
plt.rcParams["font.size"] = 10
plt.rcParams["savefig.pad_inches"] = 0.0

# WARNINGS --------------------------------------------------------------------
import warnings

warnings.filterwarnings("ignore", category=UserWarning, module="iris*")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="npstreams*")
