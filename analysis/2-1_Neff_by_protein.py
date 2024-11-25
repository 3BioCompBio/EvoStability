
"""
Generate plot between Spearman[DDG, LORw] and log(Neff)
And get average correlations by subset D, D10, D100, ...
"""

# Imports ----------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.CSV import CSV
from src.Color import Color, COLORS
from src.stats import get_spearman
from PARAMETERS import PARAMETERS

# Constants --------------------------------------------------------------------
PROP = "LORw"
PROP_ALIAS = "$LOR_w$"
SUBSETS = ["D10", "D100", "D1000", "D10000"]

# Execution --------------------------------------------------------------------

# Read Tsuboyama
tsuboyama = CSV().read(PARAMETERS.TSUBOYAMA_PATH)
tsuboyama.filter(lambda e: e["subset"] != "Dnull", do_print=True, filter_name="Only group D")
ts_log_neff = np.log10(tsuboyama.get_col("MSA_Neff", float, as_numpy=True))
ts_corr = tsuboyama.get_col(f"spearman_{PROP}", float, as_numpy=True)
print(f" * D: r={np.mean(ts_corr)}")

# Plot Param
sns.set_style("whitegrid")
plt.tick_params(axis='both', which='major', labelsize=8)

# Plot Tsuboyama
data = pd.DataFrame({'X': ts_log_neff, 'Y': ts_corr})
#sns.scatterplot(data, x="X", y="Y")
sns.regplot(
    data=data, x='X', y='Y', scatter=True,
    color=COLORS.BLUE_MILD.rgb,
    scatter_kws={'s': 4, 'alpha': 1.0},
    line_kws={'alpha': 0.7},
    #label="$\\mathcal{D}$",
)

# Labels
plt.xlabel("$N_{eff}$", fontsize=10)
plt.ylabel(f"$\\rho$({PROP_ALIAS}-$\\Delta \\Delta G$)", fontsize=10)
plt.xticks(ticks=[1, 2, 3, 4, 5], labels=["$10^1$", "$10^2$", "$10^3$", "$10^4$", "$10^5$"])
#plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
#plt.legend()

# Save
plt.show()
#plt.savefig(f"./fig/neff_{PROP}.png", dpi=400, bbox_inches='tight')
plt.clf()

# Plor correlations
for D in SUBSETS:
    d_entries = [e for e in tsuboyama if e["subset"] == D]
    ts_log_neff = np.log10([float(e["MSA_Neff"]) for e in d_entries])
    ts_corr = [float(e[f"spearman_{PROP}"]) for e in d_entries]
    print(f" * {D}: r={np.mean(ts_corr)}")