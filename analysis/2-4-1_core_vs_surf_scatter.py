
"""
Core-residues vs. Surface-residues: Generate scatter plots of correlations.
"""

# Imports ----------------------------------------------------------------------
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.stats import get_spearman
from src.CSV import CSV
from src.Color import Color, COLORS
from PARAMETERS import PARAMETERS

# Constnats --------------------------------------------------------------------
PROPERTY = "LORw"
ALIAS = "LOR$_w$"
RSA_THR = 20.0
COL_CORE = COLORS.ORANGE_NICE_1.updated_lightness(0.10)
COL_SURF = COLORS.BLUE_NICE_1.updated_lightness(-0.06)

# Execution --------------------------------------------------------------------

# Read datasets
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Get by protein correlations
neff = []
sp_core = []
sp_surf = []
for i, dataset_entry in enumerate(datasets):

    # Init
    name = dataset_entry["protein_id"]
    dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")
    dataset = CSV().read(dataset_path)
    entries_core = [e for e in dataset if float(e["RSA"]) <= RSA_THR]
    entries_surf = [e for e in dataset if float(e["RSA"]) > RSA_THR]
    ddg_core = [float(e["DDG"]) for e in entries_core]
    ddg_surf = [float(e["DDG"]) for e in entries_surf]
    val_core = [float(e[PROPERTY]) for e in entries_core]
    val_surf = [float(e[PROPERTY]) for e in entries_surf]
    sp_core.append(get_spearman(ddg_core, val_core))
    sp_surf.append(get_spearman(ddg_surf, val_surf))
    neff.append(float(dataset_entry["MSA_Neff"]))
log_neff = np.log10(np.array(neff))

# Plot
sns.set_style("whitegrid")
data_core = pd.DataFrame({'X': log_neff, 'Y': sp_core})
data_surf = pd.DataFrame({'X': log_neff, 'Y': sp_surf})
#sns.scatterplot(data, x="X", y="Y")
sns.regplot(
    data=data_surf, x='X', y='Y', scatter=True,
    color=COL_SURF.rgb,
    scatter_kws={'s': 6, 'alpha': 0.9},
    line_kws={'alpha': 0.8},
    label="Surface residues",
)
sns.regplot(
    data=data_core, x='X', y='Y', scatter=True,
    color=COL_CORE.rgb,
    scatter_kws={'s': 6, 'alpha': 0.9},
    line_kws={'alpha': 0.8},
    label="Core residues",
)
plt.axhline(0.0, linestyle="--", color="grey", alpha=0.7, linewidth=1.3)
plt.xlabel("$N_{eff}$")
plt.ylabel(f"$\\rho$({ALIAS}-$\\Delta \\Delta G$)")
plt.xlim([1.0, 5.0])
plt.xticks(ticks=[1, 2, 3, 4, 5], labels=["$10^1$", "$10^2$", "$10^3$", "$10^4$", "$10^5$"])
plt.legend()
plt.show()
plt.clf()