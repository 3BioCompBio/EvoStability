
"""
Generate plot on study about G-proteins.
"""

# Imports ----------------------------------------------------------------------
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.CSV import CSV
from src.stats import get_spearman
from src.Color import Color, COLORS
from PARAMETERS import PARAMETERS

# Constants --------------------------------------------------------------------
mutations_by_pdb_dir = PARAMETERS.DATASETS_DIR
colors_base = [[251/255, 133/255, 0/255], [0/255, 76/255, 146/255]]

# Execution --------------------------------------------------------------------

# Set colors
zw1_dataset = CSV().read(os.path.join(mutations_by_pdb_dir, f"2zw1_A_0-56_V54S.csv"))
gb4_dataset = CSV().read(os.path.join(mutations_by_pdb_dir, f"1gb4_A_1-57_F53D.csv"))
ddg_zw1, lorw_zw1 = zw1_dataset.get_col("DDG", float), zw1_dataset.get_col("LORw", float)
ddg_gb4, lorw_gb4 = gb4_dataset.get_col("DDG", float), gb4_dataset.get_col("LORw", float)

# Plotting ---------------------------------------------------------------------

color1 = Color(colors_base[0][0], colors_base[0][1], colors_base[0][2])
color1 = color1.updated_lightness(+0.00).updated_saturation(-0.15)

color2 = Color(colors_base[1][0], colors_base[1][1], colors_base[1][2])
color2 = color2.updated_lightness(+0.00).updated_saturation(-0.15)

# Scatter plot ZW1
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
data = pd.DataFrame({'X': ddg_zw1, 'Y': lorw_zw1})
sns.regplot(
    data=data, x="X", y="Y",
    color=color1.rgb,
    scatter_kws={'edgecolor': (0.2, 0.2, 0.2), "alpha": 1.0, "s": 45.0},
)
plt.xlabel("$\\Delta \\Delta G$", fontsize=18)   #DDG
plt.ylabel("LOR$_w$", fontsize=18)   #WLOR
plt.xlim(-3.5, 4)
plt.ylim(-6, 12.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.axhline(y=0, color='black', linestyle="--")  # Horizontal line at y=0
plt.axvline(x=0, color='black', linestyle="--")  # Vertical line at x=0
plt.show()

# Scatter plot GB4
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
data = pd.DataFrame({'X': ddg_gb4, 'Y': lorw_gb4})
sns.regplot(
    data=data, x="X", y="Y",
    color=color2.rgb,
    scatter_kws={'edgecolor': (0.2, 0.2, 0.2), "alpha": 1.0, "s": 45.0},
)
plt.xlabel("$\\Delta \\Delta G$", fontsize=18)   #DDG
plt.ylabel("LOR$_w$", fontsize=18)   #WLOR
plt.xlim(-3.5, 4)
plt.ylim(-6, 12.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.axhline(y=0, color='black', linestyle="--")  # Horizontal line at y=0
plt.axvline(x=0, color='black', linestyle="--")  # Vertical line at x=0
plt.show()
