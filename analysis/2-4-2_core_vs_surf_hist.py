
"""
Core-residues vs. Surface-residues: Generate historgram.
"""

# Imports ----------------------------------------------------------------------
import os.path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.CSV import CSV
from src.Color import Color, COLORS
from PARAMETERS import PARAMETERS

# Constnats --------------------------------------------------------------------
PROPERTY = "LORw"
ALIAS = "LOR$_w$"
RSA_THR = 20.0
COL_CORE = COLORS.ORANGE_NICE_1.updated_lightness(0.08)
COL_SURF = COLORS.BLUE_NICE_1.updated_lightness(-0.06)

# Execution --------------------------------------------------------------------

# Read datasets
entries_core, entries_surf = [], []
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)
for i, dataset_entry in enumerate(datasets):

    # Init
    name = dataset_entry["protein_id"]
    dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")
    dataset = CSV().read(dataset_path)
    for e in dataset:
        if float(e["RSA"]) <= RSA_THR:
            entries_core.append(e)
        else:
            entries_surf.append(e)

# Get values
ddg_core = [float(e["DDG"]) for e in entries_core]
ddg_surf = [float(e["DDG"]) for e in entries_surf]
val_core = [float(e[PROPERTY]) for e in entries_core]
val_surf = [float(e[PROPERTY]) for e in entries_surf]

# Plot Histogram
sns.set_style("whitegrid")
plt.xlabel(ALIAS)
plt.ylabel("Density")
plt.xlim([-10.0, 15.0])
plt.hist(val_surf, 50, density=True, alpha=0.8, color=COL_SURF.rgb, label="Surface residues")
plt.hist(val_core, 50, density=True, alpha=0.6, color=COL_CORE.rgb, label="Core residues")
plt.axvline(0.0              , linestyle="--", color="black",      alpha=0.7, linewidth=1.5)
plt.axvline(np.mean(val_surf), linestyle="--", color=COL_SURF.rgb, alpha=1.0, linewidth=1.5)
plt.axvline(np.mean(val_core), linestyle="--", color=COL_CORE.rgb, alpha=1.0, linewidth=1.5)
plt.legend()
plt.show()