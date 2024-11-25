
"""
Core-residues vs. Surface-residues: Generate dependencies between RSA and DDG/CI.
"""

# Imports ----------------------------------------------------------------------
import os.path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.CSV import CSV
from src.stats import get_spearman, get_pearson
from src.Bins import Bins
from src.Color import Color, COLORS
from PARAMETERS import PARAMETERS

# Constants --------------------------------------------------------------------
X_PROPERTIES = [
    {"name": "RSA", "range": [0.0, 100.0], "color": COLORS.BLUE_NICE_1.updated_lightness(-0.06).rgb, "scale": 1.0},
]
Y_PROPERTIES = [
    {"name": "DDG", "alias": "$\\Delta \\Delta G$", "steps": 15, "hlines": [0.0], "units": "(kcal/mol)", "scale": 1.0},
    {"name": "CI", "alias": "CI", "steps": 15, "hlines": [0.0], "units": "", "scale": 1.0},
]
PEARSON = "$\\rho$"

# Execution --------------------------------------------------------------------

datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

stats = CSV(["Val1", "Val2", "Pearson", "Spearman"])
for X in X_PROPERTIES:
    for Y in Y_PROPERTIES:

        # Init
        X_NAME, Y_NAME = X["name"], Y["name"]
        if X_NAME == Y_NAME: continue
        X_RANGE = X["range"]
        STEPS = Y["steps"]
        COLOR = X["color"]
        X_DELTA = X_RANGE[1] - X_RANGE[0]
        bins = Bins(X_RANGE[0], X_RANGE[1], STEPS)
        print(f"\n{X_NAME} vs. {Y_NAME} -----------------------------------------------------------")

        # Count
        pearson_arr, spearman_arr = [], []
        for dataset_entry in datasets:
            name = dataset_entry["protein_id"]
            dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")
            dataset = CSV().read(dataset_path)
            y = dataset.get_col(Y_NAME, float, as_numpy=True)
            x = np.minimum(dataset.get_col(X_NAME, float, as_numpy=True), X_RANGE[1])
            x = X["scale"] * x
            y = Y["scale"] * y
            pearson_arr.append(get_pearson(x, y))
            spearman_arr.append(get_spearman(x, y))
            bins.add_measures_list(x, y)
        bins.sort()
        print(bins.length())
        print(bins.total())
        stats.add_entry({
            "Val1": X_NAME, "Val2": Y_NAME,
            "Pearson": np.mean(pearson_arr), "Spearman": np.mean(spearman_arr)
        })
        bins_list = bins.get_bins(50)
        bins_centers = bins.get_center_list(50)

        # Plot
        sns.set_style("whitegrid")
        fig, ax = plt.subplots()

        # Set violins
        parts = ax.violinplot(
            bins_list,
            positions=bins_centers,
            widths= 0.7*X_DELTA / STEPS,
            showmeans=False, showmedians=False, showextrema=False,
            #quantiles = [[0.05, 0.25, 0.75, 0.95] for _ in bins_list],
        )
        for pc in parts['bodies']:
            pc.set_facecolor(COLOR)
            pc.set_edgecolor('black')
            pc.set_linewidth(0.35)
            pc.set_alpha(0.75)
        plt.xlim(X_RANGE[0], X_RANGE[1])

        # Set boxes
        quartile1, medians, quartile3 = bins.get_percentile_list(0.25, 50), bins.get_percentile_list(0.5, 50), bins.get_percentile_list(0.75, 50)
        print(medians)
        percentile_inf, percentile_sup = bins.get_percentile_list(0.05, 50), bins.get_percentile_list(0.95, 50),
        ax.scatter(bins_centers, medians, marker='o', color=(0.9, 0.9, 0.9), s=15, zorder=3)
        ax.vlines(bins_centers, quartile1, quartile3, color=(0.1, 0.1, 0.1), linestyle='-', lw=3)
        ax.vlines(bins_centers, percentile_inf, percentile_sup, color=(0.1, 0.1, 0.1), linestyle='-', lw=1)

        # Set other details
        for hline in Y["hlines"]:
            plt.axhline(y=hline, linestyle="--", alpha=0.7, color="black", linewidth=1.0)
        plt.ylabel(f"Average {Y['alias']} {Y['units']}")
        plt.xlabel(f"{X_NAME}")
        plt.title(f"{Y['alias']} dependency on {X_NAME} ({PEARSON}={np.mean(spearman_arr):.2f})")

        # Save
        plt.show()
        plt.clf()

# Log stats
stats.show()
