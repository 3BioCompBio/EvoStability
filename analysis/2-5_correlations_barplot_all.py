
"""
Generate bar plot of correlations of all models: evol + RSA*evol + DDG-predictors
"""

# Imports ----------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.Color import Color, COLORS
from src.CSV import CSV
from PARAMETERS import PARAMETERS

# Constants --------------------------------------------------------------------
text_size = 18
bar_width = 0.3

# Read datasets
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

independent_models = ["GEMME-ind", "WLOR", "LOR", "EVcouplings-ind", "CI"]
epistatic_models = ["pycofitness", "GEMME-epi", "ArDCA", "EVcouplings-epi"]

data_1 = {
    "LORw":             [np.mean(datasets.get_col("spearman_LORw", float)), np.mean(datasets.get_col("spearman_RSA*LORw", float))],
    "LOR":              [np.mean(datasets.get_col("spearman_LOR", float)), np.mean(datasets.get_col("spearman_RSA*LOR", float))],
    "GEMME-ind":        [np.mean(datasets.get_col("spearman_GEMME_ind", float)), np.mean(datasets.get_col("spearman_RSA*GEMME_ind", float))],
    "EVcouplings-ind":  [np.mean(datasets.get_col("spearman_EVCouplings_ind", float)), np.mean(datasets.get_col("spearman_RSA*EVCouplings_ind", float))],
}

data_2 = {
    "ArDCA":            [np.mean(datasets.get_col("spearman_ArDCA", float)), np.mean(datasets.get_col("spearman_RSA*ArDCA", float))],
    "pycofitness":      [np.mean(datasets.get_col("spearman_pycofitness", float)), np.mean(datasets.get_col("spearman_RSA*pycofitness", float))],
    "EVcouplings-epi":  [np.mean(datasets.get_col("spearman_EVCouplings_epi", float)), np.mean(datasets.get_col("spearman_RSA*EVCouplings_epi", float))],
    "GEMME-epi":        [np.mean(datasets.get_col("spearman_GEMME_epi", float)), np.mean(datasets.get_col("spearman_RSA*GEMME_epi", float))],
}

data_3 = {
    "PoPMuSiC":         [np.mean(datasets.get_col("spearman_PoPMuSiC", float))],
    "RaSP":             [np.mean(datasets.get_col("spearman_RaSP", float))],
    "PremPS":           [np.mean(datasets.get_col("spearman_PremPS", float))],
    "DDMut":            [np.mean(datasets.get_col("spearman_DDMut", float))],
    "KORPM":            [np.mean(datasets.get_col("spearman_KORPM", float))],
    "MAESTRO":          [np.mean(datasets.get_col("spearman_MAESTRO", float))],
    "DDGun3D":          [np.mean(datasets.get_col("spearman_DDGun3D", float))],
}

# Execution --------------------------------------------------------------------

# Convert the data dictionary into a list of tuples
data_list_1 = [(method, values[0], values[1]) for method, values in data_1.items()]
data_list_2 = [(method, values[0], values[1]) for method, values in data_2.items()]
data_list_3 = [(method, values[0]) for method, values in data_3.items()]

# Create a DataFrame from the list of tuples
df1 = pd.DataFrame(data_list_1, columns=['Model', 'Evolutionary models', 'Combination with RSA'])
df2 = pd.DataFrame(data_list_2, columns=['Model', 'Evolutionary models', 'Combination with RSA'])
df3 = pd.DataFrame(data_list_3, columns=['Model', 'ddG predictors'])

COL3 = Color(223, 32, 108).updated_saturation(-0.1)
COL4 = COL3.updated_lightness(0.3)

# Init figure
sns.set_style(style="whitegrid")
plt.figure(figsize=(12, 8))

pos_xticks = []

# Plot the bars
for i, row in df1.iterrows():
    pos_xticks.append(i + 0.25)
    plt.bar(i - 0.2, row['Evolutionary models'], width=bar_width, color=COLORS.ORANGE_NICE_1.rgb, label='Independent', edgecolor='black')    
    plt.bar(i + 0.2, row['Combination with RSA'], width=bar_width, color=COLORS.BLUE_NICE_1.rgb, label='RSA ⨀ Independent', edgecolor='black')

for i, row in df2.iterrows():
    pos_xticks.append(len(df1) + 0.2 + i + 0.25)
    plt.bar(len(df1) + 0.2 + i - 0.2, row['Evolutionary models'], width=bar_width, color=COLORS.ORANGE_NICE_1.rgb, label='Epistatic ', edgecolor='black')    
    plt.bar(len(df1) + 0.2 + i + 0.2, row['Combination with RSA'], width=bar_width, color=COLORS.BLUE_NICE_1.rgb, label='RSA ⨀ Epistatic', edgecolor='black')

for i, row in df3.iterrows():
    pos_xticks.append(len(df1) + len(df2) + 0.2 + i - 0.4*i + 0.25)
    plt.bar(len(df1) + len(df2) + 0.2 + i - 0.4*i, row['ddG predictors'], width=bar_width, color=COL3.rgb, label='$\\Delta \\Delta G$ predictors', edgecolor='black')

# Dotted black vertical line to separate the plots
plt.axvline(x=len(df1) - 0.4, color='black', linestyle='--')
plt.axvline(x=len(df1) + len(df2) - 0.2, color='black', linestyle='--')

# Axis labels
plt.xlabel("", fontsize=text_size)
plt.ylabel("$\\rho$(Model-$\\Delta \\Delta G$)", fontsize=text_size)

# Axis ticks
#pos_xticks = []
#for i in range(len(df1) + len(df2) + len(df3)):
#    pos_xticks.append(i+0.25)
plt.xticks(ticks=pos_xticks, labels=[x.replace("LORw", "LOR$_w$") for x in list(df1['Model'])] + [x.replace("LORw", "LOR$_w$") for x in list(df2['Model'])] + list(df3['Model']), rotation=45, fontsize=14, ha="right")
plt.yticks(fontsize=text_size)
plt.ylim(0,0.85)

plt.text(0.4, 0.65, "Independent", fontsize=18)
plt.text(5.0, 0.65, "Epistatic", fontsize=18)
plt.text(9.0, 0.65, "$\\Delta \\Delta G$ predictors", fontsize=18)


# Add legend
legend_labels = ['Evol', 'RSA $\\odot$ Evol', "$ \\Delta \\Delta G$ predictors"]
legend_handles = [plt.Rectangle((0,0),1,1, edgecolor='black', facecolor=col) for col in [COLORS.ORANGE_NICE_1.rgb, COLORS.BLUE_NICE_1.rgb, COL3.rgb]]
plt.legend(legend_handles, legend_labels, loc='upper left', fontsize=14)

# Remove vertical gridlines
plt.grid(axis='x')

# Show plot
plt.tight_layout()
plt.show()
