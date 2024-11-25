
"""
Generate bar plot of correlations of evolutionary models.
"""

# Imports ----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.CSV import CSV
from PARAMETERS import PARAMETERS

# Constants --------------------------------------------------------------------
text_size = 20
cols = [[251/255, 133/255, 0/255], [0/255, 76/255, 146/255], [142/255, 202/255, 230/255]]
bar_width = 1

# Execution --------------------------------------------------------------------

CORRELATORS = ["GEMME_ind", "pycofitness", "LORw", "GEMME_epi", "ArDCA", "LOR", "EVCouplings_epi", "EVCouplings_ind"]
colors = [cols[0], cols[1], cols[0], cols[1], cols[1], cols[0], cols[1], cols[0]]
tsuboyama = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

data = {}
for correlator in CORRELATORS:
    sp = np.mean(tsuboyama.get_col(f"spearman_{correlator}", float))
    data[correlator.replace("_", "-").replace("LORw", "LOR$_w$")] = sp

# Init figure
sns.set_style("whitegrid")
plt.figure(figsize=(10, 8))

# Plot the lines
sns.barplot(x=list(data.keys()), y=list(data.values()), palette=colors, edgecolor='black')

# Axis labels
plt.xlabel("", fontsize=text_size)
plt.ylabel("$\\rho$(Model-$\\Delta \\Delta G$)", fontsize=text_size)

# Axis ticks
pos_xticks = []
for i in range(len(data)):
    pos_xticks.append(i+0.25)
plt.xticks(ticks=pos_xticks, labels=list(data.keys()), rotation=45, fontsize=text_size, ha="right")
plt.yticks(fontsize=text_size)
plt.ylim(0,0.6)

# Add legend
legend_labels = ['Independent models', 'Epistatic models']
legend_handles = [plt.Rectangle((0,0),1,1, edgecolor='black', facecolor=col) for col in cols[0:2]]
#plt.legend(legend_handles, legend_labels, loc='center left', fontsize=text_size, bbox_to_anchor=(1, 0.5))
plt.legend(legend_handles, legend_labels, loc='best', fontsize=text_size)

# Show plot
plt.tight_layout()
plt.show()
