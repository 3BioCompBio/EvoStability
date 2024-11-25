
"""
Log all by-protein Spearman correlations on D and subsets.
"""


# Imports ----------------------------------------------------------------------
import numpy as np
from src.CSV import CSV
from PARAMETERS import PARAMETERS

# Execution --------------------------------------------------------------------

# Constants
CORRELATORS = [
    "RSA", "CI",
    "GEMME_ind", "pycofitness", "LORw", "GEMME_epi", "ArDCA", "LOR", "EVCouplings_epi", "EVCouplings_ind",
    "RSA*GEMME_ind", "RSA*pycofitness", "RSA*LORw", "RSA*GEMME_epi", "RSA*ArDCA", "RSA*LOR", "RSA*EVCouplings_epi", "RSA*EVCouplings_ind",
    "PoPMuSiC",
    "RaSP", "PremPS", "DDMut", "KORPM", "MAESTRO", "DDGun3D",
]
GROUPS = ["D", "D10", "D100", "D1000", "D10000"]

# Read data
tsuboyama = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Compute stats
stats = CSV(["Method"] + GROUPS)
for correlator in CORRELATORS:
    entry = {"Method": correlator}
    for group in GROUPS:
        group_entries = [e for e in tsuboyama if e["subset"] == group or group == "D"]
        sp = np.mean([float(e[f"spearman_{correlator}"]) for e in group_entries])
        entry[group] = sp
    stats.add_entry(entry)

# Log
stats.show(n_entries=100, round_digit=4)

