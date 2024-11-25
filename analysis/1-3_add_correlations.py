
"""
Compute all by-protein correlations.
"""

# Imports ---------------------------------------------------------------------
import os.path
from src.CSV import CSV
from src.stats import get_spearman, get_pearson
from PARAMETERS import PARAMETERS

# Execution --------------------------------------------------------------------

# Constants
CORRELATORS = [
    "RSA", "CI",
    "GEMME_ind", "pycofitness", "LORw", "GEMME_epi", "ArDCA", "LOR", "EVCouplings_epi", "EVCouplings_ind",
    "RSA*GEMME_ind", "RSA*pycofitness", "RSA*LORw", "RSA*GEMME_epi", "RSA*ArDCA", "RSA*LOR", "RSA*EVCouplings_epi", "RSA*EVCouplings_ind",
    "PoPMuSiC",
    #"RaSP", "PremPS", "DDMut", "KORPM", "MAESTRO", "DDGun3D",
]

# Read dataset
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Solve MSA
for i, dataset_entry in enumerate(datasets):

    # Init
    name = dataset_entry["protein_id"]
    dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")
    #print(f"\n * Run {i+1} / {len(datasets)}: '{name}': ...")
    dataset = CSV().read(dataset_path)
    DDG = dataset.get_col("DDG", float, as_numpy=True)

    # Assign correlations
    for correlator in CORRELATORS:
        val = dataset.get_col(correlator, float, as_numpy=True)
        sp = get_spearman(val, DDG)
        pr = get_pearson(val, DDG)
        dataset_entry[f"spearman_{correlator}"] = round(sp, 4)
        dataset_entry[f"pearson_{correlator}"] = round(pr, 4)

datasets.write(PARAMETERS.TSUBOYAMA_PATH)