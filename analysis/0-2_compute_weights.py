
"""
Compute weights with plmc.
"""

# Imports ----------------------------------------------------------------------
import os.path
from PARAMETERS import PARAMETERS
from src.CSV import CSV
from src.run_plmc_weights import run_plmc_weights

# Execution --------------------------------------------------------------------

# Read dataset
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Loop on protien
for i, entry in enumerate(datasets):

    # Init
    name = entry["protein_id"]
    msa_path = os.path.join(PARAMETERS.MSA_DIR, f"{name}.fasta")
    weights_path = os.path.join(PARAMETERS.MSA_DIR, f"{name}-weights.txt")

    # Skip existing
    print(f" * Run {i+1} / {len(datasets)}: '{name}': ...")
    if os.path.isfile(weights_path):
        print("SKIPPED: ALREADY DONE.")
        continue

    # Run
    run_plmc_weights(
        msa_path,
        weights_path,
        reuse_existing_output=True,
        n_cpu=4,
        plmc_path=PARAMETERS.PLMC_PATH
    )