
"""
Generate MSA files with JackHMMER.
"""


# Imports ----------------------------------------------------------------------
import os.path
from PARAMETERS import PARAMETERS
from src.CSV import CSV
from src.run_jackhmmer import run_jackhmmer

# Execution --------------------------------------------------------------------

# Read dataset
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Loop on protien
for i, entry in enumerate(datasets):

    # Init
    name = entry["protein_id"]
    fasta_path = os.path.join(PARAMETERS.FASTA_DIR, f"{name}.fasta")
    msa_path = os.path.join(PARAMETERS.MSA_DIR, f"{name}.fasta")

    # Skip existing
    print(f" * Run {i+1} / {len(datasets)}: '{name}': ...")
    if os.path.isfile(msa_path):
        print("SKIPPED: ALREADY DONE.")
        continue

    # Run
    run_jackhmmer(
        jackhmmer_path=PARAMETERS.JACKHMMER_PATH,
        seq_db_path=PARAMETERS.UNIREF90_PATH,
        fasta_path=fasta_path,
        output_path=msa_path,
        n=2,
        cpu=2,
        incE=0.0000001,
    )
