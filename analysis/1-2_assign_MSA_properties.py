
"""
Compute MSA-based properties like CI, LOR and LORw.
"""

# Imports ---------------------------------------------------------------------
import os.path
from src.CSV import CSV
from src.MSA import MSA
from PARAMETERS import PARAMETERS

# Execution --------------------------------------------------------------------

# Read dataset
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)
datasets.add_empty_col("sequence_length",  allow_replacement=True)
datasets.add_empty_col("MSA_Ntot",  allow_replacement=True)
datasets.add_empty_col("MSA_Neff",  allow_replacement=True)
datasets.show()

# Solve MSA
for i, dataset_entry in enumerate(datasets):

    # Init
    name = dataset_entry["protein_id"]
    msa_path = os.path.join(PARAMETERS.MSA_DIR, f"{name}.fasta")
    weights_path = os.path.join(PARAMETERS.MSA_DIR, f"{name}-weights.txt")
    dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")
    print(f" * Run {i+1} / {len(datasets)}: '{name}': ...")

    # Run MSA
    msa = MSA(msa_path, weights_path, None, print_warnings=False)
    scores = msa.get_scores()

    # Format scores
    scores = {e["aa_wt"]+str(e["res_fasta"])+e["aa_mt"]: e for e in scores}
    dataset_entry["sequence_length"] = msa.length
    dataset_entry["MSA_Ntot"] = msa.depth
    dataset_entry["MSA_Neff"] = str(round(msa.neff, 1))

    # Fill MSA-based properties
    dataset = CSV().read(dataset_path)
    dataset.add_empty_col("msa_depth_local", allow_replacement=True)
    dataset.add_empty_col("msa_depth_total", allow_replacement=True)
    dataset.add_empty_col("CI", allow_replacement=True)
    dataset.add_empty_col("LOR", allow_replacement=True)
    dataset.add_empty_col("LORw", allow_replacement=True)
    dataset.add_empty_col("RSA*LOR", allow_replacement=True)
    dataset.add_empty_col("RSA*LORw", allow_replacement=True)
    for mutation_entry in dataset:
        mut = mutation_entry["mutation_fasta"]
        scores_mut = scores[mut]
        mutation_entry["msa_depth_local"] = scores_mut["msa_aa"]
        mutation_entry["msa_depth_total"] = scores_mut["msa_depth"]
        mutation_entry["CI"] = round(scores_mut["CI"], 4)
        LOR = scores_mut["LOR"]
        LORw = scores_mut["LORw"]
        RSA = float(mutation_entry["RSA"])
        mutation_entry["LOR"] = round(LOR, 4)
        mutation_entry["LORw"] = round(LORw, 4)
        mutation_entry["RSA*LOR"] = round((1.0 - min(RSA, 100.0)/100.0) * LOR, 4)
        mutation_entry["RSA*LORw"] = round((1.0 - min(RSA, 100.0)/100.0) * LORw, 4)
    dataset.write(dataset_path)

datasets.write(PARAMETERS.TSUBOYAMA_PATH)