
"""
Compute RSA and secodary structure with MuSiC.
"""

# Imports ---------------------------------------------------------------------
import os.path
from sys import exit
from typing import Dict
from PARAMETERS import PARAMETERS
from src.PDBCat import PDBCat
from src.CSV import CSV

# Execution -------------------------------------------------------------------

# MuSiC WARNING
MUSIC_WARNING = f"""
WARNING: Missing MuSiC executable file: '{PARAMETERS.MUSIC_PATH}'.

RSA and secodary structure are computed using our in hous Software MuSiC which source code is not provided here.
If you do not have access to a MuSiC executable, you can:

    - RSA and secodary structure are already assigned in the provided datasets, so MuSiC runs are not required to replicate current analysis
    - You can run MuSiC using our PoPMuSiC Web Server (free for academic uses): https://soft.dezyme.com/login
    - You can use another software to compute RSA like DSSP: https://swift.cmbi.umcn.nl/gv/dssp/
"""
if not os.path.isfile(PARAMETERS.MUSIC_PATH):
    print(MUSIC_WARNING)
    exit()

# Read dataset
datasets = CSV().read(PARAMETERS.TSUBOYAMA_PATH)

# Loop on datasets
pdb_map: Dict[str, PDBCat] = {}
for i, dataset_entry in enumerate(datasets):

    # Init
    name = dataset_entry["protein_id"]
    print(f" * Solve RSA for PDB {i+1} / {len(datasets)}: '{name}'")
    pdb_path = os.path.join(PARAMETERS.PDB_DIR, f"{name}.pdb")
    dataset_path = os.path.join(PARAMETERS.DATASETS_DIR, f"{name}.csv")

    # Solve RSA and secodary structure
    pdb = PDBCat(pdb_path, "/tmp/", PARAMETERS.MUSIC_PATH, delete_log_files=True, delete_output_files=True)

    # Fill RSA and secodary structure values
    dataset = CSV().read(dataset_path)
    dataset.add_empty_col("RSA", allow_replacement=True)
    dataset.add_empty_col("secondary_structure", allow_replacement=True)
    for mutation_entry in dataset:
        mut = mutation_entry["mutation_pdb"]
        resid = mut[1:-1]
        res = pdb.res_map[resid]
        mutation_entry["RSA"] = str(round(res.RSA, 2))
        mutation_entry["secondary_structure"] = res.sec_str_grp
    dataset.write(dataset_path)