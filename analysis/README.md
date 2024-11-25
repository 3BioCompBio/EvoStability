
# Exploring evolution to uncover insights into protein mutational stability

Some of the statistical analysis from the paper [1] are provided here.  
Link: <https://www.biorxiv.org/content/10.1101/2024.05.28.596203v2>  
[1] Hermans P., Tsishyn M., Schwersensky M, Rooman, M. & Pucci, F. (2024). Exploring evolution to uncover insights into protein mutational stability.  

## Usage

- Install python version 3.11 or later.
- Install pip dependencies (`numpy`, `scipy`, `biopython`, `matplotlib` and `seaborn`).
- In `./PARAMETERS.py`, set the PATHS to all the required softwares and databases (described in the Method section of the main paper).
- From the root directory `./analysis/`, run all python scripts in the numbered order.


- **NOTE**: All MSA and wieghts files are pre-computed and zipped in the `../msa/` folder, so you can just unzip them and skip execution of `0-*.py` scripts.
- **NOTE**: RSA values, secondary structure and MSA-base properties like CI, LOR and LORw are already pre-computed in the dataset files, so you can skip execution of `1-*.py` scripts.
- **NOTE**: Values from MSA-based models GEMME, EVCouplings and pycofitness are pre-computed and assigned to the dataset files. Runs of these programs can be replicated using their respective source code.
- **NOTE**: Values of PoPMuSiC are pre-computed and assigned to the `./datasets/*.csv` files. Run of PoPMuSiC can be replicated using our Web Server (free for academic uses) (https://soft.dezyme.com/login).
- **NOTE**: To avoid eventual troubles with Copyright, we do not provide all computed values of structure-based DDG predictors but only their correlation with the DDG for each studied protein.
