
# Imports ----------------------------------------------------------------------
import os
import shutil
import subprocess
from typing import List, Dict

# Main -------------------------------------------------------------------------
class PDBCat():
    """
    PDBCat() container class for Residues of a PDB and assigned RSA and secondary structure values using 'MuSiC -cat'.
    Dependency: MuSiC v4.0-4.1

    Arg:
        pdb_path                                                 Path to PDB file
        output_dir          (="/tmp/")                           Directory to save temporary/output files
        music_path          (=".../MuSiC-4.0.0/music_retro")     Path to MuSiC executalbe (use music_retro and v4.0)
        cache_use           (=False)                             Reuse previously computed results if the output '.cat' file already exists
        cache_only          (=False)                             Force to reuse cache file and do not execute MuSiC
        delete_log_files    (=True)                              True to delete temporary files (logs, .pdb, .in)
        delete_output_files (=False)                             True to delete all generated files and folders (including .cat file)
        chains              (=None)                              Chains to consider (as a string) or None to consider all chains

    Usage
        my_pdb_cat = PDBCat(my_pdb_path)

    Methods:
        my_pdb_cat.get(res_id):     Residue (at "res_id" = chain + position)
        my_pdb_cat.res_arr:         [Residue ...]
        my_pdb_cat.res_map:         {res_id => Residue ...}
        my_pdb_cat.chain_map:       {"chain" => [Residue ...] ...}
    """

    # Constructor --------------------------------------------------------------
    def __init__(
            self,
            pdb_path :str,
            output_dir :str="/tmp/",
            music_path :str="/home/user/Documents/INTERACTOME/softs/MuSiC-4.0/music_retro",
            cache_use :bool=False,
            cache_only :bool=False,
            delete_log_files :bool=True,
            delete_output_files :bool=False,
            chains=None,
        ):

        # Init
        self.pdb_file_name = os.path.basename(pdb_path)
        self.pdb_name = self.pdb_file_name.split(".")[0]
        if chains is not None: self.pdb_name = f"{self.pdb_name}_{chains}"
        self.res_arr: List[Residue] = []
        self.res_map: Dict[str, Residue] = {}
        self.chain_map = {}

        # IO guardian
        if self.pdb_file_name.split(".")[-1] != "pdb" or len(self.pdb_file_name.split(".")) != 2:
            raise ValueError(f"ERROR in PDBCat({pdb_path}): Input file should be a '.pdb' and not contain dots in its name.")
        if not os.path.isdir(output_dir):
            raise ValueError(f"ERROR in PDBCat({pdb_path}): Input 'output_dir'={output_dir} do not exists.")

        # Solve names and paths
        io_dir = os.path.join(output_dir, f"music_cat_{self.pdb_name}")
        path_in_path = os.path.join(io_dir, self.pdb_name+".in")
        tmp_pdb_path = os.path.join(io_dir, self.pdb_name+".pdb")
        cat_path = os.path.join(io_dir, self.pdb_name+".cat")
        log_path = os.path.join(io_dir, "log_"+self.pdb_name+".txt")

        
        # Manage case of cache .cat output files
        if cache_only or (cache_use and os.path.exists(cat_path)):
            pass # do not run 'MuSiC -cat' if cache_only is True or if (cache_use is True and there is an output file)
        else:

            # Create IO dir if needed
            if not os.path.isdir(io_dir): os.mkdir(io_dir)

            # Copy PDB and filter PDB chains
            if chains is None:
                shutil.copyfile(pdb_path, tmp_pdb_path)
            else:
                pdb_str = _filter_pdb(pdb_path, chains)
                with open(tmp_pdb_path, "w") as f:
                    f.write(pdb_str)

            # Generate path.in file
            path_in = "\n".join([
                f"DATA    {os.path.join(os.path.dirname(music_path), 'MuSiC/Data/')}",
                f"PDB     {io_dir}/",
                f"OUTPUT  {io_dir}/",
                f"CAT     {io_dir}/\n"
            ])
            with open(path_in_path, "w") as f:
                f.write(path_in)

            # Run MuSiC cat
            music_last_folder_name = os.path.basename(os.path.dirname(music_path))
            # Adapt run to MuSiC 4.0 or 4.1 assuming version is specified in MuSiC folder name
            sidechain_parameter = "FULLATOM" if "4.1" in music_last_folder_name else ""
            music_cmd = f"{music_path} -cat {self.pdb_name} {sidechain_parameter} {self.pdb_name} -init {path_in_path} -log {self.pdb_name}"
            subprocess.run(music_cmd, shell=True, stdout=open(os.devnull, 'wb')) # Run in silent mode

        # Parse cat file
        res_arr = parse_cat(cat_path)

        # Set residues array
        self.res_arr = res_arr

        # Set residues map (resid -> residue)
        for residue in self.res_arr:
            self.res_map[residue.id] = residue

        # Set chains map
        for res in self.res_arr:
            chain = res.chain
            if chain not in self.chain_map:
                self.chain_map[chain] = []
            self.chain_map[chain].append(res)

        # Set positions
        for chain, chain_res_arr in self.chain_map.items():
            for i, res in enumerate(chain_res_arr):
                res._chain_seq_position = i+1
        for i, res in enumerate(self.res_arr):
            res._seq_position = i+1

        # Remove temporary IO dir and files
        if delete_log_files:
            if os.path.exists(path_in_path): os.remove(path_in_path)
            if os.path.exists(log_path): os.remove(log_path)
            if os.path.exists(tmp_pdb_path): os.remove(tmp_pdb_path)

        if delete_output_files:
            if os.path.exists(path_in_path): os.remove(path_in_path)
            if os.path.exists(log_path): os.remove(log_path)
            if os.path.exists(tmp_pdb_path): os.remove(tmp_pdb_path)
            if os.path.exists(cat_path): os.remove(cat_path)
            try:
                if os.path.isdir(io_dir): os.rmdir(io_dir)
            except:
                print(f"PDBCat({self.pdb_name}) WARNING: Failed to delete cat directory '{io_dir}'.")

    # Methods ------------------------------------------------------------------
    def __len__(self):
        return len(self.res_arr)

    def __contains__(self, res_id):
        return res_id in self.res_map

    def __str__(self):
        return f"<PDBCat '{self.pdb_name}' (chains='{''.join(self.chains)}') with {len(self)} residues>"

    def get(self, res_id: str):
        """Return Residue for a 'res_id' (chain + position)"""
        if res_id not in self.res_map:
            raise ValueError(f"ERROR in PDBCat.get({self.pdb_name}): non-existing position: '{res_id}'.")
        return self.res_map[res_id]

    @property
    def chains(self):
        return list(self.chain_map)

    @property
    def n_chains(self):
        return len(self.chain_map)


def parse_cat(cat_path :str):
    """Return list of Residues containing CAT information"""

    # Default output if no output file
    if not os.path.exists(cat_path):
        print(f"WARNING in parse_cat(): no output file '{cat_path}': return empty output.")
        return []

    # Read file
    with open(cat_path, "r") as f:
        lines = list(f.readlines())

    # Start parsing at RESIDUES lines
    start_i = 0
    while start_i < len(lines)-1 and lines[start_i][0:9] != "#RESIDUES":
        start_i += 1

    # Create array of residues
    res_arr = []
    for i in range(start_i+1, len(lines)):
        line = lines[i]
        if line[0] == "#": break
        res_id = line[0:6].replace(" ", "")
        chain = res_id[0]
        position = res_id[1:]
        aa_name = line[11]
        aa = AminoAcids.get(aa_name)
        rsa = round(float(line[30:40]), 2)
        asa = round(float(line[19:29]), 2)
        sec_str = line[13]
        residue = Residue(chain, position, aa, rsa, asa, sec_str)
        res_arr.append(residue)

    # Return
    return res_arr


# Dependencies -----------------------------------------------------------------

class AminoAcid():

    def __init__(self, id: int, letter: str, name: str, long_name: str):
        self._id = id
        self._letter = letter
        self._name = name
        self._long_name = long_name

    @property
    def id(self):
        return self._id

    @property
    def letter(self):
        return self._letter

    @property
    def name(self):
        return self._name

    @property
    def long_name(self):
        return self._long_name

    def __hash__(self):
        return hash(self.letter)

    def __str__(self):
        return f"Amino Acid: id: {self._id}, letter: {self._letter}, name: {self._name}, long_name: {self._long_name}"


class AminoAcids():

    _missing_value = "X"
    _AAs = [
        AminoAcid( 1, "A", "ALA", "Alanine"),
        AminoAcid( 2, "C", "CYS", "Cysteine"),
        AminoAcid( 3, "D", "ASP", "Aspartate"),
        AminoAcid( 4, "E", "GLU", "Glutamate"),
        AminoAcid( 5, "F", "PHE", "Phenylalanine"),
        AminoAcid( 6, "G", "GLY", "Glycine"),
        AminoAcid( 7, "H", "HIS", "Histidine"),
        AminoAcid( 8, "I", "ILE", "Isoleucine"),
        AminoAcid( 9, "K", "LYS", "Lysine"),
        AminoAcid(10, "L", "LEU", "Leucine"),
        AminoAcid(11, "M", "MET", "Methionine"),
        AminoAcid(12, "N", "ASN", "Asparagine"),
        AminoAcid(13, "P", "PRO", "Proline"),
        AminoAcid(14, "Q", "GLN", "Glutamine"),
        AminoAcid(15, "R", "ARG", "Arginine"),
        AminoAcid(16, "S", "SER", "Serine"),
        AminoAcid(17, "T", "THR", "Thréonine"),
        AminoAcid(18, "V", "VAL", "Valine"),
        AminoAcid(19, "W", "TRP", "Tryptophane"),
        AminoAcid(20, "Y", "TYR", "Tyrosine"),
    ]
    _not_AA = AminoAcid(0, _missing_value, _missing_value, _missing_value)
    _map = {}
    for AA in _AAs:
        _map[AA.letter] = AA
        _map[AA.name] = AA

    @classmethod
    def get(cls, aa_tag: str)-> AminoAcid:
        return cls._map.get(aa_tag.upper(), cls._not_AA)


class Residue():

    _sec_str_grps = {
        "H": "A",
        "G": "A",
        "E": "B",
        "B": "B",
        "T": "X",
        "S": "X",
        "C": "X",
    }
    _sec_str_default_grp = "-"

    def __init__(self, chain: str, position: str, aa: AminoAcid, RSA: float, ASA: float, sec_str: str):
        self._chain = chain
        self._position = position
        self._aa = aa
        self._RSA = RSA
        self._ASA = ASA
        self._sec_str = sec_str
        self._seq_position = None
        self._chain_seq_position = None

    @property
    def chain(self):
        return self._chain

    @property
    def position(self):
        return self._position

    @property
    def aa(self):
        return self._aa

    @property
    def RSA(self):
        return self._RSA

    @property
    def ASA(self):
        return self._ASA

    @property
    def sec_str(self):
        return self._sec_str

    @property
    def sec_str_grp(self):
        return self._sec_str_grps.get(self._sec_str, self._sec_str_default_grp)

    @property
    def id(self):
        return self.chain + self.position

    @property
    def seq_position(self):
        return self._seq_position

    @property
    def chain_seq_position(self):
        return self._chain_seq_position

    @property
    def aa_id(self):
        return self.aa.id

    @property
    def aa_letter(self):
        return self.aa.letter

    @property
    def aa_name(self):
        return self.aa.name

    @property
    def aa_long_name(self):
        return self.aa.long_name

    def __str__(self):
        return f"Residue at {self.id}: {self.aa_name} (RSA={self.RSA}, sec_str={self.sec_str})"


def _filter_pdb(pdb_path: str, chains: str):
    """Parses and filter a PDB structure (only ATOM lines, only first model) by chains.

    Usage: pdb_AB = _filter_pdb("./my_pdb_ABCD.pdb", "AB")

    Args:
        pdb_path ::str          Path to a PDB file
        chains ::str            Chains to keep

    Returns:
        complex ::str                 PDB string with chain in 'chains'
    """

    # Initialize parsed PDB
    pdb_name = os.path.basename(pdb_path)
    complex = []

    # Parse .pdb and filter chains
    model_counter = 0
    cryst1_line_is_missing = True
    chains_in_pdb = set()
    with open(pdb_path, 'r') as f:
        line = f.readline()
        while line:
            prefix = line[0:6]

            # Add ATOM/HETATM lines in correct interactant and in complex
            if prefix == "ATOM  " or prefix == "HETATM":
                chain = line[21]
                chains_in_pdb.add(chain)
                if chain in chains:
                    complex.append(line)

            # Detect presence of CRYST1 line
            elif prefix == "CRYST1":
                cryst1_line_is_missing = False
                complex.append(line)

            # Count model(s) in PDB and break parsing after the 1st one
            elif prefix == "MODEL ":
                model_counter += 1
                if model_counter > 1:
                    print(f"WARNING in PDBCat({pdb_name})._filter_pdb(): PDB contains more than one MODEL:  -> Only MODEL 1 will be considered.")
                    break
            else:
                complex.append(line)
            line = f.readline()

    # Create PDB strings
    complex_str = "".join(complex)

    # Add default CRYST1 line if it is missing
    if cryst1_line_is_missing:
        default_cryst1_line = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n"
        complex_str = default_cryst1_line + complex_str

    # Guardiant for coherence between interacting_chains and PDB
    for chain in chains:
        if chain not in chains_in_pdb:
            raise ValueError(f"ERROR in PDBCat({pdb_name})._filter_pdb(): chain '{chain}' (from input '{chains}') is not in PDB (with chains '{chains_in_pdb}').")

    return complex_str

# Usage ------------------------------------------------------------------------

"""
pdb_path = "/home/user/Documents/INTERACTOME/publications/DDG_evolution/DDG_evolution_algo/input/pdb_raw/1znj.pdb"
pdb = PDBCat(pdb_path, chains="AB")

res = pdb.get("A3")
print(res)

print(pdb.chains)
print(pdb)

print(pdb.chain_map["A"])

rsa_arr = [res.RSA for res in pdb.res_arr]
print(rsa_arr)
"""
