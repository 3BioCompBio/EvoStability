
"""
Small MSA management package.
Author: Matsvei Tsishyn
"""

# Imports ----------------------------------------------------------------------
#import argparse
import os
import subprocess
from typing import Union, List, Tuple, Dict
from math import log, sqrt
from collections import Counter, defaultdict
import numpy as np


# Main -------------------------------------------------------------------------
class MSA():
    """
    Class for MSA objects:
    usage: 
        msa = MSA(msa_path, weights_path, evalues_path)
    example:
        See usage examples in the end of this script
    """

    # Static properties --------------------------------------------------------
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
    GAP = '-'
    NON_STANDARD_AA = "X"
    ACCEPTED_CHARACTERS_SET = set(AMINO_ACIDS + GAP)
    N_STATES = len(ACCEPTED_CHARACTERS_SET)
    DEPTH_CUT_MODES = ["ordered", "random", "uniform"]
    SCORE_MODES = ["evalues", "matrix", "seqid", "seqid_ungapped"]

    # Constructor --------------------------------------------------------------
    def __init__(
            self, msa_path: str, weights_path: Union[None, str]="", evalues_path: Union[None, str]="",
            print_warnings: bool=True, regularisation_factor: float=0.01,
        ):

        # Guadians
        assert msa_path.endswith(".fasta"), f"ERROR in MSA('{msa_path}'): input file should be a '.fasta'."
        assert os.path.exists(msa_path), f"ERROR in MSA('{msa_path}'): MSA file do not exists."
        assert 0.0 < regularisation_factor < 1.0, f"ERROR in MSA('{msa_path}'): impossible regularization factor."

        # Fill basic properties
        self.msa_path = msa_path
        self.weights_path = weights_path
        self.evalues_path = evalues_path
        self.weights_are_set = False
        self.evalues_are_set = False
        self.print_warnings = print_warnings
        self.regularisation_factor = regularisation_factor
        self.file_name = os.path.basename(self.msa_path)
        self.name = self.file_name.removesuffix(".fasta")
        self.sorted_by = "file_order"

        # Read MSA
        self._set_msa(msa_path)
        self._set_weights(weights_path)
        self._set_evalues(evalues_path)
        assert len(self.msa) > 0, f"ERROR in MSA('{msa_path}'): input file contain no sequences."
        assert len(self.msa[0]) > 0, f"ERROR in MSA('{msa_path}'): input file contains sequence of length 0."
        assert MSA.GAP not in self.msa[0], f"ERROR in MSA('{msa_path}'): gaps ('-') detected in template sequence."

        # Init counts
        self.col_count, self.col_frequency, self.col_frequency_reg = None, None, None # counts, frequencies and regularized frequences by colums
        self.col_count_weighted, self.col_frequency_reg_weighted = None, None         # wieghted counts and frequencies
        self.col_aas, self.col_gaps = None, None,                                     # array of counts of aas and gaps by columns
        self.tot_count, self.tot_frequency, self.tot_frequency_reg = None, None, None # counts, frequencies and regularized frequences in total
        self.tot_aas, self.tot_gaps = None, None                                      # total counts of aas and gaps
        self.tot_weights = None
        self._init_counts()

    # Properties ---------------------------------------------------------------
    def __str__(self) -> str:
        return f"MSA('{self.name}')"

    @property
    def length(self) -> int:
        return len(self.msa[0])

    @property
    def depth(self) -> int:
        return len(self.msa)

    @property
    def seq(self) -> str:
        return self.msa[0]
    
    @property
    def ntot(self) -> int:
        return self.length
    
    @property
    def neff(self) -> float:
        return self.tot_weights
    
    def d_group(self) -> str:
        log_neff = np.log10(self.neff)
        n = 0
        while n <= log_neff:
            n += 1
        n -= 1
        if n == 0:
            return "Dnull"
        return f"D{10**n}"
    
    def get_col(self, i: int) -> str:
        return "".join([seq[i] for seq in self.msa])
    
    # Scores -------------------------------------------------------------------

    def CI(self) -> List[float]:
        return [
            sqrt(sum([(col_frequency_reg_i[aa] - self.tot_frequency_reg[aa])**2 for aa in MSA.AMINO_ACIDS]))
            for col_frequency_reg_i in self.col_frequency_reg
        ]

    def LOR(self) -> List[Dict[str, float]]:
        scores = []
        for aa_wt, col_frequency_reg_i in zip(self.seq, self.col_frequency_reg):
            freq_wt = col_frequency_reg_i[aa_wt]
            lor_wt = MSA._get_lor(freq_wt)
            res_scores = {}
            for aa_mt in MSA.AMINO_ACIDS:
                if aa_wt == aa_mt:
                    continue
                freq_mt = col_frequency_reg_i[aa_mt]
                lor_mt = MSA._get_lor(freq_mt)
                res_scores[aa_mt] = lor_wt - lor_mt
            scores.append(res_scores)
        return scores
    
    def LORw(self) -> List[Dict[str, float]]:
        if not self.weights_are_set:
            self._print_warning(f"Access to LORw but weights are not set (all weights are equal to 1.0 by default).")
        scores = []
        for aa_wt, col_frequency_reg_i in zip(self.seq, self.col_frequency_reg_weighted):
            freq_wt = col_frequency_reg_i[aa_wt]
            lor_wt = MSA._get_lor(freq_wt)
            res_scores = {}
            for aa_mt in MSA.AMINO_ACIDS:
                if aa_wt == aa_mt:
                    continue
                freq_mt = col_frequency_reg_i[aa_mt]
                lor_mt = MSA._get_lor(freq_mt)
                res_scores[aa_mt] = lor_wt - lor_mt
            scores.append(res_scores)
        return scores
    
    def LR(self) -> List[Dict[str, float]]:
        scores = []
        for aa_wt, col_frequency_reg_i in zip(self.seq, self.col_frequency_reg):
            freq_wt = col_frequency_reg_i[aa_wt]
            lr_wt = MSA._get_logratio(freq_wt)
            res_scores = {}
            for aa_mt in MSA.AMINO_ACIDS:
                if aa_wt == aa_mt:
                    continue
                freq_mt = col_frequency_reg_i[aa_mt]
                lr_mt = MSA._get_logratio(freq_mt)
                res_scores[aa_mt] = lr_wt - lr_mt
            scores.append(res_scores)
        return scores
    
    def LRw(self) -> List[Dict[str, float]]:
        if not self.weights_are_set:
            self._print_warning(f"Access to LORw but weights are not set (all weights are equal to 1.0 by default).")
        scores = []
        for aa_wt, col_frequency_reg_i in zip(self.seq, self.col_frequency_reg_weighted):
            freq_wt = col_frequency_reg_i[aa_wt]
            lr_wt = MSA._get_logratio(freq_wt)
            res_scores = {}
            for aa_mt in MSA.AMINO_ACIDS:
                if aa_wt == aa_mt:
                    continue
                freq_mt = col_frequency_reg_i[aa_mt]
                lr_mt = MSA._get_logratio(freq_mt)
                res_scores[aa_mt] = lr_wt - lr_mt
            scores.append(res_scores)
        return scores
    
    @staticmethod
    def _get_lor(freq: float) -> float:
        return log(freq / (1.0 - freq))
    
    @staticmethod
    def _get_logratio(freq: float) -> float:
        return log(freq)
    
    def get_scores(self) -> List[Dict]:
        scores_list = []
        depth = self.depth
        n_aa_list = self.col_aas
        ci_list = self.CI()
        lor_list = self.LOR()
        lorw_list = self.LORw()
        lr_list = self.LR()
        lrw_list = self.LRw()
        for i, (aa_wt, n_aa, ci, lor, lorw, lr, lrw, freq) in enumerate(zip(self.seq, n_aa_list, ci_list, lor_list, lorw_list, lr_list, lrw_list, self.col_frequency)):
            for aa_mt in MSA.AMINO_ACIDS:
                if aa_wt == aa_mt:
                    continue
                scores = {
                    "aa_wt": aa_wt,
                    "aa_mt": aa_mt,
                    "freq_wt": freq[aa_wt],
                    "freq_mt": freq[aa_mt],
                    "res_fasta": i+1,
                    "msa_aa": n_aa,
                    "msa_depth": depth,
                    "CI": ci,
                    "LOR": lor[aa_mt],
                    "LORw": lorw[aa_mt],
                    "LR": lr[aa_mt],
                    "LRw": lrw[aa_mt],
                }
                scores_list.append(scores)
        return scores_list

    def write_scores(self, save_path: str, sep: str=",", missing: str="XXX"):
        scores = self.get_scores()
        header = list(scores[0].keys())
        header_str = sep.join(header)
        scores_str = [sep.join([str(e.get(p, missing)) for p in header]) for e in scores]
        output_str = "\n".join([header_str] + scores_str)
        with open(save_path, "w") as f:
            f.write(output_str)
        return self

    # Methods ------------------------------------------------------------------
    def trim_gapped_colums(
            self,
            gap_ratio_thr: float,
            mapping_path: Union[str, None]=None,
            print_logs: bool=False,
        ) -> Dict[int, int]:
        """Trim columns from the MSA which gap_ratio is greater than gap_ratio_thr. Return the fasta-notation ids mappings."""

        # Guardians
        assert 0.0 < gap_ratio_thr < 1.0, f"ERROR in MSA().trim_gapped_colums(): gap_ratio_thr='{gap_ratio_thr}' should be sticktly between 0 and 1."

        # Find ids to trim
        d = self.depth
        ids_to_trim = []
        for i, col_gaps in enumerate(self.col_gaps):
            gap_ratio = col_gaps / d
            if gap_ratio > gap_ratio_thr:
                ids_to_trim.append(i)
        if print_logs:
            print(f"Trim columns of {self} for gap_ratio > {gap_ratio_thr} ({len(ids_to_trim)} trimmed cols): {len(self.seq)} -> {len(self.seq) - len(ids_to_trim)}.")

        # Trim and return mapping
        return self._trim_by_ids(ids_to_trim, mapping_path=mapping_path)

    def _trim_by_ids(self, ids_to_trim: List[int], mapping_path: Union[str, None]=None) -> Dict[int, int]:
        """Trim columns from the MSA by deleting the [ids_to_trim]. Return the fasta-notation ids mappings."""

        # Init
        l = len(self.seq)
        ids_to_trim_set = set(ids_to_trim)

        # Guardians            
        assert len(ids_to_trim) == len(ids_to_trim_set), f"ERROR in MSA()._trim_by_ids(): repetitions are not allowed in ids_to_trim."
        for id in ids_to_trim:
            assert 0 <= id < l, f"ERROR in MSA()._trim_by_ids(): id='{id}' out or range([{0}, {l-1}])."
        
        # Find ids to keep and generate mapping
        ids_mapping = {}
        ids_to_keep = []
        new_id = 0
        for i in range(l):
            if i not in ids_to_trim_set:
                ids_to_keep.append(i)
                ids_mapping[i+1] = new_id+1
                new_id += 1

        # Trim each sequences
        for j, seq in enumerate(self.msa):
            new_seq = []
            for i in ids_to_keep:
                new_seq.append(seq[i])
            self.msa[j] = "".join(new_seq)

        # Reinitialize counts
        self._init_counts()

        # Save mapping
        if mapping_path is not None:
            self.write_trim_mapping(ids_mapping, mapping_path)

        # Return
        return ids_mapping
    
    def write_trim_mapping(self, ids_mapping: Dict[int, int], mapping_path: str) -> str:

        # Guardians
        assert os.path.isdir(os.path.dirname(mapping_path)), f"ERROR in MSA().write_trim_mapping(): directory of mapping_path='{mapping_path}' does not exists."
        assert mapping_path.endswith(".txt"), f"ERROR in MSA().write_trim_mapping(): mapping_path='{mapping_path}' sould be a '.txt' file."

        # Save
        mapping_str = "\n".join([f"{i1} {i2}" for i1, i2 in ids_mapping.items()])
        with open(mapping_path, "w") as fs:
            fs.write(mapping_str)
        return mapping_str
    
    @staticmethod
    def read_trim_mapping(mapping_path: str) -> Dict[int, int]:

        # Guardians
        assert os.path.isfile(mapping_path), f"ERROR in MSA().read_trim_mapping: mapping_path='{mapping_path}' does not exists."

        # Read
        with open(mapping_path, "r") as fs:
            lines = [line.removesuffix("\n").split() for line in fs.readlines()]
        return {int(i1): int(i2) for i1, i2 in lines}

    def cut_depth(
            self, depth: int, mode: str,
            do_weights_update: bool=True,
            print_logs: bool=False,
        ):
        
        # Init and Guardiants
        assert mode in self.DEPTH_CUT_MODES, f"ERROR in MSA.cut_depth({depth}, {mode}): mode should be in {self.DEPTH_CUT_MODES}."
        assert 0 < depth <= self.depth, f"ERROR in MSA.cut_depth({depth}): Depth should be above zero and below current depth."
        d1 = self.depth
        
        # Ordered cut --------------------------------------------
        if mode == "ordered":
            if self.sorted_by == "file_order":
                self._print_warning(f"cut_depth by mode='ordered' but sequences are sorted by 'file_order': think about sorting the MSA first.")
            ids = np.arange(0, len(self.msa))[:depth]

        # Random cut ---------------------------------------------
        elif mode == "random":
            ids = np.arange(1, len(self.msa))
            np.random.shuffle(ids)
            ids = [x for x in ids[:depth-1]]
            ids.append(0) # Always keep first sequence

        # Uniform cut --------------------------------------------
        elif mode == "uniform":
            if self.sorted_by == "file_order":
                self._print_warning(f"cut_depth by mode='uniform' but sequences are sorted by 'file_order': think about sorting the MSA first.")
            ids = np.arange(1, len(self.msa))
            L = len(self.msa) - 1
            ids = [int(i*(L/depth)) for i in range(depth)]
            ids.append(0) # Always keep first sequence

        # Cut
        self._cut_depth_by_ids(ids)

        # Weights Update
        if do_weights_update:
            self._update_weights(print_logs=print_logs)
        else:
            self._print_warning(f"weights are no longer coherent after .cut_depth().")

        # Reinitialize counts
        self._init_counts()

        # Log
        d2 = self.depth
        if print_logs:
            print(f"{self}.cut_depth() by mode='{mode}': {d1} -> {d2}.")

        return self
    
    def cut_threshold(
            self, thr: float, mode: str="matrix", mutation_matrix=None,
            open_gap_score :float= -10.0, extend_gap_score :float= -0.5, tail_gap_score :float=-0.5,
            do_weights_update: bool=True,
            keep_closest_sequences: bool=True,
            print_logs: bool=False,
        ):
        
        # Init and Guardians
        assert mode in self.SCORE_MODES, f"ERROR in MSA.cut_threshold(): mode='{mode}' should be among {self.SCORE_MODES}."
        d1 = self.depth

        # Define compare function
        is_greater = lambda x: x >= thr
        is_less = lambda x: x <= thr
        if mode == "evalues":
            accept_score = is_less if keep_closest_sequences else is_greater
        else:
            accept_score = is_greater if keep_closest_sequences else is_less

        # Get scores
        scores = self.get_msa_pairwise_scores(mode, mutation_matrix, open_gap_score, extend_gap_score ,tail_gap_score)

        # Find ids to keep
        ids = [0] # Always keep the taret sequence
        for id, score in enumerate(scores):
            if accept_score(score):
                ids.append(id + 1) # Shift +1 since scores ignore the target sequence

        # Cut
        self._cut_depth_by_ids(ids)

        # Weights Update
        if do_weights_update:
            self._update_weights(print_logs=print_logs)
        else:
            self._print_warning(f"weights are no longer coherent after .cut_threshold().")

        # Reinitialize counts
        self._init_counts()

        # Log
        d2 = self.depth
        if print_logs:
            inversed_mode_note = "" if keep_closest_sequences else ", INVERTED_MODE (keep farest sequences)"
            mode_str = f"{mode}({mutation_matrix.name})" if mode == "matrix" else mode
            print(f"{self}.cut_threshold() by mode='{mode_str}', thr={thr}{inversed_mode_note}: {d1} -> {d2}.")

        return self

    def write(self, save_path: str, print_logs: bool=False):

        # Check if save_path is a path to a file or a directory
        if os.path.isdir(save_path):
            save_path = os.path.join(save_path, f"{self.name}.fasta")
        else:
            assert save_path.endswith(".fasta"), f"ERROR in MSA.write('{save_path}'): output path should end with a '.fasta'."

        # Log
        if print_logs:
            print(f"{self}.write() to '{save_path}'.")
        
        # Generate MSA fasta file
        msa_fasta = []
        for name, seq in zip(self.seq_names, self.msa):
            msa_fasta.append(f">{name}")
            msa_fasta.append(seq)

        # Write
        with open(save_path, "w") as fs:
            fs.write("\n".join(msa_fasta) + "\n")
        return self
    
    
    # Sequences Scores ---------------------------------------------------------
    def compute_pairwise_score(
            self, seq1: str, seq2: str,
            mode: str="matrix", mutation_matrix=None,
            open_gap_score :float= -10.0, extend_gap_score :float= -0.5, tail_gap_score :float=-0.5,
        ) -> float:

        # Guardians
        MODES = [mode for mode in self.SCORE_MODES if mode != "evalues"]
        assert mode in MODES, f"ERROR in MSA.compute_pairwise_score(): mode='{mode}' should be among {MODES}."
        assert len(seq1) == len(seq2), f"ERROR in MSA.compute_pairwise_score(): seq1 ({len(seq1)}) and seq2 ({len(seq2)}) should have same length."
        assert len(seq1) > 0, f"ERROR in MSA.compute_pairwise_score(): sequences should be of length > 0."

        # Mode: seqid/seqid_ungapped
        if mode == "seqid" or mode == "seqid_ungapped":
            aa_matches_arr = [int(aa1 == aa2) for aa1, aa2 in zip(seq1, seq2) if aa1 != self.GAP and aa2 != self.GAP]
            tot_matches = sum(aa_matches_arr)
            if mode == "seqid":
                return tot_matches / len(seq1)
            else:
                return tot_matches / len(aa_matches_arr)
            
        # Mode: Substitution Matrix
        assert mutation_matrix is not None, f"ERROR in MSA.compute_pairwise_score(): if mode='{mode}' then set a mutation_matrix."

        # Count score
        is_previous_gap = False
        score = 0.0
        for c1, c2 in zip(seq1, seq2):
            if c1 == self.GAP or c2 == self.GAP:
                score += extend_gap_score if is_previous_gap else open_gap_score
                is_previous_gap = True
            else:
                score += mutation_matrix.get(c1, c2)
                is_previous_gap = False

        # Adjust for tail gaps
        if seq1[0] == self.GAP or seq2[0] == self.GAP:
            score = score + tail_gap_score - open_gap_score
        if seq1[-1] == self.GAP or seq2[-1] == self.GAP:
            score = score + tail_gap_score - open_gap_score
        return score / len(seq1)
    
    def get_msa_pairwise_scores(
            self, mode: str="matrix", mutation_matrix=None,
            open_gap_score :float= -10.0, extend_gap_score :float= -0.5, tail_gap_score :float=-0.5,
            include_target_seq: bool=False,
        ) -> List[float]:

        # Guardians
        assert mode in self.SCORE_MODES, f"ERROR in MSA.get_msa_pairwise_scores(): mode='{mode}' should be among {self.SCORE_MODES}."

        # Set if the score of target sequence is included in the array
        start_id = 0 if include_target_seq else 1

        # Mode: e-values
        if mode == "evalues":
            if not self.evalues_are_set:
                self._print_warning(f"Access to evalues but evalues are not set (all evalues are equal to 0.0 by default.)")
            return self.evalues[start_id:]
        
        # Other modes
        return [
            self.compute_pairwise_score(self.seq, query_seq, mode, mutation_matrix, open_gap_score, extend_gap_score, tail_gap_score)
            for query_seq in self.msa[start_id:]
        ]
    
    def get_msa_avg_pairwise_score(
            self, mode: str="matrix", mutation_matrix=None,
            open_gap_score :float= -10.0, extend_gap_score :float= -0.5, tail_gap_score :float=-0.5
        ) -> float:
        scores = self.get_msa_pairwise_scores(mode, mutation_matrix, open_gap_score, extend_gap_score, tail_gap_score)
        return np.mean(scores)
    
    def sort(
            self, mode: str="matrix", mutation_matrix=None,
            open_gap_score :float= -10.0, extend_gap_score :float= -0.5, tail_gap_score :float=-0.5,
            print_logs: bool=False,
        ):

        # Guardians
        assert mode in self.SCORE_MODES, f"ERROR in MSA.sort(): mode='{mode}' should be among {self.SCORE_MODES}."

        # Log
        if print_logs:
            mode_str = f"{mode}({mutation_matrix.name})" if mode == "matrix" else mode
            print(f"{self}.sort() by mode='{mode_str}'.")

        # Sort
        sort_direction = 1.0 if mode == "evalues" else -1.0
        scores = self.get_msa_pairwise_scores(mode, mutation_matrix, open_gap_score, extend_gap_score, tail_gap_score)
        ordered_ids = np.argsort(sort_direction * np.array(scores))   # Get order by best score
        new_msa, new_seq_names, new_weights, new_evalues = [self.msa[0]], [self.seq_names[0]], [self.weights[0]], [self.evalues[0]] # Keep target sequence first
        for i in ordered_ids:
            i_shift = i + 1
            new_msa.append(self.msa[i_shift])
            new_seq_names.append(self.seq_names[i_shift])
            new_weights.append(self.weights[i_shift])
            new_evalues.append(self.evalues[i_shift])
        self.msa = new_msa
        self.seq_names = new_seq_names
        self.weights = new_weights
        self.evalues = new_evalues
        self.sorted_by = mode
        return self

    # Dependencies -------------------------------------------------------------

    # Parsers
    def _set_msa(self, msa_path: str) -> Tuple[List[str], List[str]]:
        msa, seq_names = [], []
        with open(msa_path, "r") as f:
            lines = [line for line in f.readlines() if len(line) > 1]
        for line in lines:
            if line.startswith(">"):
                seq_names.append(self._read_sequence_name(line))
            else:
                msa.append(self._read_sequence(line))
        assert len(msa) == len(seq_names), f"ERROR in MSA('{self.name}'): unequal number of sequences and ids."
        assert len(set([len(seq) for seq in msa])) <= 1, f"ERROR in MSA('{self.name}'): unequal length of sequences in MSA."
        self.msa = msa
        self.seq_names = seq_names
        return msa, seq_names

    def _read_sequence(self, line: str) -> str:
        sequence = []
        unknown_characters = []
        for char in line:
            if char in MSA.ACCEPTED_CHARACTERS_SET:
                sequence.append(char)
            elif char == MSA.NON_STANDARD_AA:
                sequence.append(MSA.GAP)
            elif char == "\n":
                continue
            else:
                sequence.append(MSA.GAP)
                unknown_characters.append(char)
        if len(unknown_characters) > 0:
            self._print_warning(f"Presence of {len(unknown_characters)} unknown characters {set(unknown_characters)} in sequence.")
        return "".join(sequence)

    def _read_sequence_name(self, line: str) -> str:
        return line[1:].replace("\n", "")
    
    def _update_weights(self, tmp_dir: str="/tmp/", print_logs: bool=False):

        # Guardians
        assert os.path.isdir(tmp_dir), f"ERROR in {self}.update_weights(): directory tmp_dir='{tmp_dir}' does not exists."

        # Find plmc executable path
        PLMS_PATHS_LIST = [
            "/Users/pauline/Desktop/Software/plmc-master/bin/plmc", # Pauline Laptop
            "/home/user/Documents/INTERACTOME/publications/DDG_evolution/DDG_evolution_algo/BIN/plmc/bin/plmc", # Matvei Laptop
        ]
        plmc_path = None
        for plmc_path_current in PLMS_PATHS_LIST:
            if os.path.isfile(plmc_path_current):
                plmc_path = plmc_path_current
                break
        if plmc_path is None:
            print("\nPLMC paths list: ")
            for plmc_path_current in PLMS_PATHS_LIST:
                print(f" * '{plmc_path_current}'")
            raise ValueError(f"ERROR in run_plmc: no path found for plmc executable.")
        
        # Init paths
        weights_path = os.path.join(tmp_dir, f"tmp_{self.name}-{self.depth}-weights.txt")
        msa_path = os.path.join(tmp_dir, f"tmp_{self.name}-{self.depth}-msa.fasta")
        
        # Save MSA in temporary fasta file
        self.write(msa_path)

        # Run plmc
        command = f"{plmc_path} --fast --save-weights {weights_path} {msa_path}"
        subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Parse plmc weights
        self._set_weights(weights_path)

        # Remove tmp files
        if os.path.isfile(msa_path): os.remove(msa_path)
        if os.path.isfile(weights_path): os.remove(weights_path)

        # Log
        if print_logs:
            print(f"{self}.update_weights(): weights updates.")

        # Return
        return self

    def _set_weights(self, weights_path: Union[None, str]) -> List[float]:

        # Case when not weights are provided
        if weights_path is None:
            self._print_warning("input weights_path is None: all wieghts are set as 1.0.")
            weights = [1.0 for _ in self.msa]
            self.weights = weights
            return weights

        # Case when weights_path does not exists
        if not os.path.isfile(weights_path):
            raise ValueError(f"ERROR in MSA('{self.name}'): weights_path='{weights_path}' file does not exists: set an existing path or set as None to ignore weights.")

        # Base case
        with open(weights_path, "r") as fs:
            weights = [line for line in fs.readlines()]
        try:
            weights = [float(w) for w in weights]
        except:
            raise ValueError(f"ERROR in MSA('{self.name}'): weights_path='{weights_path}' error in parsing weights file: some lines are not convertible to float.")
        if len(weights) != len(self.msa):
            print(f" * MSA depth: {len(self.msa)}")
            print(f" * Weights length: {len(weights)}")
            raise ValueError(f"ERROR in MSA('{self.name}'): weights_path='{weights_path}': unequal length of weights and sequences.")
        self.weights = weights
        self.weights_are_set = True
        return weights
    
    def _set_evalues(self, evalues_path: Union[None, str]) -> List[float]:

        # Case when not evalues are provided
        if evalues_path is None:
            self._print_warning("input evalues_path is None: all e-values are set as 0.0.")
            evalues = [0.0 for _ in self.msa]
            self.evalues = evalues
            return evalues
        
        # Case when evalues_path does not exists
        if not os.path.isfile(evalues_path):
            raise ValueError(f"ERROR in MSA('{self.name}'): evalues_path='{evalues_path}' file does not exists: set an existing path or set as None to ignore e-values.")

        # Base case
        evalues_map = {}
        with open(evalues_path, "r") as fs:
            lines = [line.split() for line in fs.readlines() if len(line) > 1 and not line.startswith("#")]
            for line in lines:
                seq_id = line[0]
                value = float(line[4])
                evalues_map[seq_id] = value
        evalues_arr = [0.0]
        for seq_name in self.seq_names[1:]:
            seq_name_short = seq_name.split("/")[0]
            evalue = evalues_map[seq_name_short]
            evalues_arr.append(evalue)
        self.evalues = evalues_arr
        self.evalues_are_set = True
        return evalues_arr

    # Counts
    def _init_counts(self):

        # By column counts
        self.col_count = []
        self.col_count_weighted = []
        for i in range(self.length):
            col = self.get_col(i)
            self.col_count.append(Counter(col))
            self.col_count_weighted.append(self._get_weighted_count(col))

        # By col frequencies
        self.col_frequency = [self._to_frequency(counter) for counter in self.col_count]
        self.col_frequency_reg = [self._to_frequency_reg(counter) for counter in self.col_count]
        self.col_frequency_reg_weighted = [self._to_frequency_reg(counter) for counter in self.col_count_weighted]
        self.col_gaps = [count[MSA.GAP] for count in self.col_count]
        self.col_aas = [self.depth - n_gaps for n_gaps in self.col_gaps]
        
        # Total
        self.tot_count = Counter({
            char: sum([counter[char] for counter in self.col_count])
            for char in MSA.ACCEPTED_CHARACTERS_SET
        })
        self.tot_frequency = self._to_frequency(self.tot_count)
        self.tot_frequency_reg = self._to_frequency_reg(self.tot_count)
        self.tot_gaps = sum(self.col_gaps)
        self.tot_aas = sum(self.col_aas)
        self.tot_weights = sum(self.weights)

        return self

    def _to_frequency_reg(self, counter) -> defaultdict:
        total = sum(counter.values())
        if total == 0:
            return defaultdict(float, {char: self._regularize(0.0) for char in MSA.ACCEPTED_CHARACTERS_SET})
        return defaultdict(float, {char: self._regularize(counter[char] / total) for char in MSA.ACCEPTED_CHARACTERS_SET})

    def _to_frequency(self, counter) -> defaultdict:
        total = sum(counter.values())                       # To count gaps in the raw frequencies
        #total = sum(counter.values()) - counter[self.GAP]   # To not count gaps in the raw frequences
        if total == 0:
            return defaultdict(float, {char: 0.0 for char in MSA.ACCEPTED_CHARACTERS_SET})
        return defaultdict(float, {char: counter[char] / total for char in MSA.ACCEPTED_CHARACTERS_SET})
    
    def _regularize(self, freq: float) -> float:
        return (self.regularisation_factor / MSA.N_STATES) + (1.0 - self.regularisation_factor) * freq
    
    def _get_weighted_count(self, aa_array: list) -> Dict[str, float]:
        counter = {aa: 0.0 for aa in self.ACCEPTED_CHARACTERS_SET}
        for aa, weight in zip(aa_array, self.weights):
            counter[aa] += weight
        return counter
    
    def _cut_depth_by_ids(self, ids):
        ids_set = set(ids)
        new_msa, new_seq_names, new_weights, new_evalues = [], [], [], []
        for i, (seq, seq_name, weight, evalue) in enumerate(zip(self.msa, self.seq_names, self.weights, self.evalues)):
            if i in ids_set:
                new_msa.append(seq)
                new_seq_names.append(seq_name)
                new_weights.append(weight)
                new_evalues.append(evalue)
        self.msa = new_msa
        self.seq_names = new_seq_names
        self.weights = new_weights
        self.evalues = new_evalues
        return self

    # Logs
    def _print_warning(self, warning_str: str):
        if self.print_warnings:
            print(f"WARNING in {self}: {warning_str}")
        return self


# Usage Examples ---------------------------------------------------------------

if __name__ == "__main__":

    # Imports ------------------------------------------------------------------
    import matplotlib.pyplot as plt

    # Constants ----------------------------------------------------------------
    PROTEIN = "2l6q_A_2-56"
    MSA_PATH =     f"../input/msa/UNIREF90/JM1/{PROTEIN}.fasta"
    WEIGHT_PATH =  f"../input/plmc/UNIREF90/JM1/{PROTEIN}-weights.txt"
    EVALUES_PATH = f"../input/msa/UNIREF90/JM1/{PROTEIN}.txt"
    
    # Execution ----------------------------------------------------------------
    
    print("\n * Initialize MSA: ")
    # You can omit the weights and the evalues paths if needed by setting them as None
    msa = MSA(MSA_PATH, None, None, print_warnings=True)

    # Or you can include them in the MSA
    msa = MSA(MSA_PATH, WEIGHT_PATH, EVALUES_PATH, print_warnings=True)

    # Access MSA properties
    print("\n * Properties: ")
    print(msa) # MSA name
    print(msa.seq) # Target (first) sequence
    print(f"l={msa.length}, d={msa.depth}") # Sequences-length and MSA-depth
    print("\n * Loop on sequences: ")
    for i in range(5): # Loop on msa sequences
        print()
        print(f"SEQ[{i}]: weight={msa.weights[i]:.3f}, e-value={msa.evalues[i]}")
        print("name:", msa.seq_names[i])
        print(msa.msa[i])
    print("\nSEQ[5]: ...")

    # Get CI and LOR scores
    print("\n * CI and LOR: ")
    cilor = msa.get_scores()
    for i in range(3):
        print(cilor[i])
    print("...")

    # Get scores (for each sequences)
    # NOTE:
    #   by default msa.get_msa_pairwise_scores() omits the score of the target sequence to itself
    #   you can still set [include_target_seq] as True to also get the score to itself so that "len(scores) == msa.depth"
    # NOTE:
    #   in contract to all other modes, mode = "evalues" increases as the sequences are far from the target sequence
    # NOTE:
    #   for mode="matrix", please specify a MutationMatrix and you can also change values of "open_gap_score", "extend_gap_score" and "tail_gap_score"
    print("\n * Scores: ")
    scores_seqid = msa.get_msa_pairwise_scores(mode="seqid") # Sequence-Identity ratio
    scores_ungapped = msa.get_msa_pairwise_scores(mode="seqid_ungapped") # Sequence-Identity ratio but not counting the gaps (only compare aligner parts of the sequences)
    evalues = msa.get_msa_pairwise_scores(mode="evalues")
    
    plt.title(f"MSA scores: {PROTEIN}")
    plt.xlabel("Sequences")
    plt.ylabel("Scores")
    plt.plot(scores_seqid, label="seqid")
    plt.plot(scores_ungapped, label="seqid_ungapped")
    plt.legend()
    plt.show()
    plt.clf()

    plt.title(f"MSA evalues: {PROTEIN}")
    plt.xlabel("Sequences")
    plt.ylabel("e-values")
    plt.plot(evalues, label="e-values")
    plt.legend()
    plt.show()
    plt.clf()

    # Sort MSA: You can sort MSA by the previously defined scores
    # NOTE:
    #   sorting is done such that sequences are sorted from the closest to the farest sequences to target
    #   so sorting order for mode="evalues" is opposed to other soring orders
    print("\n * Sort:")
    msa.sort(mode="seqid", print_logs=True)
    msa.sort(mode="seqid_ungapped", print_logs=True)
    msa.sort(mode="evalues", print_logs=True)

    # Cut MSA by depth: modes in ["ordered", "random", "uniform"]
    # NOTE:
    #   if you want msa.cut_depth(mode="ordered") to have meaning with respect to the score of your choice,
    #   do a msa.sort() first, then the remaining sequences will be the closest to the target sequence according to the score
    # WARNING: weights are not longer coherent after msa.cut_depth()
    # NOTE:
    #   after each cut() of MSA, weights are no longer coherent (since some sequences where deleted).
    #   by default weights are updated with a new plmc run, however this may take some time
    #   if you do not use weights like in LORw and do not want to loose time, you can deactivate it
    #   by setting do_weights_update=False.
    print("\n * Cut MSA by depth: ")
    msa.sort(mode="seqid", print_logs=True) # First sort the MSA by the "seqid" (or another score)
    depth = 80
    msa.cut_depth(depth, mode="ordered", print_logs=True, do_weights_update=False) # Then cut_depth(mode="ordered") to only keep the closest sequences according to "seqid"

    # Cut MSA by threshold
    # NOTE:
    #   keep the closest sequences to target with respect to a threshold and to the score,
    #   thus for "evalues", it will be the smaller scores and for all other scores, it will be the largest scores
    # WARNING: weights are not longer coherent after msa.cut_threshold()
    print("\n * Cut MSA by threshold: ")
    thr = 0.8
    msa.cut_threshold(thr, mode="seqid", print_logs=True) # here weights are updated by default

    # Save MSA (with all its modifications)
    print("\n * Save MSA: ")
    msa.write(f"/tmp/{PROTEIN}.fasta", print_logs=True)
