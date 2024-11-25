
# Imports ----------------------------------------------------------------------
import os
import uuid
import subprocess
from typing import List, Union
from Bio import SeqIO

# Run JackHMMER ----------------------------------------------------------------
def run_jackhmmer(
        jackhmmer_path :str, seq_db_path :str,
        fasta_path :str, output_path: str,
        run_remove_first_seq_gaps: bool=True, run_zip_file: bool=False,
        n: int=1, cpu: int=1, incE: Union[None, float]=None,
        los_steps: bool=True,
    ) -> Union[None, str]:
    """
    Run jackHMMER on sequence from 'fasta_path' with DB 'seq_db_path' -> 'output_path' ('.ali' or '.fasta').

    - Automatically converts '.ali' (stockholm format) -> '.fasta' format if 'output_path' is a '.fasta' file.
    - Remove all positions in the MSA that are a gap in the first (target) sequence using perl executable 'a2m2aln' if required.
    - Zip the output file if required.
    - Deletes all intermediate files except the final output
    - Returns final output path (not always same as 'output_path' if file is zipped) or None if execution failed.
    """

    # Guardians
    error_log = f"ERROR in run_jackhmmer('{fasta_path}')"
    if not os.path.isfile(jackhmmer_path):
        print(f"{error_log}: jackhmmer_path='{jackhmmer_path}' does not exists.")
        return None
    if not os.path.isfile(seq_db_path):
        print(f"{error_log}: seq_db_path='{seq_db_path}' does not exists.")
        return None
    if not seq_db_path.endswith(".fasta"):
        print(f"{error_log}: seq_db_path='{seq_db_path}' should end with '.fasta'.")
        return None
    if not os.path.isfile(fasta_path):
        print(f"{error_log}: fasta_path='{fasta_path}' does not exists.")
        return None
    if not fasta_path.endswith(".fasta"):
        print(f"{error_log}: fasta_path='{fasta_path}' should end with '.fasta'.")
        return None
    if not is_nonempty_file(fasta_path):
        print(f"{error_log}: fasta_path='{fasta_path}' is empty.")
        return None
    if not output_path.endswith(".ali") and not output_path.endswith(".fasta"):
        print(f"{error_log}: output_path='{output_path}' should end with '.ali' or '.fasta'.")
        return None
    if not os.path.isdir(os.path.dirname(output_path)):
        print(f"{error_log}: directory of output_path='{output_path}' does not exists.")
        return None
    if os.path.abspath(fasta_path) == os.path.abspath(output_path):
        print(f"{error_log}: fasta_path and output_path ='{fasta_path}' should be different.")
        return None
    if n not in [1, 2, 3, 4, 5]:
        print(f"{error_log}: number of jackHMMER iterations n={n} should be in [1, 2, 3, 4, 5].")
        return None
    if not cpu > 0:
        print(f"{error_log}: number or cpu={cpu} should be > 0.")
        return None
    
    # Determine output format
    run_convert_ali_to_fasta = False
    output_path_ali = output_path
    if output_path.endswith(".fasta"):
        run_convert_ali_to_fasta = True
        output_path_ali = output_path.removesuffix(".fasta") + ".ali"
    
    # Parameters coherence guardian
    if not run_convert_ali_to_fasta and run_remove_first_seq_gaps:
        print(f"{error_log}: impossible to run_remove_first_seq_gaps if not run_convert_ali_to_fasta.")
        return None
    
    # Define jackHMMER cmd
    jackhmmer_cmd = f"{os.path.abspath(jackhmmer_path)} -N {n} -A {os.path.abspath(output_path_ali)} --cpu {cpu}"
    if incE is not None:
        jackhmmer_cmd = f"{jackhmmer_cmd} --incE {incE}"
    jackhmmer_cmd = f"{jackhmmer_cmd} {os.path.abspath(fasta_path)} {os.path.abspath(seq_db_path)}"
    if los_steps:
        print(f"\nRun JackHMMER with: \n$ {jackhmmer_cmd}")

    # CPU excessive usage warning
    if cpu > 3:
        print(f"WARNING in run_jackhmmer('{fasta_path}'): number of CPU(={cpu}) above 3 do not speed-up jackHMMER execution.")

    # DB warning (in case of decoy DB usage)
    if "_cutted" in os.path.basename(seq_db_path):
        print(f"WARNING in run_jackhmmer('{fasta_path}'): seq_db_path='{seq_db_path}' is probably a decoy dataset.")

    # Run jackhmmer
    jackhmmer_status = subprocess.run(
        args=jackhmmer_cmd.split(), check=False,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
    )

    # Detect jackHMMER errors
    if jackhmmer_status.returncode != 0:
        print(f"{error_log}: jackHMMER execution failed.")
        return None
    if not is_nonempty_file(output_path_ali):
        print(f"{error_log}: jackHMMER output file '{output_path_ali}' is empty.")
        return None
    
    # Convert .ali -> .fasta if required
    if run_convert_ali_to_fasta:
        output_path = convert_ali_to_fasta(output_path_ali, output_path, delete_input=True)
        if output_path is None: # Fail the run if convert_ali_to_fasta failed
            print(f"{error_log}: convertion '.ali' -> '.fasta' FAILED.")
            return None
        
    # Remove first sequence gaps if required
    if run_remove_first_seq_gaps:
        output_path = remove_first_seq_gaps(output_path)
        if output_path is None: # Fail the run if remove_first_seq_gaps failed
            print(f"{error_log}: Remove first sequence gaps FAILED.")
            return None
        
    # Zip alignment if required
    if run_zip_file:
        output_path = zip_file(output_path, delete_input=True)
        if output_path is None: # Fail the run if zip_file failed
            print(f"{error_log}: Zip of output file FAILED.")
            return None
        
    # Return
    return output_path


# MSA files processing functions -----------------------------------------------

def convert_ali_to_fasta(input_path: str, output_path: str, delete_input: bool=False) -> Union[None, str]:
    """
    Convert '.ali' (stickholm format) file to a '.fasta' file.
        * Then, deletes input file if required.
        * Returns output_path or None if execution failed.
    Source: https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
    """

    # Guardians
    error_log = f"ERROR in run_jackhmmer -> convert_ali_to_fasta('{input_path}', '{output_path}')"
    if not os.path.isfile(input_path):
        print(f"{error_log}: input file does not exists.")
        return None
    if not is_nonempty_file(input_path):
        print(f"{error_log}: input file is empty.")
        return None
    if not input_path.endswith(".ali"):
        print(f"{error_log}: input file should end with '.ali'.")
        return None
    if not output_path.endswith(".fasta"):
        print(f"{error_log}: output file should end with '.fasta'.")
        return None

    # Run convertion
    try:
        records = SeqIO.parse(input_path, "stockholm")
    except Exception as error:
        print(f"{error_log}: input file parsing failed.")
        print(error)
        return None
    try:
        SeqIO.FastaIO.FastaWriter(output_path, wrap=None).write_file(records)
    except Exception as error:
        print(f"{error_log}: file convertion + writing failed.")
        print(error)
        return None
    
    # Detect errors
    if not is_nonempty_file(output_path):
        print(f"{error_log}: converted output file is empty.")
        return None
    
    # Delete initial '.ali' file if required
    if delete_input:
        if os.path.isfile(input_path):
            os.remove(input_path)

    # Return
    return output_path

def remove_first_seq_gaps(input_path: str) -> Union[None, str]:
    """
    Remove positions that corresponds to a gap in the target (first) sequence in all sequences of a MSA.
        * Stream the read/write IO, so can handle very large files without RAM problems
        * Overwrites input file
        * Returns input_path or None if execution failed
    """

    # Guardians
    error_log = f"ERROR in run_jackhmmer -> remove_first_seq_gaps('{input_path}')"
    if not os.path.isfile(input_path):
        print(f"{error_log}: input file does not exists.")
        return None
    if not is_nonempty_file(input_path):
        print(f"{error_log}: input file is empty.")
        return None
    if not input_path.endswith(".fasta"):
        print(f"{error_log}: input file should end with '.fasta'.")
        return None
    
    # Init temporary file
    tmp_output_path = generate_tmp_path(input_path, "fasta")
    if os.path.isfile(tmp_output_path): os.remove(tmp_output_path)
    
    # Run convertion
    with open(input_path, "r") as fs_input, open(tmp_output_path, "a") as fs_output:

        # Init Fasta Stream
        fasta_stream = FastaStream(fs_input)

        # Process target sequence
        target_seq = fasta_stream.get_next()
        if target_seq is None:
            print(f"{error_log}: no sequences in input file.")
            return None
        if not target_seq.is_valid():
            print(f"{error_log}: Invalid target sequence")
            target_seq.print()
            return None
        non_gap_position = target_seq.get_non_gap_position()
        if sum(non_gap_position) == 0:
            print(f"{error_log}: Invalid target sequence: contains only gaps")
            target_seq.print()
            return None
        target_seq_len = len(target_seq)
        
        # Stream and filter all sequences
        current_seq = target_seq
        while current_seq:
            if target_seq_len != len(current_seq):
                print(f"{error_log}: Length of current sequence [{fasta_stream.current_id}] (l={len(current_seq)}) do not match length of target sequence (l={target_seq_len}):")
                print("[Target Sequence]:")
                target_seq.print()
                print("[Current Sequence]:")
                current_seq.print()
                if os.path.isfile(tmp_output_path): os.remove(tmp_output_path)
                return None
            current_seq.remove_positions(non_gap_position)
            fs_output.write(str(current_seq))
            current_seq = fasta_stream.get_next()

    # Detect errors: empty output file
    if not is_nonempty_file(tmp_output_path):
        print(f"{error_log}: converted output file is empty.")
        if os.path.isfile(tmp_output_path): os.remove(tmp_output_path)
        return None
    
    # Overwrite initial file
    os.rename(tmp_output_path, input_path)

    # Return
    return input_path

def zip_file(input_path: str, delete_input: bool=False) -> Union[None, str]:
    """
    Zip a file and returns .zip file's path or None if execution failed.
        * Then deletes initial file if required.
    """

    # Guardian
    error_log = f"ERROR in run_jackhmmer -> zip_file('{input_path}')"
    if not os.path.isfile(input_path):
        print(f"{error_log}: input file does not exists.")
        return None

    # Init paths
    dir_name = os.path.dirname(input_path)
    filename = os.path.basename(input_path)
    zip_filename = f"{filename}.zip"
    zip_path = os.path.join(dir_name, zip_filename)

    # Zip
    zip_cmd = f"zip {zip_filename} {filename}"
    zip_status = subprocess.run(
        args=zip_cmd.split(), check=False, cwd=dir_name,
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
    )

    # Detect errors
    if zip_status.returncode != 0:
        print(f"{error_log}: zip execution failed.")
        return None
    if not os.path.isfile(zip_path):
        print(f"{error_log}: '.zip' file does not exists.")
        return None
    
    # Delete initial file
    if delete_input:
        if os.path.isfile(input_path):
            os.remove(input_path)

    # Return
    return zip_path



# Dependencies: MSA and FASTA management ---------------------------------------

class FastaSequence():
    """Contained class for a single fasta sequence"""

    def __init__(self, header: str, seq: str):
        self.header = header
        self.seq = seq

    def __len__(self) -> int:
        return len(self.seq)
    
    def __str__(self) -> str:
        return f"{self.header}\n{self.seq}\n"
    
    def is_valid(self) -> bool:
        return self.header.startswith(">") and len(self) > 0

    def get_non_gap_position(self) -> List[bool]:
        return [char != "-" for char in self.seq]

    def remove_positions(self, keep_positions: List[bool]) -> None:
        self.seq = "".join([char for char, keep in zip(self.seq, keep_positions) if keep])
    
    def print(self, max_print_len: int=50) -> None:
        print(f"{self.header} (length={len(self)})\n{self.seq[0:max_print_len]}...")

class FastaStream():
    """Allow to stream fasta sequence one by one from a fasta file"""

    def __init__(self, fs):
        self.fs = fs
        self.current_id = 0
        self.current_line = self._next_line()

    def _next_line(self) -> Union[None, str]:
        line = self.fs.readline()
        if line is None:
            return None
        self.current_id += 1
        return line.removesuffix("\n")
    
    def is_current_line_header(self) -> bool:
        return self.current_line.startswith(">")

    def get_next(self) -> Union[None, FastaSequence]:
        if self.current_line is None:
            return None
        header = self.current_line
        seq_arr = []
        self.current_line = self._next_line()
        while self.current_line:
            if self.is_current_line_header():
                break
            seq_arr.append(self.current_line)
            self.current_line = self._next_line()

        seq = "".join(seq_arr)
        return FastaSequence(header, seq)
    
# Dependencies: FileSystem -----------------------------------------------------

def list_files(input_dir :str, ext :str=Union[None, str]) -> List[str]:
    """List of files in 'input_dir', if 'ext' is not None, only returns '.ext' files."""
    files = [file for file in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, file))]
    if ext is not None:
        files = [file for file in files if file.endswith(f".{ext}")]
    return files

def is_nonempty_file(input_path: str) -> bool:
    """Check if 'input_path' is an existing non-empty file."""
    if not os.path.isfile(input_path):
        return False
    with open(input_path, "r") as fs:
        line = fs.readline()
        return len(line) > 0
    
def find_path(paths_list: List[str], comment: str="") -> str:
    """Find first existing file among a list of candidate paths."""
    for candidate_path in paths_list:
        if os.path.isfile(candidate_path):
            print(f"\nFound candidate path for {comment}:")
            print(f" -> '{candidate_path}'\n")
            return candidate_path
    print(f"No candidate path found for {comment}")
    raise IOError(f"NO EXISTING PATH FOUND for {comment}.")

def generate_tmp_path(template_path: str, ext: str) -> str:
    """Generate random and unique tempory path starting from an initial path"""
    dir_name = os.path.dirname(template_path)
    file_name = os.path.basename(template_path).removesuffix(f".{ext}")
    file_uuid = str(uuid.uuid4()) # Get random and unique id to avoid file clashed for multiple runs
    return os.path.join(dir_name, f"{file_name}_tmp_{file_uuid}.{ext}")