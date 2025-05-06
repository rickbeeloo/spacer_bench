import os
import random
import subprocess
import tempfile
import polars as pl
import json
import re
import pysam
import pysam.samtools
import glob
from tqdm import tqdm
import parasail as ps
import pyfastx as pfx
from needletail import parse_fastx_file #, NeedletailError, normalize_seq
# from needletail import reverse_complement as needletail_reverse_complement
import hashlib
# import time

pl.Config.set_tbl_cols(-1)


def read_fasta_needletail(fasta_file):
    seqs=[]
    seq_ids=[]
    for record in parse_fastx_file(fasta_file):
        seqs.append(record.seq)
        seq_ids.append(record.id)
    return seq_ids, seqs

def read_fasta_polars(fasta_file, idcol="contig_id", seqcol="contig_seq"):
    # read fasta file with needletail
    seq_ids, seqs = read_fasta_needletail(fasta_file)
    return pl.DataFrame({idcol: seq_ids, seqcol: seqs})

def create_fasta_index_needletail(fasta_file, fai_output=None):
    """Create a FASTA index (.fai) file similar to samtools faidx. Currently BORKED.
    
    Format: NAME LENGTH OFFSET LINEBASES LINEWIDTH
    - NAME: Name of sequence
    - LENGTH: Length of sequence
    - OFFSET: Offset of first base in sequence (including header)
    - LINEBASES: Number of bases per line
    - LINEWIDTH: Number of bytes per line (including newline)
    
    Args:
        fasta_file (str): Path to input FASTA file
        fai_output (str, optional): Path to output .fai file. If None, uses fasta_file + '.fai'
    
    Returns:
        str: Path to created index file
    """

    if fai_output is None:
        fai_output = fasta_file + '.fai'
        
    index_entries = []
    offset = 0
    
    with open(fasta_file, 'rb') as f:
        for record in parse_fastx_file(fasta_file):
            # Get sequence name
            name = record.id
            
            # Get sequence length
            seq_len = record.seq.__sizeof__()+1
            
            # Calculate line metrics from first sequence line
            seq_lines = record.seq.split('\n')
            if len(seq_lines) > 0:
                linebases = len(seq_lines[0])
                linewidth = linebases + 1  # Add 1 for newline
            else:
                linebases = seq_len
                linewidth = seq_len + 1
                
            # Calculate offset of first base
            header_line = f">{name}\n"
            seq_offset = offset + len(header_line)
            
            # Add entry
            index_entries.append(f"{name}\t{seq_len}\t{seq_offset}\t{linebases}\t{linewidth}")
            
            # Update offset for next sequence
            # Account for header, sequence, and newlines
            offset = f.tell()
            
    # Write index file
    with open(fai_output, 'w') as f:
        f.write('\n'.join(index_entries))
        
    return fai_output

def generate_random_sequence(length):
    return ''.join(random.choice('ATCG') for _ in range(length))

def apply_mismatches(sequence, n_mismatches):
    mismatch_positions = random.sample(range(len(sequence)), n_mismatches)
    sequence_list = list(sequence)
    for pos in mismatch_positions:
        original_base = sequence_list[pos]
        new_base = random.choice([b for b in 'ATCG' if b != original_base])
        sequence_list[pos] = new_base
    return ''.join(sequence_list)

def reverse_complement(sequence):
    """reverse complement non+ambiguous bases."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'U': 'A', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join(complement[base] for base in reversed(sequence))

def simulate_data(contig_length_range, spacer_length_range, n_mismatch_range, sample_size_contigs, sample_size_spacers, insertion_range, contigs=None, spacers=None, prop_rc=0.5, debug=False):
    if contigs is None: # if not provided, generate contigs
        if debug:
            print("Generating contigs...")
            import time
            start_time = time.time()
        # Pre-generate all random lengths at once
        contig_lengths = [random.randint(*contig_length_range) for _ in range(sample_size_contigs)]
        total_length = sum(contig_lengths)
        
        # Generate one large sequence and split it
        large_sequence = ''.join(random.choices('ATCG', k=total_length))
        
        # Split into individual contigs
        contigs = {}
        used_positions_per_contig = {}
        current_pos = 0
        
        for i in range(sample_size_contigs):
            length = contig_lengths[i]
            contig_id = f"contig_{i}"
            contigs[contig_id] = large_sequence[current_pos:current_pos + length]
            used_positions_per_contig[contig_id] = set()
            current_pos += length
    else:
        if debug: print("Reading contigs from file...")
        contigs = read_fasta(contigs)
        used_positions_per_contig = {contig_id: set() for contig_id in contigs.keys()}
        sample_size_contigs = len(contigs)
        if debug: print(f"Read {sample_size_contigs} contigs")
        
    if spacers is None: # if not provided, generate spacers
        if debug: print("Generating spacers...")
        spacers = {}
        # Generate spacers
        for i in tqdm(range(sample_size_spacers), desc="Generating spacers", total=sample_size_spacers):
            spacer_length = random.randint(*spacer_length_range)
            spacer = generate_random_sequence(spacer_length)
            spacer_id = f"spacer_{i}"
            spacers[spacer_id] = spacer
            # if debug and i % 1000 == 0: print(f"Generated {i} spacers...")
    else:
        if debug: print("Reading spacers from file...")
        spacers = read_fasta(spacers)   
        sample_size_spacers = len(spacers)
        if debug: print(f"Read {sample_size_spacers} spacers")
        
    ground_truth = []
    if insertion_range == (0, 0):
        if debug: print("No insertions requested, returning empty ground truth")
        return contigs, spacers, pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "strand": pl.Boolean, "mismatches": pl.UInt32})
    
    if debug: print("Starting spacer insertions...")
    # insert every spacer into random contigs
    for i in tqdm(range(sample_size_spacers), desc="Processing spacer", total=sample_size_spacers):
        n_insertions = random.randint(*insertion_range)
        spacer_id = f"spacer_{i}"
        spacer = spacers[spacer_id]
        if debug and i % 100 == 0: print(f"Processing spacer {i}, making {n_insertions} insertions")
        
        # For each insertion, randomly select a contig
        for j in range(n_insertions):
            target_contig_id = random.choice(list(contigs.keys()))
            target_contig = contigs[target_contig_id]
            
            # Decide whether to reverse complement
            is_rc = random.random() < prop_rc
            spacer_to_insert = reverse_complement(spacer) if is_rc else spacer
            
            # Apply mismatches
            n_mismatches = random.randint(*n_mismatch_range)
            spacer_with_errors = apply_mismatches(spacer_to_insert, n_mismatches)
            
            # Find non-overlapping position for substitution
            while True: # this could probably be made more efficient by selecting non occupied positions. but with a small number of insertions (and small insertion size, and large contig size) I guess the chacnes of overlap are low.
                # Ensure we don't try to substitute beyond contig bounds
                max_pos = len(target_contig) - len(spacer_with_errors)
                if max_pos < 0:  # Skip if spacer is longer than contig
                    if debug: print(f"Warning: spacer {spacer_id} is longer than contig {target_contig_id}, skipping")
                    continue
                    
                start_pos = random.randint(0, max_pos)
                end_pos = start_pos + len(spacer_with_errors)
                position_range = set(range(start_pos, end_pos))
                
                if not position_range.intersection(used_positions_per_contig[target_contig_id]):
                    used_positions_per_contig[target_contig_id].update(position_range)
                    break
            
            # Substitute the spacer into the contig
            contig_list = list(target_contig)
            contig_list[start_pos:end_pos] = spacer_with_errors
            contigs[target_contig_id] = ''.join(contig_list)
            
            # Record the ground truth
            ground_truth.append([
                spacer_id,
                target_contig_id,
                start_pos,
                end_pos,
                is_rc,
                n_mismatches
            ])

    if debug: print("Creating ground truth DataFrame...")
    ground_truth = pl.DataFrame(
        ground_truth,
        schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "strand": pl.Boolean, "mismatches": pl.UInt32}
    )
    if debug:
        end_time = time.time()
        print(f"Simulation complete in {end_time - start_time} seconds")
    return contigs, spacers, ground_truth


from rust_simulator import Simulator
def generate_simulation_id(params):
    """Generate a unique simulation ID based on input parameters.
    
    Args:
        params (dict): Dictionary of simulation parameters
        
    Returns:
        str: 8-character hex string unique to these parameters
    """
    # Sort the parameters to ensure consistent hashing
    param_str = json.dumps(params, sort_keys=True)
    # Generate hash
    hash_obj = hashlib.sha256(param_str.encode())
    # Return first 8 characters of hash
    return hash_obj.hexdigest()[:8]

def simulate_data_rust(contig_length_range, spacer_length_range, n_mismatch_range, 
                      sample_size_contigs, sample_size_spacers, insertion_range, 
                      contigs=None, spacers=None, prop_rc=0.5, debug=False, threads=None, verify=True, results_dir=None, id_prefix=None):
    
    if contigs is not None or spacers is not None:
        return simulate_data(contig_length_range, spacer_length_range, n_mismatch_range, 
                           sample_size_contigs, sample_size_spacers, insertion_range, 
                           contigs, spacers, prop_rc, debug)
    
    if threads is None:
        import multiprocessing
        threads = multiprocessing.cpu_count()
    
    if debug: 
        print(f"Using {threads} threads for Rust simulation...")
        import time
        start_time = time.time()
    
    # Generate simulation ID if not provided
    if id_prefix is None:
        sim_params = {
            'contig_length_range': contig_length_range,
            'spacer_length_range': spacer_length_range,
            'n_mismatch_range': n_mismatch_range,
            'sample_size_contigs': sample_size_contigs,
            'sample_size_spacers': sample_size_spacers,
            'insertion_range': insertion_range,
            'prop_rc': prop_rc
        }
        id_prefix = generate_simulation_id(sim_params)
        print(f"Using generated simulation ID: {id_prefix}")
    else:
        print(f"Using provided ID prefix: {id_prefix}")
    
    simulator = Simulator()
    contigs, spacers, ground_truth = simulator.simulate_data(
        tuple(contig_length_range),
        tuple(spacer_length_range),
        tuple(n_mismatch_range),
        sample_size_contigs,
        sample_size_spacers,
        tuple(insertion_range),
        prop_rc,
        threads,
        verify,
        results_dir,
        id_prefix
    )
    
    if debug:
        end_time = time.time()
        print(f"Simulation time when using Rust: {end_time - start_time} seconds")
    
    # Convert the ground truth data types before creating DataFrame
    processed_ground_truth = [
        [
            row[0],                   # spacer_id (already has prefix from Rust)
            row[1],                   # contig_id (already has prefix from Rust)
            int(row[2]),              # start (int)
            int(row[3]),              # end (int)
            row[4].lower() == "true", # strand (bool)
            int(row[5])               # mismatches (int)
        ]
        for row in ground_truth
    ]
    
    ground_truth_df = pl.DataFrame(
        processed_ground_truth,
        schema={
            "spacer_id": pl.Utf8,
            "contig_id": pl.Utf8,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "strand": pl.Boolean,
            "mismatches": pl.UInt32
        }
    )
    if verify:
         verrify_df = verify_simulated_data(contigs, spacers, ground_truth_df)
         if verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")).height > 0:
            print(verrify_df.filter(pl.col("alignment_test") > pl.col("mismatches")))
            raise ValueError("Simulated data is not correct")
         else:
             print("Simulated data is correct")
    return contigs, spacers, ground_truth_df

def verify_simulated_data(contigs, spacers, ground_truth, return_fraction=False, return_bool=False,return_positive_only=False):
    """Verify the simulated data by checking the ground truth against the contigs and spacers."""
    if not isinstance(ground_truth, pl.DataFrame):
        if isinstance(ground_truth, str):
            ground_truth = pl.read_csv(ground_truth, separator="\t")
        else:
            raise ValueError("ground_truth must be a pandas DataFrame or a string")
    # load seqs
    if isinstance(contigs, str):
        ground_truth = populate_pldf_withseqs_needletail(ground_truth,seqfile=contigs,seqcol="contig_seq",idcol="contig_id",trim_to_region=False,reverse_by_strand_col=False)
    else:
        contig_df = pl.DataFrame(data= {"contig_id": list(contigs.keys()), "contig_seq": list(contigs.values()) })
        ground_truth = ground_truth.join(contig_df, on="contig_id", how="left")
    if isinstance(spacers, str):
        ground_truth = populate_pldf_withseqs_needletail(ground_truth,seqfile=spacers,seqcol="spacer_seq",idcol="spacer_id",trim_to_region=False,reverse_by_strand_col=False)
    else:
        spacer_df = pl.DataFrame(data= {"spacer_id": list(spacers.keys()), "spacer_seq": list(spacers.values()) })
        ground_truth = ground_truth.join(spacer_df, on="spacer_id", how="left")

    # check if the ground truth is correct
    verify_df = test_alignment_polars(ground_truth,ignore_region_strands=False)
    if return_bool:
        return verify_df.filter(pl.col("alignment_test") != pl.col("mismatches")).height == ground_truth.height
    if return_fraction:
        error_df = verify_df.filter(pl.col("alignment_test") != pl.col("mismatches"))
        return error_df.height / ground_truth.height
    elif return_positive_only:
        return verify_df.filter(pl.col("alignment_test") == pl.col("mismatches"))
    return verify_df


def write_fasta(sequences, filename):
    with open(filename, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")

def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip()
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()
    return sequences

def is_it_there(subsequence, sequence):
    return sequence.contains(subsequence)

def is_it_there_fromtbl(mama_seq, baby_seq,start_pos,end_pos):
    return mama_seq[start_pos:end_pos].contains(baby_seq)

def run_tool(tool, results_dir, max_runs=1,warmups=1):
    create_bash_script(tool, results_dir, max_runs=max_runs,warmups=warmups)
    subprocess.run(f'{results_dir}/bash_scripts/{tool["script_name"]}', shell=True)
    
    # time_info = {}
    with open(f"{results_dir}/raw_outputs/{tool['script_name']}.json", 'r') as f:
        data = json.load(f)
        # time_info["User time (seconds)"] = data['results'][0]['mean']
        # time_info["system time (seconds)"] = data['results'][0]['system']
    return data

def create_bash_script(tool, results_dir, max_runs=1,warmups=1):
    script_name = tool['script_name']
    mamba_env = tool.get('mamba_env', None)
    hyperfine_command =  ["hyperfine",
               "--warmup", str(warmups),
               "--max-runs", str(max_runs),
                "--prepare", f'"micromamba deactivate || true && rm -rf  {results_dir}/raw_outputs/tmp* || true"',
               "--output", f"{results_dir}/raw_outputs/hyperfine_output_{tool['script_name']}.txt", "--export-json", f"{results_dir}/raw_outputs/{tool['script_name']}.json"]
    
    with open(f"{results_dir}/bash_scripts/{tool['script_name']}", 'w') as f:
        if mamba_env:
            f.write(f"eval \"$(micromamba shell hook --shell bash)\"\n")
            f.write(f"micromamba activate {mamba_env}\n")
            f.write(f"{' '.join(hyperfine_command)} ")
            f.write(f"'{' '.join(tool['command'])}'")
        else:
            f.write("#!/bin/bash\n")    
            hyperfine_command =  ["hyperfine",
               "--warmup", "1",
                "--prepare", f'"rm -rf  {results_dir}/raw_outputs/tmp* || true"',
               "--output", f"{results_dir}/raw_outputs/hyperfine_output_{tool['script_name']}.txt", "--export-json", f"{results_dir}/raw_outputs/{tool['script_name']}.json"]
            f.write(f"{' '.join(hyperfine_command)} ")
            f.write(f"'{' '.join(tool['command'])}'")

        os.chmod(f"{results_dir}/bash_scripts/{script_name}", 0o755)

def clean_everything(results_dir):
    os.system(f"rm -rf {results_dir}/raw_outputs/")
    os.system(f"rm -rf {results_dir}/bash_scripts/")
    os.system(f"rm -rf {results_dir}/results/")
    os.system(f"rm -rf {results_dir}/simulated_data/")

def clean_before_rerun(tool_name, results_dir):
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*tsv")
    os.system(f"rm -rf {results_dir}/raw_outputs/{tool_name}*sam")
    os.system(f"rm -rf {results_dir}/bash_scripts/{tool_name}*")
    os.system(f"rm -rf {results_dir}/raw_outputs/tmp*")

def insert_spacer(contig, spacer, position):
    # deprecated
    return contig[:position] + spacer + contig[position:]

def get_aln_len_from_cigar(cigar):
    return sum(int(num) for num in re.findall(r'\d+', cigar))

def fix_cigar_for_pysamtools(cigar):
    # if any character is not preceded by a digit, append to it 1
    # handle both standard CIGAR operators (A-Z) and extended operators (=,X)
    t = re.sub(r'([A-Z|=|X])(?!\d)([A-Z|=|X])', r'\g<1>1\g<2>', cigar, count=0)
    if re.search(r'[A-Z|=|X][A-Z|=|X]', t):  # Check for consecutive operators
        return fix_cigar_for_pysamtools(t)
    else:
        return t

def get_mismatches_from_cigar(cigar):
    return sum(int(num[:-1]) for num in re.findall(r'\d+X', cigar) or [0])

def parse_samVext(sam_file, max_mismatches=5):
    """Parse SAM file using the version 1.4 cigar strings, where = for matches, X for mismatches, I for insertions, D for deletions."""
    results = []
    with open(sam_file, 'r') as f:
        for line in f:
            # print(line)
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "*":  # unmapped
                continue
            if len(fields) >= 3:
                spacer_id = str(fields[0])
                flag = int(fields[1])
                is_reverse = bool(flag & 16)
                strand = is_reverse  # Changed from '-' if is_reverse else '+'
                contig_id = str(fields[2])
                start = int(fields[3]) - 1  # SAM is 1-based
                cigar = fields[5]
                seq_length = get_aln_len_from_cigar(cigar)
                mismatches = -1 # default value
                if len(fields) >= 12: # has tags
                    # get which tag is NM:i:``
                    nm_tag = [x for x in fields[11:] if x.startswith('NM:i:')]
                    if len(nm_tag) > 0:
                        mismatches = int(nm_tag[0].split(':')[2])
                        
                if mismatches == -1: # if we don't find NM:i: we try to get it from the cigar string.
                    mismatches = get_mismatches_from_cigar(cigar)

                if mismatches > max_mismatches:
                    continue
                end = start + seq_length 
                results.append([
                    str(spacer_id),
                    str(contig_id),
                    str(strand),
                    int(start),
                    int(end),
                    int(mismatches)
                ])

    if len(results) == 0:
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32}) # empty dataframe, so we can still join with other results
    
    # Create lists for each column
    spacer_ids, contig_ids, strands, starts, ends, mismatches = zip(*results)
    
    df = pl.DataFrame({
        "spacer_id": spacer_ids,
        "contig_id": contig_ids,
        "strand": strands,
        "start": starts,
        "end": ends,
        "mismatches": mismatches
    })
    return df.unique()

def parse_samV1(sam_file, max_mismatches=5):
    """Parse SAM file and return standardized coordinates (1-based).
    SAM format uses 1-based coordinates. 
    Filter out rows with more than max_mismatches mismatches. (i.e. retain up to (including) max_mismatches mismatches)"""
    results = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "*":  # unmapped
                continue
            if len(fields) >= 3:
                spacer_id = fields[0]
                flag = int(fields[1])
                is_reverse = bool(flag & 16)
                strand = is_reverse  # Changed from '-' if is_reverse else '+'
                contig_id = fields[2]
                start = int(fields[3]) - 1  # SAM is 1-based
                cigar = fields[5]
                seq_length = get_aln_len_from_cigar(cigar)
                mismatches = -1 # default value, to raise suspicion later if we don't find NM:i: and the sam version is not v1.4
                if len(fields) >= 12: # has tags
                    # get which tag is NM:i:``
                    nm_tag = [x for x in fields[11:] if x.startswith('NM:i:')]
                    if len(nm_tag) > 0:
                        mismatches = int(nm_tag[0].split(':')[2])
                        # Mummer doesn't report terminal gaps (clipping) as mismatches, so we will test if the alignment is clipped and if so, add the number of clipped bases as mismatches
                        clipping = re.findall(r'\d+S', cigar)
                        if clipping:
                            # get the number of clipped bases
                            clipped_bases = sum(int(num[:-1]) for num in clipping)
                            mismatches += clipped_bases # can not figure out how did mummer calculate this as 0 or 1 mis: 88755144\t16\tIMGVR_UViG_2571042239_000005|2571042239|2571059887|117-31391\t14196\t30\t3S26M\t*\t0\t0\ttattaacctcctttctagctaccaaataa\t*\tNM:i:1\tMD:Z:5a14215

                        
                if mismatches > max_mismatches: 
                    continue
                # Calculate end position based on CIGAR string, probably wrong...
                end = start + seq_length 
                # For reverse complement, coordinates are already correct in SAM
                results.append([spacer_id, contig_id, strand, start, end, mismatches]) # maybe later add mismatches to the results
    if len(results) == 0:
        return pl.DataFrame(schema=["spacer_id", "contig_id", "strand", "start", "end", "mismatches"]) # empty dataframe, so we can still join with other results
    return pl.DataFrame(results, schema=["spacer_id", "contig_id", "strand", "start", "end", "mismatches"]).unique()

def order_columns_to_match(df1_to_order, df2_to_match):
    return df1_to_order.select(df2_to_match.columns)

def cast_cols_to_match(df1_to_cast, df2_to_match):
    for col in df2_to_match.columns:
        df1_to_cast = df1_to_cast.with_columns(pl.col(col).cast(df2_to_match.schema[col]))
    return df1_to_cast

def vstack_easy(df1_to_stack, df2_to_stack):
    df2_to_stack = cast_cols_to_match(df2_to_stack, df1_to_stack)
    df2_to_stack = order_columns_to_match(df2_to_stack, df1_to_stack)
    return df1_to_stack.vstack(df2_to_stack)

def parse_blastn_custom(blastn_file, max_mismatches=5, spacer_lendf=None, **kwargs):
    """Parse a custom BLAST or mmseqs output format 
    (query,target,nident,alnlen,mismatch,qlen,gapopen,qstart,qend,tstart,tend,evalue,bits)
    and return standardized coordinates (1-based).
    BLAST format uses 1-based coordinates. 
    Filter out rows with more than max_mismatches mismatches. (i.e. retain up to (including) max_mismatches mismatches)"""
    try:
        results = pl.read_csv(
            blastn_file, 
            separator="\t", 
            has_header=False,
            infer_schema_length=100000,
            # n_threads=kwargs.get("threads", 1),
            new_columns=["spacer_id", "contig_id", "nident", "alnlen", "mismatch", "qlen", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits"]
        )
        if any(results["qend"] < results["qstart"]): # mmseqs reports reverse on the subject, not the query.
            print("qend < qstart detected, assuming tend and tstart are reversed")
            results = results.with_columns(
                pl.when(pl.col("qend") < pl.col("qstart"))
                .then(pl.struct(
                    tstart=pl.col("tend"),
                    tend=pl.col("tstart")
                ))
                .otherwise(pl.struct(
                    tstart=pl.col("tstart"),
                    tend=pl.col("tend")
                ))
                .struct.field(["tstart", "tend"])
            )
        # results = results.with_columns(
        #     pl.col("spacer_id").cast(pl.Utf8)
        # )
    except Exception as e:
        print(f"Failed to read BLAST file {blastn_file}: {e}, returning empty dataframe")
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})
    
    results = results.with_columns(
        pl.col("spacer_id").cast(pl.Utf8)
    )
    results = results.with_columns(
        # Determine strand based on coordinate ordering
        pl.when(pl.col("tstart") <= pl.col("tend"))
        .then(pl.lit(False))  # Changed from '+' to False
        .otherwise(pl.lit(True))  # Changed from '-' to True
        .alias('strand'),
        pl.min_horizontal(["tstart", "tend"]).alias("start_min")-1,
        pl.max_horizontal(["tstart", "tend"]).alias("end_max"),
        # Calculate mismatches as edit distance (non-identical positions + gaps)
         (pl.col("qlen")-pl.col("nident")).alias("mismatches") # Need to consider setting this as max(qlen,alnlen), as mmseqs can report seems to allow gaps in the query to be counted in the aligned length.
    )
    # print(results.filter(pl.col("mismatches") != pl.col("mismatch")).height)
    
    # Filter out rows with more than max_mismatches mismatches
    results = results.with_columns(pl.col("mismatches") + pl.col("gapopen").alias("mismatches"))
    results = results.filter(max_mismatches >= pl.col("mismatches"))
    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.rename({"length": "spacer_length"})

    results = results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        pl.col("start_min").alias("start"),
        pl.col("end_max").alias("end"),
        "mismatches"
    ).unique()
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    return results

def parse_lexicmap(tsv_file, max_mismatches=5, spacer_lendf=None, **kwargs):
    """Parse LexicMap TSV output format and return standardized coordinates.
    
    LexicMap output columns (1-based positions):
    1. query      - Query sequence ID
    2. qlen       - Query sequence length
    3. hits       - Number of subject genomes
    4. sgenome    - Subject genome ID
    5. sseqid     - Subject sequence ID
    6. qcovGnm    - Query coverage per genome
    7. hsp        - HSP number
    8. qcovHSP    - Query coverage per HSP
    9. alenHSP    - Aligned length in HSP
    10. pident    - Percentage identity
    11. gaps      - Gaps in HSP
    12. qstart    - Query start
    13. qend      - Query end
    14. sstart    - Subject start
    15. send      - Subject end
    16. sstr      - Subject strand
    17. slen      - Subject length
    18+ Optional columns with -a flag (cigar, qseq, sseq, align)
    """
    try:
        results=pl.read_csv(tsv_file, separator="\t", has_header=True,
                            new_columns=["spacer_id", "qlen", "hits", "sgenome", "contig_id", "qcovGnm", "hsp", "qcovHSP", "alenHSP", "pident", "gaps", "qstart", "qend", "sstart", "end", "strand", "slen","cigar", "qseq", "sseq", "align"],
                            infer_schema_length=100000)
    except Exception as e:
        print(f"Failed to create read from file {tsv_file}: {e}, returning empty dataframe")
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})
            
    results = results.filter(pl.col("alenHSP") >= 17, pl.col("gaps") == 0)
    results = results.with_columns(pl.col("spacer_id").cast(pl.Utf8))
    results = results.with_columns(
        pl.col("align").str.count_matches("|",literal=True).alias("matches"))
    
    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.with_columns(
        (pl.col("length")-pl.col("matches")).alias("mismatches"),
        (pl.col("sstart")-1).alias("start")
    )
    results = results.filter(pl.col("mismatches") <= max_mismatches)
    results = results.rename({"length": "spacer_length"})
    results = results.select("spacer_id", "contig_id", "spacer_length", "strand", "start", "end", "mismatches")
    results = results.with_columns(
        pl.when(pl.col("strand") == "-")
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("strand")
    )
        
    return results.select("spacer_id", "contig_id", "spacer_length", "strand", "start", "end", "mismatches")

def parse_samV1_with_lens(sam_file, spacer_lendf, max_mismatches=5, **kwargs):
    """Parse SAM file using spacer lengths from reference table to compute mismatches.
    spacer_lens_df should be a polars DataFrame with columns: spacer_id, length"""
    results = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "*":  # unmapped
                continue
            if len(fields) >= 3:
                spacer_id = fields[0]
                # Get query sequence length from spacer_lens_df
                spacer_in_df = spacer_lendf.filter(pl.col("spacer_id") == spacer_id)
                if spacer_in_df.height == 0:
                    # print(f"Spacer {spacer_id} not found in spacer_lendf")
                    continue
                flag = int(fields[1])
                is_reverse = bool(flag & 16)
                strand = is_reverse  # Changed from '-' if is_reverse else '+'
                contig_id = fields[2]
                start = int(fields[3]) - 1  # SAM is 1-based
                cigar = fields[5]
                seq_length = get_aln_len_from_cigar(cigar)
                spacer_len = spacer_in_df["length"].unique().item()
                
                # Count matches from CIGAR string
                matches = sum(int(num) for num in re.findall(r'(\d+)[=M]', cigar))
                
                # Calculate mismatches as spacer length minus matches
                mismatches = spacer_len - matches
                
                if mismatches > max_mismatches:
                    continue
                    
                end = start + seq_length
                results.append([
                    str(spacer_id),
                    str(contig_id),
                    str(strand),
                    int(start),
                    int(end),
                    int(mismatches)
                ])

    if len(results) == 0:
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})
    
    return pl.DataFrame(
        results,
        schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32}
    ).unique()

def parse_samVext_with_lens(sam_file, spacer_lendf, max_mismatches=5, **kwargs):
    """Parse SAM file using spacer lengths and extended CIGAR (=X) to compute mismatches."""
    results = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == "*":  # unmapped
                continue
            if len(fields) >= 3:
                spacer_id = str(fields[0])
                # Get query sequence length from spacer_lens_df
                spacer_in_df = spacer_lendf.filter(pl.col("spacer_id") == spacer_id)
                if spacer_in_df.height == 0:
                    # print(f"Spacer {spacer_id} not found in spacer_lendf")
                    continue
                flag = int(fields[1])
                is_reverse = bool(flag & 16)
                strand = is_reverse  # Changed from '-' if is_reverse else '+'
                contig_id = str(fields[2])
                start = int(fields[3]) - 1  # SAM is 1-based
                cigar = fields[5]
                seq_length = get_aln_len_from_cigar(cigar)
                
                # Count matches from extended CIGAR string
                matches = sum(int(num) for num in re.findall(r'(\d+)[=M]', cigar))
                
                spacer_len = spacer_in_df["length"].unique().item()
                
                # Calculate mismatches as spacer length minus matches
                mismatches = spacer_len - matches

                if mismatches > max_mismatches:
                    continue
                    
                end = start + seq_length
                results.append([
                    str(spacer_id),
                    str(contig_id),
                    str(strand),
                    int(start),
                    int(end),
                    int(mismatches)
                ])

    if len(results) == 0:
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})

    return pl.DataFrame(
        results,
        schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32}
    ).unique()

def parse_samVn_with_lens_polars(sam_file, spacer_lendf, max_mismatches=5, sam_version="1"):
    """Parse SAM file using spacer lengths for qlen and numerics in CIGAR for alignment length
    uses polars to speed up things"""
    # First read mandatory columns
    mandatory_sam_cols = ["spacer_id", "flag", "contig_id", "start", "MAPQ", "cigar", "rnext", "pnext", "tlen", "seq", "qual"]
    
    try:
        # Read first 11 mandatory columns
        sam_df = pl.read_csv(
            sam_file, 
            separator="\t", 
            comment_prefix="@", 
            has_header=False,
            new_columns=mandatory_sam_cols,
            columns=list(range(11)),
            truncate_ragged_lines=False,
            infer_schema_length=1000
            )
        # add index column
        sam_df = sam_df.with_row_index(name="index")
        sam_df = sam_df.drop(["rnext", "pnext", "MAPQ", "qual"])
        # get the index of the unmapped lines
        unmmaped_lines_idx = sam_df.filter(pl.col("contig_id") == "*").get_column("index")
        if unmmaped_lines_idx.len() > 0:
            print(f"Warning: {unmmaped_lines_idx.len()} unmapped lines found in {sam_file}")
            sam_df = sam_df.filter(pl.col("contig_id") != "*")
    
        # Now read optional tags as a single column
        tags_df = pl.read_csv(
            sam_file,
            separator="%", # risky...
            comment_prefix="@",
            has_header=False,
            new_columns=["tags"],
            truncate_ragged_lines=False
            )
        tags_df = tags_df.with_row_index(name="index")
        if unmmaped_lines_idx.len() > 0:
            tags_df = tags_df.filter(~pl.col("index").is_in(unmmaped_lines_idx))
        
        # tags_df = tags_df.with_columns(
        #     pl.col("tags").str.split( "\t").list.gather_every(n=1,offset=11).alias("tags")
        # )
        # tags_df = tags_df.with_columns(
        #     pl.col("tags").list.len().alias("tag_lengths")
        # )
        nm_tag = "nm:i:" if tags_df["tags"][1].count("NM:i:") == 0 else "NM:i:"
        tags_df = tags_df.with_columns(
            pl.col("tags").str.extract(f"{nm_tag}(\d+)").cast(pl.Int32).alias("nm")
        )
        # Join with main dataframe
        sam_df =sam_df.with_columns(tags_df["nm"].alias("nm"))
        sam_df = sam_df.drop(["index"])
    except Exception as e:
        print(f"Failed to read SAM file {sam_file}: {e}, returning empty dataframe")
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})

    # Continue with existing processing
    sam_df = sam_df.with_columns(
        (pl.col("start")-1).alias("start")
    )
    sam_df = sam_df.with_columns([
        pl.when(pl.col("flag") == 16)
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias('strand')
    ])
    sam_df = sam_df.drop("flag")

    sam_df = sam_df.with_columns([
        pl.col("cigar").map_elements(lambda x: sum(int(num) for num in re.findall(r'(\d+)', x)), return_dtype=pl.UInt32).alias("aligned_length"),
        pl.col("cigar").map_elements(lambda x: sum(int(num) for num in re.findall(r'(\d+)S', x)), return_dtype=pl.UInt32).alias("soft_clipped_length")
    ])

    # Join with spacer lengths
    sam_df = spacer_lendf.join(sam_df, on="spacer_id", how="inner")
    
    # Calculate mismatches based on length and matches
    sam_df = sam_df.with_columns([
        ((pl.col("length")-pl.col("aligned_length"))+pl.col("soft_clipped_length")+pl.col("nm")).alias("mismatches")
    ])

    sam_df = sam_df.with_columns(
       (pl.col("start")+pl.col("aligned_length"))
        .alias("end")
    )
    sam_df = sam_df.rename({"length": "spacer_length"})
    sam_df = sam_df.filter(pl.col("mismatches") <= max_mismatches)
    sam_df = sam_df.filter(pl.col("nm") <= max_mismatches)
    sam_df = sam_df.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches"
    ).unique()
    
    return sam_df

def parse_samVn_with_lens_pysam(sam_file, spacer_lendf, max_mismatches=5, threads=4, ref_file=None, **kwargs):
    """Parse SAM file using pysam and spacer lengths to compute mismatches"""
    try: # if no headers, we need to find the smallest sam file in the output file directory and copy its header
        sam_file = add_sqlines_sam(sam_file,ref_file)
        results = []
        
        # Create a dictionary for quick spacer length lookups
        spacer_lens = dict(zip(spacer_lendf["spacer_id"], spacer_lendf["length"]))
        
        with pysam.AlignmentFile(sam_file, "r", check_sq=False) as samfile:
            for read in tqdm(samfile, desc="Parsing SAM file"):
                if read.is_unmapped:
                    continue

                # spacer_len = read.query_length # 213768.22it/s
                # read.cigarstring = fix_cigar_for_pysamtools(read.cigarstring) # if missing integars in the cigar string...
                spacer_len = spacer_lens.get(read.query_name) #  199485, with the "fix_cigar_for_pysamtools" ~99k.

                # spacer_len = read.infer_read_length() #
                # ## Get hard clipped length
                # hard_clipped = sum(length for op, length in read.cigartuples 
                #                  if op == 5)  # H operation
                # # if hard_clipped > 0:
                # #     print(f"Hard clipped {hard_clipped} for {read.query_name}")
                # #     # break
                # if spacer_len == 0: # some alignments get 0 length, maybe because some tools do not output the sequence for secondary/ambiguous alignments. I suspect pysam get's the query length from the nchar of the seq.
                #     # print (f"spacer_len == 0 for {read.query_name}, extracting alignment info from cigar")
                # # Get spacer length from reference table instead of query sequence
                #     spacer_len = spacer_lens.get(read.query_name)
                # else:
                #     spacer_len = spacer_len + hard_clipped
                #     # print(f"{spacer_len}")
                #     # break
                if spacer_len is None:
                    print(f"Warning: spacer {read.query_name} not found in reference table")
                    break
                # Get alignment infto_string()
                # strand = '-' if read.is_reverse else '+'
                start = read.reference_start  # pysam uses 0-based coordinates
                
                # Get alignment length excluding soft clips
                # query_alignment_length = read.query_alignment_length
                query_alignment_length = len(read.get_reference_positions(full_length=False))
                # aligned_length = sum(length for op, length in read.cigartuples 
                #                   if op in [0, 7, 8])  # M, =, X operations
                
                # # # Get soft clipped length
                soft_clipped = sum(length for op, length in read.cigartuples 
                                 if op == 4)  # S operation

                # if soft_clipped > 0:
                #     print(f"Soft clipped {soft_clipped} for {read.query_name}")
                #     break
                # Get edit distance from NM tag  (counts non-identical positions including gaps)
                nm = (read.get_tag('NM') if read.has_tag('NM') else 0) + soft_clipped
                
                # Calculate mismatches based on length and matches
                mismatches = (spacer_len - query_alignment_length) + nm
                if mismatches <0:
                    print("mismatches < 0")
                    break
                if mismatches > max_mismatches or nm > max_mismatches:
                    continue

                results.append([
                    str(read.query_name),
                    str(read.reference_name),
                    int(spacer_len),
                    read.is_reverse,
                    int(start),
                    int(read.reference_end), # read.reference_end
                    int(mismatches)
                ])

    except Exception as e:
        print(f"Failed to read SAM file {sam_file}: {e}, returning empty dataframe")
        return pl.DataFrame(
            schema={
                "spacer_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "spacer_length": pl.UInt32,
                "strand": pl.Boolean,
                "start": pl.UInt32,
                "end": pl.UInt32,
                "mismatches": pl.UInt32
            }
        )
    return pl.DataFrame(
        results,
        schema={
            "spacer_id": pl.Utf8,
            "contig_id": pl.Utf8,
            "spacer_length": pl.UInt32,
            "strand": pl.Boolean,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "mismatches": pl.UInt32
        }
    ).unique()


def plot_matrix(matrix, title, filename, hor_axis="pairs in tool", vert_axis="pairs not in tool", return_chart=True):
    import altair as alt

    # Rename the tool1 column to 'in_tool'
    matrix = matrix.rename({"tool1": "in_tool"})
    
    # Calculate total counts for each tool to use for sorting
    tool_cols = [col for col in matrix.columns if col != 'in_tool']
    matrix = matrix.with_columns(
        pl.sum_horizontal(tool_cols).alias("total")
    )
    
    # Sort the matrix by total counts (descending)
    matrix = matrix.sort("total", descending=True)
    sorted_tools = matrix["in_tool"].to_list()
    
    # Drop the total column after sorting
    matrix = matrix.drop("total")
    
    # Reorder the columns to match the sorted rows
    ordered_cols = ["in_tool"] + sorted_tools
    matrix = matrix.select([col for col in ordered_cols if col in matrix.columns])
    
    # Convert to percentages if requested
    # if as_percent:
    #     # Get total matches for each tool
    #     tool_totals = matrix.group_by('in_tool').agg(
    #         pl.len().alias('total_matches')
    #     ).to_dict(as_series=False)
    #     tool_totals = dict(zip(tool_totals['in_tool'], tool_totals['total_matches']))
        
    #     # Convert counts to percentages
    #     for tool in sorted_tools:
    #         if tool != 'in_tool':
    #             matrix = matrix.with_columns(
    #                 (pl.col(tool) * 100 / tool_totals[tool]).alias(tool)
    #             )
    
    # Convert the matrix to long form for Altair
    melted = matrix.melt(
        id_vars=["in_tool"],
        value_vars=[col for col in matrix.columns if col != 'in_tool'],
        variable_name="not_in_tool",
        value_name="count"
    )
    
    # Create the heatmap with Altair
    chart = alt.Chart(melted).mark_rect().encode(
        x=alt.X('not_in_tool:N', 
                title=vert_axis,
                sort=sorted_tools,
                axis=alt.Axis(labelAngle=-45,labelFontSize=14, titleFontSize=16)),
        y=alt.Y('in_tool:N', 
                title=hor_axis,
                sort=sorted_tools,
                axis=alt.Axis(labelFontSize=14, titleFontSize=16)),
        color=alt.Color('count:Q', 
                       scale=alt.Scale(scheme='blues'),
                       legend=alt.Legend(title="Count",labelFontSize=14, titleFontSize=16)),
        tooltip=['in_tool', 'not_in_tool', 'count']
    ).properties(
        title=title,
        width=600,
        height=600
    )
    
    # Add text labels on the cells
    text = chart.mark_text(baseline='middle').encode(
        text=alt.Text('count:Q', format='d'),
        color=alt.condition(
            alt.datum.count > melted['count'].mean(),
            alt.value('white'),
            alt.value('black')
        )
    )
    
    # Combine the heatmap and text
    final_chart = (chart + text)
    
    # Save both HTML and SVG versions
    final_chart.save(filename + '.html')
    
    # Configure for SVG export with higher DPI
    svg_chart = final_chart.configure_view(
        strokeWidth=0.1,  # Remove border
    ).properties(
        width=600,    # Double size for higher resolution
        height=600
        
    )
    
    # Save as SVG with higher DPI
    svg_chart.save(filename + '.svg', scale_factor=1.0)
    
    return final_chart if return_chart else None


def parse_tab(output_file, max_mismatches=5, spacer_lendf=None, **kwargs): # expecting tab file with columns: spacer_id, contig_id, start, end, strand - because Antonio insists.
    try:
        results = pl.read_csv(output_file, separator="\t", has_header=False, 
                              new_columns=["contig_id","spacer_id", "start", "end","strand"],
                              infer_schema_length=100000)
    except Exception as e:
        print(f"Failed to read TAB file {output_file}: {e}, returning empty dataframe")
        return pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Utf8, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32})
    
    results = results.with_columns(
        pl.lit(0).alias("mismatches") # we assume this function is only used for spacer_containment, so we don't have mismatches
    )
    results = results.with_columns(
        pl.col("spacer_id").cast(pl.Utf8)
    )
    results = spacer_lendf.join(results, on="spacer_id", how="inner")
    results = results.rename({"length": "spacer_length"})
    results = results.with_columns(
        pl.when(pl.col("strand") == "-")
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias("strand")  # Convert string strand to boolean
    )
    return results.select(
        "spacer_id",
        "contig_id",
        "spacer_length",
        "strand",
        "start",
        "end",
        "mismatches"
    ).unique()

def parse_hyperfine_output(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Failed to read hyperfine results for {json_file}: {e}")
        return None
    
    # Extract the first result (assuming one script-ran per file)
    result = data['results'][0]
    
    # Extract relevant timing information
    time_info = pl.DataFrame({
        "mean_time": result['mean'],
        "stddev_time": result['stddev'],
        "median_time": result['median'],
        "user_time": result['user'],
        "system_time": result['system'],
        "min_time": result['min'],
        "max_time": result['max'],
        "command": result['command'],
        "tool": json_file.split("/")[-1].replace(".sh.json", "")
    })
    return time_info

def read_hyperfine_results(tools, results_dir):
    results = pl.DataFrame(schema={"mean_time": pl.Float64, "stddev_time": pl.Float64, "median_time": pl.Float64, "user_time": pl.Float64, "system_time": pl.Float64, "min_time": pl.Float64, "max_time": pl.Float64, "command": pl.Utf8, "tool": pl.Utf8})
    for tool in tools.values():
        try:
            tool_results = parse_hyperfine_output(f"{results_dir}/raw_outputs/{tool['script_name']}.json")
            results = results.vstack(tool_results)
        except Exception as e:
            print(f"Failed to read hyperfine results for {tool['name']}: {e}")
    return results

def run_tools(tools, results_dir, max_runs=1,warmups=1):
    for tool in tools.values():
        try:
            print(f"\nRunning {tool['name']}...")
            # extra_args = tool.get('extra_args', {})
            clean_before_rerun(tool['name'], results_dir)
            tool['time_info'] = run_tool(tool, results_dir, max_runs=max_runs,warmups=warmups)
            print(f"\nSuccessfully ran tool: {tool['name']}")
        except Exception as e:
            print(f"{tool['name']} failed: {e}")

def get_seq_from_fastx(seqfile, seq_ids, return_dict=False, return_df=True,idcol="contig_id", seqcol="contig_seq"):
    if isinstance(seq_ids, str):
        seq_ids = [seq_ids]
    if isinstance(seq_ids, pl.DataFrame):
        seq_ids = seq_ids["contig_id"].to_list()
    # make a tmp file in memory with the seq_ids
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as tmp_file:
        tmp_file.write("\n".join(seq_ids))
        tmp_file_path = tmp_file.name
    seqs = subprocess.run(f"pyfastx extract {seqfile} -l '{tmp_file_path}' | seqkit seq -w 0", shell=True, capture_output=True).stdout.decode("utf-8").replace(">","").splitlines()
    tmp_file.close()
    os.remove(tmp_file_path)
    if not return_df and not return_dict:
        return seqs
    found_ids = [x.removeprefix(">") for x in seqs[::2]]
    seq_dicts = dict(zip(found_ids, seqs[1::2]))
    if return_dict:
        return seq_dicts
    elif return_df:
        return pl.DataFrame({idcol: found_ids, seqcol: seqs[1::2]})


def soft_fetch_fastx(fa, id, start=None, end=None, strand=None):
    if start != None and end != None and strand != None:
        try:
            return fa.fetch(id, (start, end), strand)
        except Exception as e:
            return None
    else:
        try:
            return fa[id].seq   
        except Exception as e:
            return None # if the id is not found, return None

def populate_pldf_from_fastx(seqfile, pldf : pl.DataFrame, trim_to_region=False, idcol="contig_id",  start_col="start", end_col="end", strand_col="strand", seqcol="contig_seq",drop_missing=True, **kwargs):
    """ Get sequences from a fasta file using pyfastx.
    sequences will be extracted from the start-1 to the end position.
    If strand_col is provided, the extracted sequences will be reverse complemented if the strand is boolean True.
    For some reason, this is slower than the populate_pldf_withseqs function.
    """
    # og_cols = pldf.columns
    fa = pfx.Fasta(seqfile)
    cols = start_col is not None and end_col is not None and strand_col is not None
    if ((not cols) or (trim_to_region)):
        print("No start, end, or strand column provided, or trim_to_region is True, returning the entire sequences")
        tmpdf = pldf.select([idcol]).unique()
        tmpdf = tmpdf.with_columns(
            pl.lit(None).alias(start_col),
            pl.lit(None).alias(end_col),
            pl.lit(None).alias(strand_col)
        )
        tmpdf = tmpdf.with_columns(
            pl.col(idcol).map_elements(lambda x: soft_fetch_fastx(fa, x, None, None, None))
        )
        return tmpdf
    # temporary dataframe with the seq_ids, start, end, and strand taht we will export as a bed file (0 based coordinates, so we need to subtract 1 from the start and end)
    tmpdf = pldf.select([idcol, start_col, end_col, strand_col]).unique()
    tmpdf = tmpdf.with_columns(
        pl.when(pl.col(strand_col))
        .then(pl.lit("-"))
        .otherwise(pl.lit("+"))
        .alias("strand4bed"),
        (pl.col(start_col)-1).alias("start_0based")
    )
    # the following is not the best way to fetch sequences... Need to think of a better way to do this...
    seqs = [None]*tmpdf.height # initialize the list of sequences the same size as the tmpdf
    for ix in tqdm(range(tmpdf.height), desc="Fetching sequences",total=tmpdf.height): # 9998/9998 [00:24<00:00, 410.95it/s, 00:15<00:00, 645.96it/s when list is initialized with None
        row = tmpdf[ix]
        seq = soft_fetch_fastx(fa, row[idcol].item(), row["start_0based"].item(), row[end_col].item(), row["strand4bed"].item())
        seqs[ix] = seq
    tmpdf = tmpdf.with_columns(pl.Series(seqs).alias(seqcol))
    # # benchmarked this against the following:
    # time_start = time.time()
    # tmpdf = tmpdf.with_columns(
    # pl.struct(pl.col(idcol), pl.col("start_0based"), pl.col(end_col), pl.col("strand4bed")).map_elements(lambda x: soft_fetch_fastx(fa, x[idcol], x["start_0based"], x[end_col], x["strand4bed"]),strategy="threading", return_dtype=pl.Utf8).alias(seqcol)
    #     )
    # time_end = time.time() # 17.89 seconds (no threading) 16.66 seconds (threading)
    # print(f"Time taken: {time_end - time_start:.2f} seconds")
    
    tmpdf = tmpdf.drop(["strand4bed", "start_0based"])

    if drop_missing:
        tmpdf = tmpdf.filter(~pl.col(seqcol).is_in([None, "", "nan"]))
    pldf = pldf.join(tmpdf, on=[idcol, start_col, end_col, strand_col], how="left")
    return pldf


# # # test
# one_mismatch = pl.read_parquet('results/real_data/results/one_mismatch.parquet')
# seqfile = 'results/real_data/results/matched_contigs.fna'
# pldf = one_mismatch[range(10000)]
# idcol = "contig_id"
# seqcol = "contig_seq"
# start_col = "start"
# end_col = "end"
# strand_col = "strand"
# test = populate_pldf_from_fastx (seqfile=seqfile, pldf=pldf, idcol=idcol, start_col=start_col, end_col=end_col, strand_col=strand_col, seqcol=seqcol)

def populate_pldf_withseqs(pldf : pl.DataFrame, seqfile, chunk_size=2000, drop_missing=True, trim_to_region=False, reverse_by_strand_col=False, idcol="contig_id",seqcol="contig_seq", strand_col="strand", **kwargs):
    # Note! if a "strand" column is present, and reverse_by_strand_col is True, it will be used to get the reverse complement of the sequence.
    # If drop_missing is True, it will drop any contig_ids that are not present in the sequence file. If no seqids were found, it will return the input dataframe.
    # NOTE: trim_to_region and reverse_by_strand_col are bugged when using slice, even though it's probably much faster.
    # get the unique contig_ids and strands
    stranded = strand_col in pldf.columns and reverse_by_strand_col
    # if stranded:
    #     nr_contigids = pldf[[idcol,strand_col]].unique()
    #     # print("rev true")
    # else:
    nr_contigids = pl.DataFrame({idcol: pldf[idcol].unique()})
        # print("rev false")
    # drop nulls nas or empty strings
    nr_contigids = nr_contigids.filter(~pl.col(idcol).is_in([None, "", "nan"]))
    print(f"Total of unique {idcol} ids to fetch: {nr_contigids.height}")
    nr_contigs_list = nr_contigids[idcol].unique().to_list()
    minipldf = pl.DataFrame({
        idcol: nr_contigs_list,
        seqcol: [None] * len(nr_contigs_list)
    })
    # Process chunks and update DataFrame
    for i in tqdm(range(0, len(nr_contigs_list), chunk_size), desc="Populating df with sequences", total=len(nr_contigs_list)//chunk_size):
        chunk = nr_contigs_list[i:i+chunk_size]
        # Get sequences for current chunk
        chunk_seqs = get_seq_from_fastx(seqfile=seqfile, seq_ids=chunk, return_df=True, idcol=idcol, seqcol=seqcol)
        
        # Update minipldf with new sequences using join
        minipldf = (minipldf
                    .join(chunk_seqs, on=idcol, how="left")
                    .with_columns(
                        # Coalesce keeps existing values where not null, uses new values only where null
                        pl.coalesce([pl.col(f"{seqcol}"), pl.col(f"{seqcol}_right")]).alias(seqcol)
                    )
                    .drop(f"{seqcol}_right"))
    
    print(f"Estimated size of temporary dataframe: {minipldf.estimated_size('mb'):.2f} MB")
    # count the number of nulls in the seqs_list
    nulls = [x for x in minipldf[seqcol].to_list() if x is None]
    print(f"Number of unfound ids: {len(nulls)}")
    if minipldf.height == 0:
        print("No sequences found, returning input as is")
        return pldf
    
    if drop_missing:
        print("Dropping unfound ids")
        pldf = pldf.filter(pl.col(idcol).is_in(minipldf[idcol].unique().to_list()))

    pldf = pldf.join(minipldf, on=idcol, how="left")
        
    if trim_to_region and "start" in pldf.columns and "end" in pldf.columns:
        print("Trimming to region")
        # Benchmarking with slice
        # time_start = time.time()
        # pldf = pldf.with_columns(
            # pl.col(seqcol).str.slice(length= pl.col("end")-pl.col("start"),offset=pl.col("start")-1).alias(seqcol)
        # )
        # time_end = time.time()
        # print(f"Time taken: {time_end - time_start:.2f} seconds") # Time taken: 0.15 seconds
        # Benchmarking with map_elements
        # time_start = time.time()
        pldf = pldf.with_columns(
            pl.struct(pl.col(seqcol), pl.col("start"), pl.col("end")).map_elements(lambda x: str(x[seqcol][x["start"]:x["end"]]),return_dtype=pl.Utf8).alias(seqcol)
        )
        # time_end = time.time()
        # print(f"Time taken: {time_end - time_start:.2f} seconds") # Time taken: 0.06 seconds

    if stranded:
        print("Reversing stranded sequences")
        #get the strand from the nr_contigids
        pldf = pldf.with_columns(
            pl.when(pl.col("strand"))
            .then(pl.col(seqcol).map_elements(reverse_complement, return_dtype=pl.Utf8))
            .otherwise(pl.col(seqcol))
            .alias(seqcol)
        )
        
    return pldf


    #debug
    # test = tools_results[range(100000)]
    # test2= populate_pldf_withseqs(seqfile=seqfile,chunk_size=12000, pldf=test, drop_missing=True, reverse_by_strand_col=False)
    # pldf=tools_results[range(10000)]
    # and add a few a null, empty string, or nan
    # pldf_null = pldf[range(4)].with_columns(
    #     pl.lit(None).alias("contig_id"),
    # )
    # pldf_empty = pldf[range(4)].with_columns(
    #     pl.lit("").alias("contig_id"),
    # )
    # pldf = pldf.vstack(pldf_null).vstack(pldf_empty)
    
    # pldf = pldf.vstack()
def populate_pldf_withseqs_needletail(pldf : pl.DataFrame, seqfile, chunk_size=20000000, trim_to_region=True, reverse_by_strand_col=True, idcol="contig_id",seqcol="contig_seq", start_col="start", end_col="end", strand_col="strand"): 
    merge_cols = [idcol]
    if reverse_by_strand_col:
        merge_cols.append(strand_col)
    if trim_to_region:
        merge_cols.extend([start_col, end_col])
    
    print(f"Initial pldf shape: {pldf.shape}")
    minipldf = pldf.select(merge_cols).unique()
    print(f"Unique entries in minipldf: {minipldf.shape}")
    
    minipldf = minipldf.filter(~pl.col(idcol).is_in([None, "", "nan"]))
    print(f"After filtering nulls: {minipldf.shape}")
    
    minipldf = minipldf.with_columns(pl.lit(None).alias(seqcol))
    
    seqs = []
    seq_ids = []
    
    # Get actual sequence count from file
    if seqfile.endswith(".gz"):
        seq_count =int(subprocess.run(f"zgrep -c '>'  {seqfile} ", shell=True,capture_output=True, text=True).stdout.strip())
    else:
        seq_count =int(subprocess.run(f"grep -F '>'  {seqfile} -c ", shell=True,capture_output=True, text=True).stdout.strip())
    # seq_count = 0
    # for _ in parse_fastx_file(seqfile):
    #     seq_count += 1
    print(f"Actual number of sequences in file: {seq_count}")
    
    # Reset file iterator
    index = 0
    for record in parse_fastx_file(seqfile):
        seqs.append(record.seq)
        seq_ids.append(record.id)
        index += 1
        
        # Process chunk when we hit chunk_size or end of file
        if len(seqs) >= chunk_size or index == seq_count:
            print(f"\nProcessing chunk {index}/{seq_count}")
            print(f"Number of sequences in chunk: {len(seqs)}")
            
            chunk_seqs = pl.DataFrame({
                idcol: seq_ids,
                seqcol: seqs
            })
            
            chunk_seqs = chunk_seqs.join(minipldf.select(merge_cols), on=idcol, how="inner") #this join get's the info columns (start, end, strand) if needed, only for the entires in this chunk that are in the minipldf.
            
            if trim_to_region:
                print("Trimming sequences")
                # print(chunk_seqs.columns)
                chunk_seqs = chunk_seqs.with_columns(
                    pl.struct(pl.col(seqcol), pl.col(start_col), pl.col(end_col))
                    .map_elements(lambda x: str(x[seqcol][x[start_col]:x[end_col]]) if x[seqcol] is not None else None, return_dtype=pl.Utf8)
                    .alias(seqcol)
                )
            
            if reverse_by_strand_col:
                print("Reversing sequences")
                # print(chunk_seqs.columns)
                chunk_seqs = chunk_seqs.with_columns(
                    pl.when(pl.col(strand_col))
                    .then(pl.col(seqcol).map_elements(lambda x: reverse_complement(x) if x is not None else None, return_dtype=pl.Utf8))
                    .otherwise(pl.col(seqcol))
                    .alias(seqcol)
                )
            
            print("Joining with nascent df")
            minipldf = minipldf.join(chunk_seqs, on=merge_cols, how="left")
            minipldf = minipldf.with_columns(
                pl.coalesce([pl.col(seqcol), pl.col(f"{seqcol}_right")]).alias(seqcol)
            ).drop(f"{seqcol}_right")
            
            print(f"Null count in seqcol after chunk: {minipldf[seqcol].null_count()}")
            
            seqs = []
            seq_ids = []
            # get count for remaining nulls, if zero, break - should be useful when fetching just a few sequences from a large file, at least if the needed seqs are closer to the start of the input fasta.
            if minipldf[seqcol].null_count() == 0:
                break
    
    print("\nFinal merge with original df")
    pldf = pldf.join(minipldf, on=merge_cols, how="left")
    print(f"Final null count in seqcol: {pldf[seqcol].null_count()}")
    
    return pldf


def test_alignment_from_faidx(spacers_file, contigs_file, alignment, **kwargs):
    """Test if an alignment matches between spacer and contig files using pyfastx. 
    alignment is a single row extracted from a pl.df results  (spacer_id, contig_id, spacer_length, start, end, strand, mismatches, tool)
    """
    
    contig_seq = get_seq_from_fastx(contigs_file, alignment['contig_id'][0])
    # Get the sequences
    spacer_seq = get_seq_from_fastx(spacers_file, alignment['spacer_id'][0])
    # Extract the aligned region from contig
    aligned_region = contig_seq[alignment['start'].item():alignment['end'].item()]
        
    # Compare sequences accounting for strand
    if alignment['strand']:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)
    
    alignment = ps.nw_stats_scan(spacer_seq,aligned_region,open=10,extend=10,matrix=ps.nuc44)
    mismatches = len(spacer_seq) - (alignment.matches)
    if mismatches != alignment['mismatches']:
        # print(f"Mismatch between pyfastx and alignment: {mismatches} != {alignment['mismatches']}")
        return False
    return True

def prettify_alignment(spacer_seq, contig_seq, strand=False, start=None, end=None, return_ref_region=False,return_query_region=False,return_comp_str=False, gap_cost=10,extend_cost=5,cost_matrix=ps.nuc44, **kwargs):
    """given a spacer and contig sequence and coordinates, return | for a match and . for a mismatch."""    
    if start is not None and end is not None:
        aligned_region = contig_seq[start:end]
    else:
        aligned_region = contig_seq
        
    if strand:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)
    # alignment = ps.nw_stats_scan(spacer_seq,aligned_region,open=10,extend=10,matrix=ps.nuc44)
    # alignment_string = ""
    # for i in range(len(spacer_seq)):
    #     if spacer_seq[i] == aligned_region[i]:
    #         alignment_string += "|"
    #     else:
    #         alignment_string += "."
    # return f"{spacer_seq}\n{alignment_string}\n{aligned_region}"
    ali=ps.nw_trace(spacer_seq,aligned_region,open=gap_cost,extend=extend_cost,matrix=cost_matrix).traceback
    if return_ref_region:
        return ali.ref
    if return_query_region:
        return ali.query
    if return_comp_str:
        return ali.comp
    return f"{ali.query}\n{ali.comp}\n{ali.ref}"

def test_alignment(spacer_seq, contig_seq, strand=False, start=None, end=None,gap_cost=10,extend_cost=5,  gaps_as_mismatch=True, **kwargs):
    if start is not None and end is not None:
        aligned_region = contig_seq[start:end]  
    else:
        aligned_region = contig_seq
    if strand:  # True means reverse strand
        aligned_region = reverse_complement(aligned_region)
        
    if gaps_as_mismatch:
        alignment = ps.nw_trace(spacer_seq,aligned_region,open=gap_cost,extend=extend_cost,matrix=ps.nuc44)
        mismatches = alignment.traceback.comp.count(" ") + alignment.traceback.comp.count("-") +alignment.traceback.comp.count(".")
        # mismatches = len(spacer_seq) - (alignment.matches)
    else:
        alignment = ps.nw_stats_scan(spacer_seq,aligned_region,open=gap_cost,extend=extend_cost,matrix=ps.nuc44)
        mismatches = len(spacer_seq) - (alignment.matches)

    return mismatches



def test_row(row, ignore_region_strands=False):
    if ignore_region_strands:
        return test_alignment(
            row['spacer_seq'],
            row['contig_seq'], 
            None,
            None,
            None,
        gaps_as_mismatch=True
    )
    else:
        return test_alignment(
            row['spacer_seq'],
            row['contig_seq'], 
            row['strand'],
            row['start'],
            row['end'],
            gaps_as_mismatch=True
        )
    
def test_alignment_polars(results : pl.DataFrame, ignore_region_strands=False, **kwargs):
    """Test if an alignment matches between a spacer and contig.
    Input is a pl.df with columns: spacer_id, contig_id, spacer_length, start, end, strand, mismatches, tool, spacer_seq, contig_seq
    """
    if "spacer_seq" not in results.columns or "contig_seq" not in results.columns:
        raise ValueError("spacer_seq and contig_seq columns are required")

    results = results.with_columns(
            pl.struct(pl.col("spacer_seq"), pl.col("contig_seq"), pl.col("strand"), pl.col("start"), pl.col("end")).map_elements(lambda x: test_row(x, ignore_region_strands=ignore_region_strands), return_dtype=pl.Int32()).alias("alignment_test")
    )
    
    return results

def read_results(tools, max_mismatches=5, spacer_lendf=None, sam_version="1",ref_file=None, threads=4):
    results_df = pl.DataFrame(schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32, "tool": pl.Utf8})
    for tool in tools.values():
        try:
            print(f"\nReading results for {tool['name']}...")
            parse_function = tool.get('parse_function')
            #if the file exists, get the size of it 
            if os.path.exists(tool['output_file']):
                file_size = os.path.getsize(tool['output_file'])
                print(f"File size: {file_size/1024/1024:.2f} MB")
                if file_size == 0:
                    print(f"File {tool['output_file']} is empty, skipping")
                    continue
            else:
                print(f"File {tool['output_file']} does not exist, skipping")
                continue
            print(" output file ", tool['output_file'], "exists")
            parsed_results = parse_function(tool['output_file'], max_mismatches=max_mismatches) if spacer_lendf is None else parse_function(tool['output_file'], max_mismatches=max_mismatches, spacer_lendf=spacer_lendf, threads=threads, ref_file=ref_file)
            parsed_results = parsed_results.with_columns(pl.lit(tool['name']).alias('tool'))
            tool_df = pl.DataFrame(data=parsed_results, orient="row",
                                   schema={"spacer_id": pl.Utf8, "contig_id": pl.Utf8, "spacer_length": pl.UInt32, "strand": pl.Boolean, "start": pl.UInt32, "end": pl.UInt32, "mismatches": pl.UInt32, "tool": pl.Utf8})
            # print(parsed_results)
            results_df = results_df.vstack(tool_df)
        except Exception as e:
            print(f"Failed to read results for {tool['name']}: {e}")
    return results_df

def add_sqlines_sam(sam_file, ref_file=None):
    first_line, second_line = os.popen(f"head -n 2 {sam_file}").read().splitlines()
    if second_line.startswith("@SQ"):
        print("sam file looks right")
        return sam_file
    else:
        print("sam file looks wrong, adding SQ lines")
        print("first_line of sam file", first_line)
        print("second_line of sam file", second_line)
        
    if ref_file is not None:
        print("Reference file provided, using it to add SQ lines")
        tmp = pysam.FastaFile(ref_file)
        pysam.AlignmentFile(filename = sam_file+".sqlines", mode="wh", reference_names=tmp.references, reference_lengths=tmp.lengths)
    else:
        print("No reference file provided, using the smallest sam file present in the output file directory")
        present_sam_files = glob.glob(os.path.join(os.path.dirname(sam_file), "*.sam"))
        present_sam_files.remove(os.path.normpath(sam_file))
        # drop mummer4 and minimap2 sam files
        present_sam_files = [thing for thing in present_sam_files if "mummer" not in thing and "minimap" not in thing] # these don't print SQ lines when batching.
        present_sam_files = [thing for thing in present_sam_files if os.path.getsize(thing) > 1]
        if len(present_sam_files) == 0:
            print("No non-empty sam files found, returning None")
            return None
        smallest_sam_file = min(present_sam_files, key=os.path.getsize)
        # write header lines to temp new file
        with open(sam_file+".sqlines", 'w') as sqlines_file:
            has_sq_lines = os.popen(f"head -n 2 {smallest_sam_file}").read().splitlines()[1].startswith("@SQ") 
            if not has_sq_lines:
                smallest_sam_file = min(present_sam_files.remove(smallest_sam_file), key=os.path.getsize) # one try again
            with open(smallest_sam_file, 'r') as f2:
                for line in f2:
                    if line.startswith('@SQ'):
                        sqlines_file.write(line)
        
    with open(sam_file+"new", 'w') as new_sam_file:
        new_sam_file.write(first_line)
        new_sam_file.write("\n")
    # subprocess.run(f"tail -n+3 {sam_file} >> {temp_sam_file}", shell=True) # faster than reading line by line                     
        with open(sam_file+".sqlines", 'r') as sqlines_file: # add the SQ lines to the new sam file
            for line in sqlines_file:
                new_sam_file.write(line)
        new_sam_file.write(second_line.replace("@PG ID:", "@PG\tID:").replace(" PN:", "\tPN:").replace(" VN:", "\tVN:").replace(" CL:", "\tCL:")) # add the PG line to the new sam file
        new_sam_file.write("\n")
        with open(sam_file, 'r') as old_sam_file:
            for line in old_sam_file:
                if not line.startswith('@'): # keep only the alignments, mummer 4.01 seems to output multiple PG lines.
                    new_sam_file.write(line)
    os.rename(sam_file+"new", sam_file)
    os.remove(sam_file+".sqlines")
    return sam_file

def prefilter_sam_with_sambamba(input_sam, output_sam,max_mismatches=5, threads=1):
    """Prefilter SAM file using sambamba to remove unmapped reads and reads with gaps"""
    # Filter expression:
    # - not unmapped: keep only mapped reads
    # - [NM] <= max_mismatches: keep reads with max_mismatches or fewer mismatches
    # - cigar !~ /[ID]/: exclude reads with insertions or deletions in CIGAR string
    
    filter_expression = f'"not (unmapped) and [NM] <= {max_mismatches}"'
    if not os.path.exists(input_sam):
        print(f"File {input_sam} does not exist, skipping")
        return None
    if os.path.getsize(input_sam) == 0:
        print(f"File {input_sam} is empty, skipping")
        return None
    
    input_sam = add_sqlines_sam(input_sam)

    command = [
        "sambamba", "view",
        "-F", filter_expression,
        "--show-progress",
        "-t", str(threads),
        "-S",  # SAM input
        "-f", "sam",  # SAM output
        "-h",
        "-o", output_sam,
        input_sam
    ]
    subprocess.run(" ".join(command), check=True,shell=True)
    return output_sam


def create_comparison_matrix(tools_results, n_mismatches):
    
    # Filter for specific number of mismatches
    tmp = tools_results.filter(pl.col('mismatches') == n_mismatches)
    
    # Get unique tools
    tools = tmp.select('tool').unique()
    tools_list = tools.get_column('tool').to_list()
    
    # Create empty matrix dataframe
    matrix_data = []
    
    for tool1 in tools_list:
        row_data = []
        # Get unique pairs for tool1
        tool1_pairs = tmp.filter(pl.col('tool') == tool1).select(["contig_id", "spacer_id"]).unique()
        
        for tool2 in tools_list:
            if tool1 == tool2:
                # Diagonal will show total number of pairs for the tool
                row_data.append(tool1_pairs.height)
            else:
                # Get unique pairs for tool2
                tool2_pairs = tmp.filter(pl.col('tool') == tool2).select(["contig_id", "spacer_id"]).unique()
                
                # Count pairs in tool1 that are not in tool2
                diff_count = (tool1_pairs
                    .join(tool2_pairs, 
                          on=['contig_id', 'spacer_id'], 
                          how='anti')
                    .height)
                row_data.append(diff_count)
        
        matrix_data.append(row_data)
    
    # Create matrix dataframe
    matrix_df = pl.DataFrame(
        matrix_data,
        schema=tools_list,
    ).with_columns(pl.Series(name='tool', values=tools_list))
    
    return matrix_df


def compare_aligner(ground_truth: pl.DataFrame, aligner_results: pl.DataFrame, return_classification=False, max_mismatches=5, soft_false_positive_threshold=0):
    # Convert string strands to boolean if needed
    if ground_truth["strand"].dtype == pl.Utf8:
        ground_truth = ground_truth.with_columns(
            (pl.col("strand") == "-").alias("strand")
        )
    if aligner_results["strand"].dtype == pl.Utf8:
        aligner_results = aligner_results.with_columns(
            (pl.col("strand") == "-").alias("strand")
        )
        
    aligner_results = aligner_results.with_row_index()
    
    aligner_results = aligner_results.filter(pl.col("mismatches") <= max_mismatches)
    
    true_positives = ground_truth.join(
        aligner_results,
        on=["spacer_id", "contig_id", "strand", "start", "end"],
        how="inner",
        suffix="_aligner")
    
    true_positives = true_positives.unique(subset=["spacer_id", "contig_id", "strand", "start", "end", "mismatches", "tool", "mismatches_aligner"], keep="last")
    # Find false positives (predictions not in ground truth). Trying to see if the aligner found any insertions (i.e. pairs of spacer_id and contig_id) that are not in the ground truth
    false_positives = aligner_results.join(
        ground_truth, 
        on=["spacer_id", "contig_id", "strand", "start", "end"], 
        how="anti")
    
    # get false positives:
    # try to see if see if these are soft false positives or totally off  (i.e. if the start and end positions are close to the ground truth)
    # these can happen if the aligner used gaps or if the mismatches are terminal causing the alignment to be off by a few bases.
    for row in false_positives.iter_rows(named=True):
        index = row["index"]
        spacer_id = row["spacer_id"]
        contig_id = row["contig_id"]
        strand = row["strand"]
        start = row["start"]
        end = row["end"]
        ground_truth_same_pair = ground_truth.filter(pl.col("spacer_id") == spacer_id, pl.col("contig_id") == contig_id)
        if len(ground_truth_same_pair) == 0:
            print(f"No ground truth for {spacer_id} {contig_id}")
            continue
        if ground_truth_same_pair.height > 0:
            diff_start = abs(start - ground_truth_same_pair["start"].to_list()[0])
            diff_end = abs(end - ground_truth_same_pair["end"].to_list()[0])
            if diff_start <= soft_false_positive_threshold and diff_end <= soft_false_positive_threshold:
                # get the ground truth row's mismatches
                ground_truth_mismatches = ground_truth_same_pair["mismatches"].to_list()[0]
                row["mismatches_aligner"] = row["mismatches"]
                row["mismatches"] = ground_truth_mismatches
                row = order_columns_to_match(pl.from_dict(row), true_positives)
                # remove from false positives, and add to true positives
                false_positives = false_positives.filter(pl.col("index") != index)
                true_positives = vstack_easy(true_positives, row)
                # print(f"Soft false positive: {spacer_id} {contig_id} diff start: {diff_start} diff end: {diff_end}")    
                
    # Calculate metrics
    total_true_positives = true_positives.height
    total_false_positives = false_positives.height
    total_ground_truth = ground_truth.height
    
    # Calculate precision and recall
    Precision = total_true_positives / (total_true_positives + total_false_positives) if (total_true_positives + total_false_positives) > 0 else 0
    Recall = total_true_positives / total_ground_truth if total_ground_truth > 0 else 0
    
    # Calculate F1 score
    F1 = 2 * (Precision * Recall) / (Precision + Recall) if (Precision + Recall) > 0 else 0
    
    if return_classification:
        # Add a classification column to the aligner results
        aligner_results = aligner_results.with_columns(
            pl.when(pl.col("index").is_in(true_positives["index"]))
            .then(pl.lit("true_positive"))
            .otherwise(pl.lit("false_positive"))
            .alias("classification")
        )
        return aligner_results
    else:
        return {
        "tool": aligner_results['tool'].unique()[0],
        "true_positives": total_true_positives,
        "false_positives": total_false_positives,
        "ground_truth_total": total_ground_truth,
        "precision": Precision,
        "recall": Recall,
        "f1_score": F1
        }
        
def compare_results(tools, ground_truth, tools_results,soft_false_positive_threshold=0):
    comparison_results = []

    for tool in tools.values():
        try:
            result = tools_results.filter(pl.col("tool") == tool['name'])
            comp=compare_aligner(ground_truth, result,soft_false_positive_threshold=soft_false_positive_threshold)
            # temp2=compare_aligner(ground_truth, result, return_classification=True)
            comparison_results.append(comp)
        except Exception as e:
            print(f"Failed to compare results for {tool['name']}: {e}")

    # Create final results dataframe
    results_df = pl.DataFrame(comparison_results)
    # print(results_df)
    return results_df

def get_parse_function(func_name):
    return globals()[func_name]

class DebugArgs:
    def __init__(self):
        self.n_mismatch_range = [0, 0]
        self.contig_length_range = [1200, 12000]
        self.spacer_length_range = [19, 60]
        self.sample_size_contigs = 154
        self.sample_size_spacers = 7
        self.insertion_range = [5, 10]
        self.threads = 1
        self.contigs = None
        self.spacers = None
        self.prop_rc = 0.5  # Add default value for prop_rc
        self.max_mismatches = 5
        self.max_runs = 1
        
        
        
# # interactive debug:
# args = DebugArgs()

# os.makedirs("results", exist_ok=True) # make sure parent results directory exists

# contigs, spacers, ground_truth = simulate_data(
#     contigs=args.contigs,
#     spacers=args.spacers,
#     contig_length_range=args.contig_length_range,
#     spacer_length_range=args.spacer_length_range,
#     n_mismatch_range=args.n_mismatch_range,
#     sample_size_contigs=args.sample_size_contigs,
#     sample_size_spacers=args.sample_size_spacers,
#     insertion_range=args.insertion_range,
#     prop_rc=args.prop_rc
#     )
# args.sample_size_contigs = len(contigs)
# args.sample_size_spacers = len(spacers)

# results_dir = f"results/run_t_{args.threads}_nc_{args.sample_size_contigs}_ns_{args.sample_size_spacers}_ir_{args.insertion_range[0]}_{args.insertion_range[1]}_lm_{args.n_mismatch_range[0]}_{args.n_mismatch_range[1]}_prc_{args.prop_rc}"
# args.results_dir = results_dir # add to the args object for simplicity

# os.makedirs(results_dir, exist_ok=True)
# os.makedirs(f"{results_dir}/simulated_data", exist_ok=True)
# os.makedirs(f"{results_dir}/raw_outputs", exist_ok=True)
# os.makedirs(f"{results_dir}/bash_scripts", exist_ok=True)
# write_fasta(contigs, f"{results_dir}/simulated_data/simulated_contigs.fa")
# write_fasta(spacers, f"{results_dir}/simulated_data/simulated_spacers.fa")
# ground_truth.write_csv(f"{results_dir}/simulated_data/ground_truth.tsv", separator="\t")
# spacer_lendf = pl.DataFrame({"spacer_id": spacers.keys(), "length": [len(seq) for seq in spacers.values()]})
# tools = populate_tools(args, spacer_lendf=spacer_lendf)

# # Run tools and get results
# run_tools(tools, results_dir, max_runs=args.max_runs)
# tools_results = read_results(tools, max_mismatches=args.max_mismatches, spacer_lendf=spacer_lendf)
# tools_results.write_csv(f"{results_dir}/tools_results.tsv", separator="\t")

# performance_results = compare_results(tools, ground_truth, tools_results,soft_false_positive_threshold=0)

# def convert_strand_for_display(df):
#     """Convert boolean strand to +/- for display purposes"""
#     return df.with_columns(
#         pl.when(pl.col("strand"))
#         .then(pl.lit("-"))
#         .otherwise(pl.lit("+"))
#         .alias("strand")
#     )
