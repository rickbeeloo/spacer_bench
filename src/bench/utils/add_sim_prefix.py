#!/usr/bin/env python

import os
import polars as pl
import argparse
import subprocess
from pathlib import Path
from .functions import generate_simulation_id, read_fasta, write_fasta

def has_prefix(file_path, sim_id):
    """Check if file already has simulation ID prefix."""
    if not os.path.exists(file_path):
        print(f"File does not exist: {file_path}")
        return False
        
    # For FASTA files
    if file_path.endswith('.fa'):
        with open(file_path) as f:
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                has_prefix = first_line[1:].startswith(sim_id)
                print(f"Checking FASTA file {file_path}: {'has prefix' if has_prefix else 'no prefix'}")
                return has_prefix
    
    # For TSV files
    elif file_path.endswith('.tsv'):
        try:
            df = pl.scan_csv(file_path, separator="\t").collect()
            if 'spacer_id' in df.columns:
                has_prefix = df['spacer_id'][0].startswith(sim_id)
                print(f"Checking TSV file {file_path} (spacer_id): {'has prefix' if has_prefix else 'no prefix'}")
                return has_prefix
            elif 'contig_id' in df.columns:
                has_prefix = df['contig_id'][0].startswith(sim_id)
                print(f"Checking TSV file {file_path} (contig_id): {'has prefix' if has_prefix else 'no prefix'}")
                return has_prefix
        except Exception as e:
            print(f"Error checking TSV file {file_path}: {e}")
            return False
    
    # For SAM files
    elif file_path.endswith('.sam'):
        try:
            with open(file_path) as f:
                for line in f:
                    if line.startswith('@SQ'):
                        has_prefix = f"SN:{sim_id}_" in line
                        print(f"Checking SAM file {file_path}: {'has prefix' if has_prefix else 'no prefix'}")
                        return has_prefix
                    elif not line.startswith('@'):
                        # Check first alignment line
                        fields = line.split('\t')
                        if len(fields) >= 3:
                            has_prefix = fields[0].startswith(sim_id) and fields[2].startswith(sim_id)
                            print(f"Checking SAM file {file_path}: {'has prefix' if has_prefix else 'no prefix'}")
                            return has_prefix
        except Exception as e:
            print(f"Error checking SAM file {file_path}: {e}")
            return False
            
    return False

def update_raw_output(file_path, sim_id):
    """Update IDs in raw output files using sed."""
    if not os.path.exists(file_path):
        print(f"File does not exist: {file_path}")
        return
        
    file_ext = os.path.splitext(file_path)[1]
    
    # Create backup
    backup_file = f"{file_path}.bak"
    os.system(f"cp {file_path} {backup_file}")
    
    try:
        if file_ext == '.sam':
            print(f"Updating SAM file: {file_path}")
            # Update sequence headers in @SQ lines
            subprocess.run(f"sed -i '/^@SQ/s/\\(SN:\\)\\(contig_[0-9]\\+\\)/\\1{sim_id}_\\2/' {file_path}", shell=True)
            # Update query names in alignment lines (first field)
            subprocess.run(f"sed -i '/^[^@]/s/^\\([^\\t]*\\)/{sim_id}_\\1/' {file_path}", shell=True)
            # Update reference names in alignment lines (contig IDs)
            subprocess.run(f"sed -i 's/\\(\\t\\)\\(contig_[0-9]\\+\\)/\\1{sim_id}_\\2/g' {file_path}", shell=True)
            # Clean up any double prefixes
            subprocess.run(f"sed -i 's/{sim_id}_{sim_id}_/{sim_id}_/g' {file_path}", shell=True)
            
        elif file_ext in ['.tsv', '.blast6']:
            print(f"Updating {file_ext} file: {file_path}")
            # Update contig IDs
            subprocess.run(f"sed -i 's/\\(^\\|\\t\\)\\(contig_[0-9]\\+\\)/\\1{sim_id}_\\2/g' {file_path}", shell=True)
            # Update spacer IDs
            subprocess.run(f"sed -i 's/\\(^\\|\\t\\)\\(spacer_[0-9]\\+\\)/\\1{sim_id}_\\2/g' {file_path}", shell=True)
            # Clean up any double prefixes
            subprocess.run(f"sed -i 's/{sim_id}_{sim_id}_/{sim_id}_/g' {file_path}", shell=True)
            
        print(f"Successfully updated: {file_path}")
        os.remove(backup_file)
        
    except Exception as e:
        print(f"Error updating {file_path}: {e}")
        print(f"Restoring from backup...")
        os.system(f"mv {backup_file} {file_path}")

def clean_prefixes(file_path):
    """Remove any existing prefixes (hash_) from file contents."""
    if not os.path.exists(file_path):
        return
        
    print(f"Cleaning prefixes from: {file_path}")
    backup_file = f"{file_path}.bak"
    os.system(f"cp {file_path} {backup_file}")
    
    try:
        file_ext = os.path.splitext(file_path)[1]
        
        if file_ext == '.fa':
            # Clean FASTA headers
            subprocess.run(f"sed -i 's/^>[0-9a-f]{{8}}_/>/' {file_path}", shell=True)
            
        elif file_ext == '.tsv':
            # Clean prefixes from contig and spacer IDs
            subprocess.run(f"sed -i 's/\\(^\\|\\t\\)[0-9a-f]{{8}}_\\(contig_[0-9]\\+\\|spacer_[0-9]\\+\\)/\\1\\2/g' {file_path}", shell=True)
            
        elif file_ext == '.sam':
            # Clean SAM headers
            subprocess.run(f"sed -i '/^@SQ/s/SN:[0-9a-f]{{8}}_\\(contig_[0-9]\\+\\)/SN:\\1/' {file_path}", shell=True)
            # Clean query names (first field)
            subprocess.run(f"sed -i '/^[^@]/s/^[0-9a-f]{{8}}_\\(spacer_[0-9]\\+\\)/\\1/' {file_path}", shell=True)
            # Clean contig IDs
            subprocess.run(f"sed -i 's/\\t[0-9a-f]{{8}}_\\(contig_[0-9]\\+\\)/\\t\\1/g' {file_path}", shell=True)
            
        elif file_ext == '.blast6':
            # Clean prefixes from contig and spacer IDs
            subprocess.run(f"sed -i 's/\\(^\\|\\t\\)[0-9a-f]{{8}}_\\(contig_[0-9]\\+\\|spacer_[0-9]\\+\\)/\\1\\2/g' {file_path}", shell=True)
            
        print(f"Successfully cleaned prefixes from: {file_path}")
        os.remove(backup_file)
        
    except Exception as e:
        print(f"Error cleaning prefixes from {file_path}: {e}")
        print(f"Restoring from backup...")
        os.system(f"mv {backup_file} {file_path}")

def clean_directory(results_dir):
    """Clean all prefixes from files in a results directory."""
    print(f"\nCleaning prefixes from directory: {results_dir}")
    
    # Clean FASTA files
    contigs_file = f"{results_dir}/simulated_data/simulated_contigs.fa"
    spacers_file = f"{results_dir}/simulated_data/simulated_spacers.fa"
    
    if os.path.exists(contigs_file):
        clean_prefixes(contigs_file)
    if os.path.exists(spacers_file):
        clean_prefixes(spacers_file)
    
    # Clean ground truth and tools results
    ground_truth_file = f"{results_dir}/simulated_data/ground_truth.tsv"
    tools_results_file = f"{results_dir}/tools_results.tsv"
    
    if os.path.exists(ground_truth_file):
        clean_prefixes(ground_truth_file)
    if os.path.exists(tools_results_file):
        clean_prefixes(tools_results_file)
    
    # Clean raw output files
    raw_outputs_dir = f"{results_dir}/raw_outputs"
    if os.path.exists(raw_outputs_dir):
        print(f"\nCleaning raw output files in: {raw_outputs_dir}")
        for file_name in os.listdir(raw_outputs_dir):
            if file_name.endswith(('.sam', '.tsv', '.blast6')):
                file_path = os.path.join(raw_outputs_dir, file_name)
                clean_prefixes(file_path)

def add_prefix_to_files(results_dir):
    """Add simulation ID prefix to all files in a simulation results directory."""
    # Parse parameters from directory name
    dir_name = os.path.basename(results_dir)
    params = parse_dir_name(dir_name)
    sim_id = generate_simulation_id(params)
    print(f"\nProcessing directory: {dir_name}")
    print(f"Generated simulation ID: {sim_id}")
    
    # First clean any existing prefixes
    clean_directory(results_dir)
    
    # Update FASTA files
    contigs_file = f"{results_dir}/simulated_data/simulated_contigs.fa"
    spacers_file = f"{results_dir}/simulated_data/simulated_spacers.fa"
    
    if os.path.exists(contigs_file):
        if not has_prefix(contigs_file, sim_id):
            print(f"Updating contigs file: {contigs_file}")
            contigs = read_fasta(contigs_file)
            contigs = {f"{sim_id}_{k}": v for k, v in contigs.items()}
            write_fasta(contigs, contigs_file)
            print(f"Successfully updated contigs file")
        else:
            print(f"Skipping contigs file (already has prefix): {contigs_file}")
    
    if os.path.exists(spacers_file):
        if not has_prefix(spacers_file, sim_id):
            print(f"Updating spacers file: {spacers_file}")
            spacers = read_fasta(spacers_file)
            spacers = {f"{sim_id}_{k}": v for k, v in spacers.items()}
            write_fasta(spacers, spacers_file)
            print(f"Successfully updated spacers file")
        else:
            print(f"Skipping spacers file (already has prefix): {spacers_file}")
    
    # Update ground truth file
    ground_truth_file = f"{results_dir}/simulated_data/ground_truth.tsv"
    if os.path.exists(ground_truth_file):
        if not has_prefix(ground_truth_file, sim_id):
            print(f"Updating ground truth file: {ground_truth_file}")
            ground_truth = pl.read_csv(ground_truth_file, separator="\t")
            ground_truth = ground_truth.with_columns([
                (pl.lit(sim_id) + "_" + pl.col("spacer_id")).alias("spacer_id"),
                (pl.lit(sim_id) + "_" + pl.col("contig_id")).alias("contig_id")
            ])
            ground_truth.write_csv(ground_truth_file, separator="\t")
            print(f"Successfully updated ground truth file")
        else:
            print(f"Skipping ground truth file (already has prefix): {ground_truth_file}")
    
    # Update tools results file
    tools_results_file = f"{results_dir}/tools_results.tsv"
    if os.path.exists(tools_results_file):
        if not has_prefix(tools_results_file, sim_id):
            print(f"Updating tools results file: {tools_results_file}")
            tools_results = pl.read_csv(tools_results_file, separator="\t")
            tools_results = tools_results.with_columns([
                (pl.lit(sim_id) + "_" + pl.col("spacer_id")).alias("spacer_id"),
                (pl.lit(sim_id) + "_" + pl.col("contig_id")).alias("contig_id")
            ])
            tools_results.write_csv(tools_results_file, separator="\t")
            print(f"Successfully updated tools results file")
        else:
            print(f"Skipping tools results file (already has prefix): {tools_results_file}")
    
    # Update raw output files
    raw_outputs_dir = f"{results_dir}/raw_outputs"
    if os.path.exists(raw_outputs_dir):
        print(f"\nChecking raw output files in: {raw_outputs_dir}")
        for file_name in os.listdir(raw_outputs_dir):
            if file_name.endswith(('.sam', '.tsv', '.blast6')):
                file_path = os.path.join(raw_outputs_dir, file_name)
                if not has_prefix(file_path, sim_id):
                    update_raw_output(file_path, sim_id)
                else:
                    print(f"Skipping raw output file (already has prefix): {file_path}")

def parse_dir_name(dir_name):
    """Extract simulation parameters from directory name."""
    parts = dir_name.split('_')
    params = {}
    
    for i, part in enumerate(parts):
        if part == 't':
            params['threads'] = int(parts[i+1])
        elif part == 'nc':
            params['sample_size_contigs'] = int(parts[i+1])
        elif part == 'ns':
            params['sample_size_spacers'] = int(parts[i+1])
        elif part == 'ir':
            params['insertion_range'] = [int(parts[i+1]), int(parts[i+2])]
        elif part == 'lm':
            params['n_mismatch_range'] = [int(parts[i+1]), int(parts[i+2])]
        elif part == 'prc':
            params['prop_rc'] = float(parts[i+1])

    return params

def main():
    parser = argparse.ArgumentParser(description='Add simulation ID prefix to files in results directory')
    parser.add_argument('results_dir', type=str, help='Path to results directory')
    args = parser.parse_args()

    if not os.path.exists(args.results_dir):
        print(f"Results directory {args.results_dir} does not exist")
        return
    
    for dir_name in os.listdir(args.results_dir):
        if dir_name.startswith("run_"):
            results_dir = os.path.join(args.results_dir, dir_name)
            add_prefix_to_files(results_dir)

if __name__ == "__main__":
    main()
    
    
    
# bash script to clear malformatted simulation ID prefix to all files in a results directory
"""
prefixes=("6b1148fc" "01da6dbd" "f25c897e" "6ef5b781" "4b7b4f05" "adb3760b" "6951d55c" "936dac80" "a238f23b" "493dc49e" "6a19eb7c" "59e8522e" "6b1148fc" "578a5025" "37d6023d" "9be1d7f4" "93cc9b21")


# function to return true if  a file is  relevant  (.txt, .tsv, .sam, .tab)
is_relevant_file() {
  local file="$1"
  if [[ "$file" == *.txt || "$file" == *.tsv || "$file" == *.sam || "$file" == *.tab || "$file" == *.fa ]]; then
    return true
  else
    return false
  fi
}


for run_dir in ./results/simulated/run_t_*/; do
  echo "run_dir $run_dir"
  for file in "$run_dir"/*; do
  if [ -f "$file" ]; then
    echo "$file"
    for prefix in "${prefixes[@]}"; do
      echo "prefix $prefix"
      if is_relevant_file "$file"; then
        echo "file $file"
        sed -i "s|${prefix}_||g" "$file"
      fi
    done
  fi
  
  if [ -d "$file" ]; then
    echo "$file"
    for subfile in "$file"/*; do
      if [ -f "$subfile" ]; then
        for prefix in "${prefixes[@]}"; do
          echo "prefix $prefix"
          if is_relevant_file "$subfile"; then
            echo "subfile $subfile"
            sed -i "s|${prefix}_||g" "$subfile"
          fi
        done
      fi
    done
  fi
done
done
"""