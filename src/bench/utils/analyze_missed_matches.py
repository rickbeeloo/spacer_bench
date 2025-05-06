import polars as pl
import argparse
from .functions import populate_pldf_withseqs, prettify_alignment
import os

def analyze_missed_matches(run_dir, max_entries=10, output_file=None):
    # Redirect stdout to file if specified
    if output_file:
        import sys
        sys.stdout = open(output_file, 'w')

    # Read ground truth and results
    ground_truth = pl.read_csv(f"{run_dir}/simulated_data/ground_truth.tsv", separator="\t")
    tools_results = pl.read_csv(f"{run_dir}/tools_results.tsv", separator="\t")
    
    # Debug: Print sample IDs and file paths
    print("\nDebug Information:")
    print(f"Sample spacer_ids from ground truth: {ground_truth['spacer_id'].head(5)}")
    spacers_path = f"{run_dir}/simulated_data/simulated_spacers.fa"
    contigs_path = f"{run_dir}/simulated_data/simulated_contigs.fa"
    
    # Read and check FASTA headers
    def read_fasta_headers(filepath):
        headers = []
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    headers.append(line.strip()[1:])  # Remove '>' and whitespace
        return headers

    print("\nChecking FASTA headers:")
    spacer_headers = read_fasta_headers(spacers_path)
    contig_headers = read_fasta_headers(contigs_path)
    print(f"First few spacer headers: {spacer_headers[:3]}")
    print(f"First few contig headers: {contig_headers[:3]}")
    
    # Filter ground truth for perfect matches
    perfect_matches = ground_truth.filter(pl.col("mismatches") == 0)
    
    # Get unique tools
    tools = tools_results["tool"].unique().to_list()
    
    # For each tool, find missed perfect matches
    for tool in tools:
        print(f"\nAnalyzing missed matches for {tool}:")
        print("----------------------------------------")
        
        # Get tool's results
        tool_results = tools_results.filter(pl.col("tool") == tool)
        
        # Find perfect matches that this tool missed
        missed_matches = perfect_matches.join(
            tool_results,
            on=["spacer_id", "contig_id", "start", "end", "strand"],
            how="anti"
        ).head(max_entries)
        
        if missed_matches.height == 0:
            print("No perfect matches were missed!")
            continue
            
        print(f"Found {missed_matches.height} missed perfect matches")
        print("\nSample of missed matches:")
        print(missed_matches.head())
        
        try:
            # Manual sequence lookup to debug the issue
            def get_sequence_from_fasta(fasta_path, target_id):
                with open(fasta_path, 'r') as f:
                    current_id = None
                    current_seq = []
                    for line in f:
                        if line.startswith('>'):
                            if current_id == target_id and current_seq:
                                return ''.join(current_seq)
                            current_id = line.strip()[1:]  # Remove '>' and whitespace
                            current_seq = []
                        else:
                            if current_id == target_id:
                                current_seq.append(line.strip())
                    if current_id == target_id and current_seq:
                        return ''.join(current_seq)
                return None

            # Process each missed match manually
            for row in missed_matches.iter_rows(named=True):
                print(f"\nMissed match:")
                print(f"spacer_id: {row['spacer_id']}")
                print(f"contig_id: {row['contig_id']}")
                print(f"position: {row['start']}-{row['end']}")
                print(f"strand: {row['strand']}")
                
                spacer_seq = get_sequence_from_fasta(spacers_path, row['spacer_id'])
                contig_seq = get_sequence_from_fasta(contigs_path, row['contig_id'])
                
                if spacer_seq and contig_seq:
                    print(f"spacer length: {len(spacer_seq)}")
                    print("\nAlignment:")
                    print(prettify_alignment(
                        spacer_seq,
                        contig_seq,
                        row['strand'],
                        row['start'],
                        row['end']
                    ))
                else:
                    print("WARNING: Could not find sequences in FASTA files")
                    if not spacer_seq:
                        print(f"Missing spacer sequence for {row['spacer_id']}")
                    if not contig_seq:
                        print(f"Missing contig sequence for {row['contig_id']}")
                
                # Look for any hits this tool found nearby
                nearby_hits = tool_results.filter(
                    (pl.col("spacer_id") == row["spacer_id"]) &
                    (pl.col("contig_id") == row["contig_id"])
                )
                
                if nearby_hits.height > 0:
                    print("\nTool found these hits for the same spacer-contig pair:")
                    print(nearby_hits)
                    
        except Exception as e:
            print(f"Error processing sequences: {str(e)}")
            import traceback
            print(traceback.format_exc())
            print("Continuing with next tool...")

    # Reset stdout if it was redirected
    if output_file:
        sys.stdout.close()
        sys.stdout = sys.__stdout__

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze matches missed by tools')
    parser.add_argument('--run_dir', type=str, 
                       default="results/simulated/run_t_10_nc_1000_ns_300_ir_1_200_lm_0_0_prc_0.5",
                       help='Directory containing the run results')
    parser.add_argument('--max_entries', type=int, default=10,
                       help='Maximum number of missed matches to analyze per tool')
    parser.add_argument('--output', type=str,
                       help='Optional file to write output to')
    
    args = parser.parse_args()
    analyze_missed_matches(args.run_dir, args.max_entries, args.output)