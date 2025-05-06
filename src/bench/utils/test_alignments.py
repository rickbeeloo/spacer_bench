import polars as pl
import argparse
from .functions import populate_pldf_withseqs, prettify_alignment, test_alignment_polars

def main(args):
    # Read ground truth and results
    ground_truth = pl.read_csv(f"{args.run_dir}/simulated_data/ground_truth.tsv", separator="\t")
    tools_results = pl.read_csv(f"{args.run_dir}/tools_results.tsv", separator="\t")

    # Take first n entries from ground truth
    test_entries = ground_truth.head(args.num_entries)

    # Get sequences for these entries
    test_entries = populate_pldf_withseqs(
        test_entries,
        f"{args.run_dir}/simulated_data/simulated_spacers.fa",
        idcol="spacer_id",
        seqcol="spacer_seq"
    )

    test_entries = populate_pldf_withseqs(
        test_entries,
        f"{args.run_dir}/simulated_data/simulated_contigs.fa",
        idcol="contig_id",
        seqcol="contig_seq"
    )

    # Print alignments
    print("\nGround Truth Alignments:")
    print("------------------------")
    for row in test_entries.iter_rows(named=True):
        print(f"\nspacer_id: {row['spacer_id']}")
        print(f"contig_id: {row['contig_id']}")
        print(f"position: {row['start']}-{row['end']}")
        print(f"strand: {row['strand']}")
        print(f"expected mismatches: {row['mismatches']}")
        print("\nAlignment:")
        print(prettify_alignment(
            row['spacer_seq'],
            row['contig_seq'],
            row['strand'],
            row['start'],
            row['end']
        ))

    # Get tool results for these entries
    for entry in test_entries.iter_rows(named=True):
        print(f"\nResults for {entry['spacer_id']} in {entry['contig_id']}:")
        print("----------------------------------------")
        
        # Get all tool results for this spacer-contig pair
        tool_matches = tools_results.filter(
            (pl.col("spacer_id") == entry["spacer_id"]) & 
            (pl.col("contig_id") == entry["contig_id"])
        )
        
        if tool_matches.height == 0:
            print("No tool found this match")
            continue
            
        # Add sequences
        tool_matches = populate_pldf_withseqs(
            tool_matches,
            f"{args.run_dir}/simulated_data/simulated_spacers.fa",
            idcol="spacer_id",
            seqcol="spacer_seq"
        )
        
        tool_matches = populate_pldf_withseqs(
            tool_matches,
            f"{args.run_dir}/simulated_data/simulated_contigs.fa",
            idcol="contig_id",
            seqcol="contig_seq"
        )
        
        # Test alignments
        tool_matches = test_alignment_polars(tool_matches)
        
        # Print results
        for row in tool_matches.iter_rows(named=True):
            print(f"\nTool: {row['tool']}")
            print(f"Position: {row['start']}-{row['end']}")
            print(f"Strand: {row['strand']}")
            print(f"Reported mismatches: {row['mismatches']}")
            print(f"Actual mismatches: {row['alignment_test']}")
            print("\nAlignment:")
            print(prettify_alignment(
                row['spacer_seq'],
                row['contig_seq'],
                row['strand'],
                row['start'],
                row['end']
            ))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test and visualize alignments from tools results')
    parser.add_argument('--run_dir', type=str, 
                       default="results/simulated/run_t_10_nc_1000_ns_300_ir_1_200_lm_0_0_prc_0.5",
                       help='Directory containing the run results')
    parser.add_argument('--num_entries', type=int, default=5,
                       help='Number of entries to test from ground truth')
    
    args = parser.parse_args()
    main(args)