import polars as pl
import numpy as np
import matplotlib.pyplot as plt

def create_spacer_counts_with_tools(recalc_only, tools_list, max_mismatches=3):
    """Create a DataFrame with tool performance per spacer.
    
    Args:
        recalc_only (pl.DataFrame): DataFrame with recalculated results
        tools_list (list): List of tools to analyze
        max_mismatches (int): Maximum number of mismatches to consider
    """
    # First get total occurrences per spacer across all tools
    spacer_counts = recalc_only.filter(pl.col('mismatches') <= max_mismatches)\
        .select(["spacer_id", "contig_id"])\
        .unique()\
        .group_by('spacer_id')\
        .agg(pl.count('contig_id').alias('n_occurrences'))

    # Calculate matches per tool and spacer without joining to spacer_counts yet
    tool_matches = recalc_only.filter(pl.col('mismatches') <= max_mismatches)\
        .select(['spacer_id', 'tool', 'contig_id'])\
        .unique()\
        .group_by(['spacer_id', 'tool'])\
        .agg(pl.count('contig_id').alias('tool_matches'))

    # Create a cross join of all spacers with all tools
    all_combinations = spacer_counts.select('spacer_id', 'n_occurrences')\
        .join(
            pl.DataFrame({'tool': tools_list}),
            how='cross'
        )

    # Join the actual matches and calculate fractions
    complete_fractions = all_combinations\
        .join(
            tool_matches,
            on=['spacer_id', 'tool'],
            how='left'
        )\
        .with_columns([
            pl.col('tool_matches').fill_null(0),
            (pl.col('tool_matches') / pl.col('n_occurrences')).alias('fraction')
        ])

    # Pivot to get tools as columns
    spacer_counts_with_tools = complete_fractions\
        .pivot(
            index=['spacer_id', 'n_occurrences'],
            on='tool',
            values='fraction'
        )\
        .fill_null(0)
    
    return spacer_counts_with_tools

def plot_tool_performance(spacer_counts_with_tools, tools_list, n_high_occ_bins=3, 
                         output_prefix='results/real_data/plots/tool_performance', 
                         max_bin=3, n_bins=150):
    """Plot tool performance vs number of occurrences.
    
    Args:
        spacer_counts_with_tools (pl.DataFrame): DataFrame with tool performance per spacer
        tools_list (list): List of tools to analyze
        n_high_occ_bins (int): Number of bins for high-occupancy region (>10^3)
        output_prefix (str): Prefix for output files
        max_bin (int): Maximum power of 10 for regular bins
        n_bins (int): Number of regular bins
    """
    # Define unique markers for each tool
    marker_dict = {
        'bowtie1': 'o',      # circle
        'bowtie2': 's',      # square
        'minimap2': '^',     # triangle up
        'bbmap_skimmer': 'D',  # diamond
        'spacer_containment': 'v',  # triangle down
        'mummer4': '<',      # triangle left
        'lexicmap': '>',     # triangle right
        'mmseqs2': 'p',      # pentagon
        'blastn': 'h',       # hexagon
        'strobealign': '8',  # octagon
        'minimap2_og': '*',  # star
        'minimap2_mod': 'P', # plus
        'bbmap_skimmer_mod': 'X',  # x
        'mmseqs2map': 'd'    # thin diamond
    }
    
    # Create range bins for number of occurrences
    bins = np.logspace(np.log10(1), max_bin, n_bins)

    # Calculate mean fraction for each tool within each bin
    bin_stats = []
    for i in range(len(bins)-1):
        mask = (spacer_counts_with_tools['n_occurrences'] >= bins[i]) & \
               (spacer_counts_with_tools['n_occurrences'] < bins[i+1])
        bin_data = spacer_counts_with_tools.filter(mask)
        if bin_data.height > 0:
            stats = {
                'bin_start': bins[i],
                'bin_end': bins[i+1],
                'n_spacers': bin_data.height
            }
            for tool in tools_list:
                stats[tool] = bin_data[tool].mean()
            bin_stats.append(stats)

    # Add points for high occurrences in multiple bins
    if n_high_occ_bins > 0:
        high_occ_edges = np.logspace(3, 4, n_high_occ_bins + 1)
        for i in range(n_high_occ_bins):
            bin_start = high_occ_edges[i]
            bin_end = high_occ_edges[i + 1]
            
            if i == n_high_occ_bins - 1:
                high_occ_mask = (spacer_counts_with_tools['n_occurrences'] >= bin_start)
            else:
                high_occ_mask = (spacer_counts_with_tools['n_occurrences'] >= bin_start) & \
                               (spacer_counts_with_tools['n_occurrences'] < bin_end)
            
            high_occ_data = spacer_counts_with_tools.filter(high_occ_mask)
            if high_occ_data.height > 0:
                high_occ_stats = {
                    'bin_start': bin_start,
                    'bin_end': bin_end,
                    'n_spacers': high_occ_data.height
                }
                for tool in tools_list:
                    high_occ_stats[tool] = high_occ_data[tool].mean()
                bin_stats.append(high_occ_stats)

    # Create the plot
    plt.figure(figsize=(15, 8))

    for tool in tools_list:
        x = [(stat['bin_start'] + stat['bin_end'])/2 for stat in bin_stats]
        y = [stat[tool] for stat in bin_stats]
        plt.plot(x, y, label=tool, marker=marker_dict.get(tool, 'o'), markersize=8, linewidth=1)

    plt.xscale('log')
    plt.xlabel('Number of occurrences (log scale)')
    plt.ylabel('Mean Detection Fraction')
    plt.title(f'Tool Performance vs Number of occurrences ({n_high_occ_bins} high-occ bins)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.grid(True, which="major", ls="-", alpha=0.5)
    plt.minorticks_on()
    plt.ylim(0, 1.05)
    plt.xlim(1, 10**4)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_detailed_{n_high_occ_bins}bins.pdf', bbox_inches='tight')
    plt.close()

    # Print statistics for every n-th bin
    n_th = 5
    print(f"\nSelected Bin Statistics (every {n_th}th bin):")
    print(f"{'Bin Range':<30} {'Number of Spacers':<20}")
    print("-" * 50)
    for i, stat in enumerate(bin_stats):
        if i % n_th == 0:
            bin_range = f"{stat['bin_start']:.1f}-{stat['bin_end']:.1f}"
            print(f"{bin_range:<30} {stat['n_spacers']:<20}")
    
    # Save full statistics
    pl.DataFrame(bin_stats).write_csv(f'{output_prefix}_detailed_stats.tsv', separator='\t')
    
    return bin_stats 