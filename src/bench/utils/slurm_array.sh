#!/bin/bash
## specify an email address
#SBATCH --mail-user=uneri@lbl.gov
## specify when to send the email when job is (a)borted, (b)egins, or (e)nds
#SBATCH --mail-type=FAIL
## specify allocation - we want jgi_shared since we don't want to use the whole node for nothing
#SBATCH -A grp-org-sc-metagen
#SBATCH -q jgi_normal
## specify number of nodes
#SBATCH -N 1
#######SBATCH --exclusive
## specify number of procs
#SBATCH -c 25
## specify ram
#SBATCH --mem=50G 
## specify runtime
#SBATCH -t 24:00:00
## specify job name
#SBATCH -J spacer_matching_bench
## specify output and error file
#SBATCH -o /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/slurmsout/Slurmout-%A_%a.out
#SBATCH -e /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/slurmsout/Slurmout-%A_%a.err

## specify that we never run more than XX jobs at a time (using "%", e.g. --array=0-5%2)
#SBATCH --array=0-5%1


cd /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/
eval "$(mamba shell hook --shell bash)"
mamba activate base_env
# echo $SLURM_MEM_PER_NODE

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"

python bench.py --contig_length_range 1501 200000 \
--spacer_length_range 18 120 \
--n_mismatch_range $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID \
--sample_size_contigs 40000 \
--sample_size_spacers 1225 \
--insertion_range 1 1225 \
--threads $SLURM_CPUS_PER_TASK \
--prop_rc 0.5 \
--max_runs 1 \
--max_mismatches 7 