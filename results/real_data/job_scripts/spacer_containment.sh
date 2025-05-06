#!/bin/bash
#SBATCH --mail-user=uneri@lbl.gov
#SBATCH --mail-type=FAIL,END,BEGIN
#SBATCH -A grp-org-sc-metagen
#SBATCH -q jgi_normal
#SBATCH -c 30
#SBATCH --mem=168G 
## specify runtime
#SBATCH -t 72:00:00
## specify job name
#SBATCH -J spacer_containment   
## specify output and error file
#SBATCH -o /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/slurm_logs/spacer_containment-%A.out
#SBATCH -e /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/slurm_logs/spacer_containment-%A.err
bash /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/bash_scripts/spacer_containment.sh $SLURM_CPUS_PER_TASK

