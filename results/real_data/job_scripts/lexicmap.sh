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
#SBATCH -J lexicmap   
## specify output and error file
#SBATCH -o /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/slurm_logs/lexicmap-%A.out
#SBATCH -e /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/slurm_logs/lexicmap-%A.err
bash /clusterfs/jgi/scratch/science/metagen/neri/code/blits/spacer_matching_bench/results/real_data/bash_scripts/lexicmap.sh $SLURM_CPUS_PER_TASK

