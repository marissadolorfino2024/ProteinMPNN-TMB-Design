#!/usr/bin/env bash
  
#SBATCH --job-name=frelax_designs
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 72:00:00
#SBATCH --array=1-488

file=$(cat tasks.txt | sed -n "$SLURM_ARRAY_TASK_ID p")

/home/mdolo/turbo/mdolo/opt/ColabFold/localcolabfold/colabfold-conda/bin/colabfold_batch --msa-mode single_sequence --num-recycle 48 $file ${file::-6}_folds
