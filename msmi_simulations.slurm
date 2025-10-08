#!/bin/bash

#SBATCH --job-name=msmi_simulations
#SBATCH --mail-user=rbtucker@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=5:00:00
#SBATCH --account=philb0
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.log
#SBATCH --array=1,2

cd /home/rbtucker/msmi_simulations
module load R

Rscript sim_illness_death.R "$SLURM_ARRAY_TASK_ID"
