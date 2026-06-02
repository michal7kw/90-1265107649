#!/bin/bash

#SBATCH --job-name=6b_visualize_splicing
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6b_visualize_splicing.err"
#SBATCH --output="./logs/6b_visualize_splicing.out"

# Load conda - using R environment for visualization
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"

cd ${BASE_DIR}

Rscript scripts/6b_visualize_splicing.R