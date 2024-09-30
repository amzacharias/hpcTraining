#!/bin/bash
#SBATCH --job-name=downloadProj
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=20MB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=downloadProj.out
#SBTACH --error=downloadProj.err
# Title: Download project files from hpc github
# Author: Amanda Zacharias
# Date: 2024-09-19
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
#
# Code -------------------------------------------
echo Job started at $(date +%T)
# Dependencies
module load StdEnv/2023
# Variables
DATADIR=~/hpcTraining/scRNATutorial/0_data

# Download and unzip
wget -P $DATADIR https://www.dropbox.com/s/vop78wq76h02a2f/single_cell_rnaseq.zip?dl=1
unzip single_cell_rnaseq.zip?dl=1

# Pull out data files from unzipped folder
cp -r single_cell_rnaseq/data/* $DATADIR

# Delete unncessary files
rm -r __MACOSX/
rm single_cell_rnaseq.zip?dl=1
rm -r single_cell_rnaseq

echo Job ended at $(date +%T)