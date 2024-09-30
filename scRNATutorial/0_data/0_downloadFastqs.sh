#!/bin/bash
#SBATCH --job-name=downloadFastqs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=downloadFastqs.out
#SBTACH --error=downloadFastqs.err
# Title: Download scRNA fastq files from SRA
# Author: Amanda Zacharias
# Date: 2024-09-19
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# For the tutorial, starting with the data processed by hbc but here's an example of downloading
# the fastq files.
# It seems odd that there's only one fastq file per sample, though only one read is described in SRA.
# This approach may be suboptimal for that reason.
#
# Data from: GEO:GSE96583
# hbctraining tutorial suggests they downloaded bam files then converted to fastq, instead of directly
# downloading fastq files.
# Code -------------------------------------------
echo Job started at $(date +%T)

# Dependencies
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9

# Variables
DATADIR=~/hpcTraining/scRNATutorial/0_data
FASTQDIR=${DATADIR}/fastq
mkdir $FASTQDIR

# Read in accession IDs to array
mapfile -t ACCESSIONIDS < ${DATADIR}/SRR_Acc_List.txt

# Download Function
download_fastq () {
  echo $1
  # Check if file has already been downloaded
  if [ -f ${1}_1.fastq.gz ]; then
    echo File ${1}_1.fastq.gz exists.
  else
    echo File ${1}_1.fastq.gz does not exist
    # Downloading begins
    prefetch --output-directory ${FASTQDIR} $1
    fastq-dump --split-files --gzip --outdir $FASTQDIR $1
  fi
}

# Execute with parallelization
N=2
(
for id in ${ACCESSIONIDS[@]}; do
  ((i=i%N)); ((i++==0)) && wait
  download_fastq $id &
done
)

echo Job ended at $(date +%T)