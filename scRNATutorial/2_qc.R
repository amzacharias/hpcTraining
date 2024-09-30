#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Quality control of merged seurat object
# Author: Amanda Zacharias
# Date: 2024-09-23
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md
# Goal: only include high quality true cells and identify failed samples
# Challenges: discriminating poor quality cells from less complex cells, choosing appropriate thresholds.
# Recommendations: Use prior knowledge to inform what cell types  you expect to be present.
#
# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0

# Source -----------------------------------------------

# Pathways -----------------------------------------------
# Input ===========
rDataDir <- file.path("1_qc", "rDataDir")

# Output ===========

# Load data -----------------------------------------------
mergedSeurat <- readRDS(file = file.path(rDataDir, "mergedSeurat.rds"))

# Calculate quality metrics -----------------------------------------------
# Number of genes detected per UMI e.g., "Novelty Score" ===========
# Measure of complexity; more genes per UMI --> more complexity
mergedSeurat$log10GenesPerUMI <- log10(mergedSeurat$nFeature_RNA) / log10(mergedSeurat$nCount_RNA)

# Mitochondrial ratio ===========
# Percentage of cell reads originating from mitochondrial genes
mergedSeurat$mitoRatio <- Seurat::PercentageFeatureSet(object = mergedSeurat, pattern = "^MT-")
# Get the percentage of genes that match the regex pattern for genes starting with "MT-"