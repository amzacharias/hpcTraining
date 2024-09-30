#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Setup metadata for quality control
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
library(dplyr) # 1.1.4

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
# Input ===========
rDataDir <- file.path(projDir, "1_qc", "rDataDir")

# Output ===========

# Load data -----------------------------------------------
mergedSeurat <- readRDS(file = file.path(rDataDir, "mergedSeurat.rds"))

# Calculate quality metrics -----------------------------------------------
# Number of genes detected per UMI e.g., "Novelty Score" ===========
# Measure of complexity; more genes per UMI --> more complexity
mergedSeurat$log10GenesPerUMI <- log10(mergedSeurat$nFeature_RNA) / log10(mergedSeurat$nCount_RNA)

# Mitochondrial ratio ===========
# Percentage of cell reads originating from mitochondrial genes
#   Get the percentage of genes that match the regex pattern for genes starting with "MT-"
# The Seurat::PercentageFeatureSet() function doesn't work if you're using gene names etc.
#   Link to calculate from scratch: https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md
mergedSeurat$mitoRatio <- Seurat::PercentageFeatureSet(object = mergedSeurat, pattern = "^MT-")
# Convert to ratio by dividing by 100
mergedSeurat$mitoRatio <- mergedSeurat@meta.data$mitoRatio / 100

# Additional metadata as separate object -----------------------------------------------
metadata <- mergedSeurat@meta.data

# Add cell IDs to metadata ===========
metadata$cells <- rownames(metadata)

# Add sample column indicating ctrl vs stm ===========
# Remove all character after the first underscore, including the underscore
metadata$sample <- sub("_.+", "", metadata$cells)

# Rename columns ===========
metadata <- metadata %>%
  dplyr::rename(
    seqFolder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA
  )

# Add back to Seurat metadata ===========
mergedSeurat@meta.data <- metadata

# Save -----------------------------------------------
saveRDS(mergedSeurat, file = file.path(rDataDir, "mergedSeuratQcMetrics.rds"))
saveRDS(metadata, file = file.path(rDataDir, "metadata.rds"))
