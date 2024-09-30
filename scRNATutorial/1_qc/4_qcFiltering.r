#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Quality control filtering
# Author: Amanda Zacharias
# Date: 2024-09-30
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# Example of QC on low quality data https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/QC_bad_data.md

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
mergedSeurat <- readRDS(file = file.path(rDataDir, "mergedSeuratQcMetrics.rds"))

# Cell-level filtering -----------------------------------------------
# Here are rough guidelines; make sure your thresholds make sense for your experiment.
filteredSeurat <- mergedSeurat %>%
  subset(
    (nUMI > 500) &
      (nGene >= 250) &
      (log10GenesPerUMI > 0.80) &
      (mitoRatio < 0.20)
  )
cat("\nNumber of cells before filtering:", ncol(mergedSeurat),
    "\nNumber of cells after filtering:", ncol(filteredSeurat),
    "\nNumber of cells removed:", ncol(mergedSeurat) - ncol(filteredSeurat))
# Number of cells before filtering: 31376
# Number of cells after filtering: 29682
# Number of cells removed: 1694

# Gene-level filtering -----------------------------------------------
# Remove genes with a zero count in >10 cells ===========
# Extract counts
counts <- GetAssayData(object = filteredSeurat, layer = "counts")
# Logical matrix; TRUE = count is > 10
nonZero <- counts > 0
# Sum all TRUEs. For each row/gene, return TRUE if number of TRUEs is >= 10
keepGenes <- Matrix::rowSums(nonZero) >= 10
# Subset
filteredCounts <- counts[keepGenes, ]

cat("\nNumber of genes before filtering:", nrow(counts),
    "\nNumber of genes after filtering:", nrow(filteredCounts),
    "\nNumber of genes removed:", nrow(counts) - nrow(filteredCounts))
# Number of genes before filtering: 33538
# Number of genes after filtering: 14244
# Number of genes removed: 19294

# New Seurat object for downstream analyses -----------------------------------------------
filteredSeurat <- CreateSeuratObject(filteredCounts, meta.data = filteredSeurat@meta.data)

# Save -----------------------------------------------
saveRDS(filteredSeurat, file = file.path(rDataDir, "filteredSeurat.rds"))
saveRDS(filteredSeurat@meta.data, file = file.path(rDataDir, "filteredMetadata.rds"))

# Re-generate QC plots with the previously run script! ###############################################
