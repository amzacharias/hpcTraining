#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Load the scRNA samples and merge into one Seurat object.
# Author: Amanda Zacharias
# Date: 2024-09-19
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/03_SC_quality_control-setup.md

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0

# Source -----------------------------------------------

# Pathways -----------------------------------------------
# Input ===========
dataDir <- file.path("0_data")
sampleDataDirs <- list.files(dataDir, full.names = TRUE,pattern = "matrix")
names(sampleDataDirs) <- sub("(.+?)(\\_.*)", "\\1", basename(sampleDataDirs))

# Output ===========
rDataDir <- file.path("1_qc", "rDataDir")
dir.create(rDataDir, showWarnings = FALSE)

# Load data -----------------------------------------------
# Two options Matrix::readMM() or Seurat::Read10X(); we're choosing the latter
# Function to read in a sample =========
Read10XSample <- function(inPath) {
  #' Function to read in a 10X genomics sample and output as Seurat object
  #'
  #' @param inPath Path to input files, which are output from cellranger
  #' @return Returns a seurat object
  # Read in data as a sparse matrix
  cnts <- Read10X(data.dir = inPath)

  # Sparse matrix --> Seurat object
  #   Cells with less than 100 genes detected are not included in analysis
  #   because they are likely random barcodes w/o any cell.
  seuratObj <- CreateSeuratObject(counts = cnts, min.features = 100)

  # Return
  return(seuratObj)
}
# Execute =======
# If there were more samples, you could use a loop
ctrl <- Read10XSample(sampleDataDirs[["ctrl"]])
stim <- Read10XSample(sampleDataDirs[["stim"]])

# Explore metadata -----------------------------------------------
head(ctrl@meta.data)
#                     orig.ident nCount_RNA nFeature_RNA
# AAACATACAATGCC-1 SeuratProject       2344          874
# AAACATACATTTCC-1 SeuratProject       3125          896
# AAACATACCAGAAA-1 SeuratProject       2578          725
# AAACATACCAGCTA-1 SeuratProject       3261          979
# AAACATACCATGCA-1 SeuratProject        746          362
# AAACATACCTCGCT-1 SeuratProject       3519          866

# orig.ident: sample identity if known. default = SeuratProject
# nCount_RNA: number of UMIs per cell
# nFeature_RNA: number of genes per cell

head(stim@meta.data)
#                     orig.ident nCount_RNA nFeature_RNA
# AAACATACCAAGCT-1 SeuratProject       1221          606
# AAACATACCCCTAC-1 SeuratProject       1782          807
# AAACATACCCGTAA-1 SeuratProject       1451          605
# AAACATACCCTCGT-1 SeuratProject       1549          747
# AAACATACGAGGTG-1 SeuratProject       1303          558
# AAACATACGCGAAG-1 SeuratProject       5445         1330

# Merge samples -----------------------------------------------
# Create merged object ==========
#   Adding a prefix to the cell ids because different samples can have the same cell ids.
#   If there were many samples, you give y as a vector e.g., y = c(sample2, sample3, sample4) etc.
mergedSeurat <- merge(x = ctrl, y = ctrl, add.cell.id = c("ctrl", "stim"))

# Concatenate count matrixes together ===========
mergedSeurat <- JoinLayers(mergedSeurat)

# Inspect ===========
head(mergedSeurat@meta.data)
#                          orig.ident nCount_RNA nFeature_RNA
# ctrl_AAACATACAATGCC-1 SeuratProject       2344          874
# ctrl_AAACATACATTTCC-1 SeuratProject       3125          896
# ctrl_AAACATACCAGAAA-1 SeuratProject       2578          725
# ctrl_AAACATACCAGCTA-1 SeuratProject       3261          979
# ctrl_AAACATACCATGCA-1 SeuratProject        746          362
# ctrl_AAACATACCTCGCT-1 SeuratProject       3519          866
tail(mergedSeurat@meta.data)
#                          orig.ident nCount_RNA nFeature_RNA
# stim_TTTGCATGCGAATC-1 SeuratProject        376          172
# stim_TTTGCATGCTTCGC-1 SeuratProject       2473          892
# stim_TTTGCATGGCAGTT-1 SeuratProject        752          343
# stim_TTTGCATGGGAACG-1 SeuratProject        889          403
# stim_TTTGCATGGTCCTC-1 SeuratProject       1326          543
# stim_TTTGCATGTTCATC-1 SeuratProject       1507          510

# Save -----------------------------------------------
saveRDS(mergedSeurat, file = file.path(rDataDir, "mergedSeurat.rds"))
