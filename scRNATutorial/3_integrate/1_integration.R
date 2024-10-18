#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Integrate cells across conditions
# Author: Amanda Zacharias
# Date: 2024-10-18
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06_integration.md
# If cells cluster by sample, condition, batch, dataset, modality, etc. integration can
# improve clustering and downstream analyses.
# Integration uses shared highly variable genes to harmonize groups.
# Goal: ensure cell types of one condition/dataset align with same cells of other conditions/datasets.
# Seurat uses canonical correlation analysis (CCA) method.
#   Expects "correspondences"/shared states among at least a subset of cells across groups.
# Review on integration methods: Leucken et al. 2022, 10.1038/s41592-021-01336-8

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========
inRDataDir <- file.path(projDir, "2_normalize", "rDataDir")

## Output ===========
rDataDir <- file.path(projDir, "3_integrate", "rDataDir")
plotsDir <- file.path(projDir, "3_integrate", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
splitSeurat <- readRDS(file.path(inRDataDir, "splitSeurat.rds"))

# Specify integration will be w/ 3000 most variable genes -----------------------------------------------
# Default is 2,000 genes
integFeatures <- SelectIntegrationFeatures(object.list = splitSeurat, nfeatures = 3000)

# Prepare SCTransform obj. for integration -----------------------------------------------
splitSeurat <- PrepSCTIntegration(object.list = splitSeurat, anchor.features = integFeatures)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04s

# Perform CCA -----------------------------------------------
## Find the best buddies and filter incorrect anchors ===========
# Multiple nearest neighbors (MNN) is used to identify anchors.
# Progress bar will stay at 0% but it is running
integAnchors <- FindIntegrationAnchors(
  object.list = splitSeurat,
  normalization.method = "SCT",
  anchor.features = integFeatures
)
# Finding all pairwise anchors
#   |                                                  | 0 % ~calculating  Running CCA
# Merging objects
# Finding neighborhoods
# Finding anchors
#         Found 41211 anchors
# Filtering anchors
#         Retained 32881 anchors
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=08m 39s

## Integrate ===========
seuratIntegrated <- IntegrateData(
  anchorset = integAnchors,
  normalization.method = "SCT"
)
# [1] 1
# Warning: Layer counts isn't present in the assay object; returning NULL
# [1] 2
# Warning: Layer counts isn't present in the assay object; returning NULL
# Merging dataset 2 into 1
# Extracting anchors for merged samples
# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Integrating data
# Warning: Layer counts isn't present in the assay object; returning NULL
# Warning: Assay integrated changing from Assay to SCTAssay
# Warning: Layer counts isn't present in the assay object; returning NULL
# Warning messages:
# 1: Different cells and/or features from existing assay SCT
# 2: Different cells and/or features from existing assay SCT
# 3: Different cells and/or features from existing assay SCT

# Visualization -----------------------------------------------
# With PCA, can only show 2 dimensions.
# UMAP is better at representing global structure.
## PCA ===========
seuratIntegrated <- RunPCA(object = seuratIntegrated)
pcaPlot <- PCAPlot(seuratIntegrated, split.by = "sample")
ggsave(
  plot = pcaPlot,
  filename = "pcaPlot_aftCCC.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

## UMAP ===========
# UMAP has randomness, so a seed is necessary.
set.seed(123456)
seuratIntegrated <- RunUMAP(
  seuratIntegrated,
  dims = 1:40,
  reduced = "pca"
)
umapPlot <- DimPlot(seuratIntegrated, combine = TRUE)
ggsave(
  plot = umapPlot,
  filename = "umapPlot_aftCCC.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

## UMAP Split ===========
# Split by condition
umapSplitPlot <- DimPlot(
  seuratIntegrated,
  split.by = "sample",
  combine = TRUE
)
ggsave(
  plot = umapSplitPlot,
  filename = "umapSplitPlot_aftCCC.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

# Save -----------------------------------------------
saveRDS(seuratIntegrated, file = file.path(rDataDir, "seuratIntegrated.rds"))

# Complex integration -----------------------------------------------
# In the case of a very complex scenario, a different integration approach may be better.
# See review cited at top of script.
# `Harmony` is another tool for complex integration which has had some success.
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06c_integration_harmony.md