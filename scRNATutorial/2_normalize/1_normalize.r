#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Normalize counts and regress out unwanted variation
# Author: Amanda Zacharias
# Date: 2024-09-30
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06_SC_SCT_normalization.md
# Goal: Account for differences in seq depth and overdispersion. Identify most variant genes to indicate cell types.
# Challenges: Checking and removing unwanted variation to minimize risk of artifacts.
# Recommendations:
#   1) Know what cell types to expect, high complexity? low mt content? are cells differentiating?
#   2) Regress out number of UMIs, mt content, and cell cycle if needed

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(tidyverse) # 2.0.0
library(RCurl) # 1.98.1.16
library(cowplot) # 1.1.3

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========
cycleDataPath <- file.path(projDir, "0_data", "cycle.rda")
inRDataDir <- file.path(projDir, "1_qc", "rDataDir")

## Output ===========
rDataDir <- file.path(projDir, "2_normalize", "rDataDir")
plotsDir <- file.path(projDir, "2_normalize", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
load(cycleDataPath)
filteredSeurat <- readRDS(file = file.path(inRDataDir, "filteredSeurat.rds"))

# Normalize -----------------------------------------------
# Adjusts for primarily sequencing depth and gene length
# Normalization as a two step process:
#   1) Scaling: each UMI count * cell specific factor --> all cells have same UMI counts
#   2) Transformation: Precise method varies across simple and complex approaches.
#     - Simple: standard log normalization approach. Hafemesiter C and Satija R Genom Biology 2019
#       issue: only low/medium abundance genes are correctly normalized.
#         very low = high variance, very high = low variance.
#     - Solution: Pearson residuals for transformation: Seurat's SCTransform()
#         measurements multiplied by gene-specifc weight,
#         genes w/ high evidence of non-uniform expression across cells have more weight.
#        Therefore, giving more value to genes expressed in few cells,
#        aiding rare cell population identification.
seuratPhase <- Seurat::NormalizeData(filteredSeurat)
# Normalizing layer: counts
# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|


# We could use SCTransform() directly, but it's more sophisticated and we want something more simple for now.

# Effects of cell cycle -----------------------------------------------
# Is cell cycle a major variation source?
## Assign each cell based on G2/M and S phase markers. ===========
# Seurat::CellCycleScoring() calculates phase based on canonical markers that we provide from `cycle.rda`.
# Info on preparing a marker list:
#   https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/cell_cycle_scoring.md
seuratPhase <- CellCycleScoring(
  seuratPhase,
  g2m.features = g2m_genes,
  s.features = s_genes
)
# Warning: The following features are not present in the object: RAD51, CDC45, E2F8, BRIP1, DTL, RRM2, EXO1, not searching for symbol synonyms # nolint: line_length_linter.
# Warning: The following features are not present in the object: ANLN, NEK2, HJURP, DLGAP5, PIMREG, KIF2C, CDC25C, CKAP2L, CDCA2, not searching for symbol synonyms # nolint: line_length_linter.

## Scaling ===========
# Scale data so highly expressed genes don't dominate the PCA.
# Mean expression for each gene will be 0 and have variance of 1.
# 1. identify most variable genes
seuratPhase <- Seurat::FindVariableFeatures(
  seuratPhase,
  selection.method = "vst", # default
  nfeatures = 2000, # default
  verbose = FALSE
)
# 2. scale the counts
seuratPhase <- ScaleData(seuratPhase)
# Centering and scaling data matrix
#   |======================================================================| 100%

## Get top 15 variable genes ===========
rankedVariableGenes <- Seurat::VariableFeatures(seuratPhase)
topGenes <- rankedVariableGenes[1:15]

## Plotting variable genes ===========
# Avg expression x variance of most variable genes w/ labels of top 15.
variableFtsPlot <- Seurat::VariableFeaturePlot(seuratPhase)
variableFtsPlot <- LabelPoints(plot = variableFtsPlot, points = topGenes, repel = TRUE)
ggsave(
  plot = variableFtsPlot, path = plotsDir, filename = "VariableFtsPlot.pdf",
  width = 180, height = 90, units = "mm"
)

## PCA ===========
seuratPhase <- RunPCA(seuratPhase)
pcaPhasePlot <- DimPlot(seuratPhase, reduction = "pca", group.by = "Phase", split.by = "Phase")
ggsave(
  plot = pcaPhasePlot, path = plotsDir, filename = "pcaPhasePlot.pdf",
  width = 180, height = 90, units = "mm"
)
# No large differences apparent due to cell cycle, so don't need to regress out!

# Effects of mitochondrial expression -----------------------------------------------
# Mt expression can bias clustering, but if it
#   represents a biological function, it shouldn't be regressed out.
## Make mitoRatio categorical with quantiles ===========
# Check quantile values
summary(seuratPhase@meta.data$mitoRatio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.01529 0.02109 0.02253 0.02801 0.14464

# To categorical
seuratPhase@meta.data$mitoFr <- cut(
  seuratPhase@meta.data$mitoRatio,
  breaks = c(0.00000, 0.01529, 0.02109, 0.02801, 0.14464),
  labels = c("Low", "Medium", "Medium high", "High")
)

## Plotting as PCA ===========
pcaMitoPlot <- DimPlot(seuratPhase, reduction = "pca", group.by = "mitoFr", split.by = "mitoFr")
ggsave(
  plot = pcaMitoPlot, path = plotsDir, filename = "pcaMitoPlot.pdf",
  width = 180, height = 90, units = "mm"
)
# No obvious effect of mitochondrial expression.
# Each group is dispersed evenly across all the cells, meaning that no cells are biased
#   towards one mitochondrial expression level. In other words, cells cannot be
#   seprated by mitochondrial expression alone.
# No, I would not regress out mitochondrial fraction as unwanted variation.

# Regress out unwanted variation -----------------------------------------------
# Normalize & regress out unwanted variation during variance transformation.
# SCTranform() constructs a GLM for each gene, response = UMI counts, explanatory var = sequencing depth.
## Separate samples from different conditions ===========
splitSeurat <- Seurat::SplitObject(seuratPhase, split.by = "sample")
# If we only want the control replicates, we could slice the split ctrl object
ctrlReps <- splitSeurat[c("ctrl_1", "ctrl_2")]

## Allow large objects option ===========
# I would usually put this at the top of the script, but leaving here for the tutorial.
options(future.globals.maxSize = 4000 * 1024^2)

## Run SCTransform on objects ===========
for (set in seq_along(splitSeurat)) {
  splitSeurat[[set]] <- SCTransform(
    splitSeurat[[set]],
    vars.to.regress = c("mitoRatio"),
    vst.flavor = "v2"
  )
}
# v2 of VST was introduced in early 2022;
# has improved speed/memory, stability of parameter estimates, and variable feature identification.
# Example output message: 
# Running SCTransform on assay: RNA
# Running SCTransform on layer: counts
# vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
# Variance stabilizing transformation of count matrix of size 14244 by 14841
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 2000 genes, 5000 cells
# Found 313 outliers - those will be ignored in fitting/regularization step

# Second step: Get residuals using fitted parameters for 14244 genes
# Computing corrected count matrix for 14244 genes
# Calculating gene attributes
# Wall clock passed: Time difference of 38.75653 secs
# Determine variable features
# Regressing out mitoRatio
#   |======================================================================| 100%
# Centering data matrix
#   |======================================================================| 100%
# Getting residuals for block 1(of 3) for counts dataset
# Getting residuals for block 2(of 3) for counts dataset
# Getting residuals for block 3(of 3) for counts dataset
# Regressing out mitoRatio
#   |======================================================================| 100%
# Centering data matrix
#   |======================================================================| 100%
# Finished calculating residuals for counts
# Set default assay to SCT

# If there are many cell, it may be useful to output more than 3,000 of the
# top most variable genes. Change with th `variable.features.n` parameter.

# The "Set default assay to SCT" message means that we will use the data after
# SCT on downstream analyses by default.
splitSeurat$ctrl@assays
# $RNA
# Assay (v5) data with 14244 features for 14841 cells
# Top 10 variable features:
#  HBB, CCL3, CXCL10, HBA2, CCL4L2, CCL4, TXN, PPBP, CCL2, CCL7
# Layers:
#  counts, data, scale.data

# $SCT
# SCTAssay data with 14244 features for 14841 cells, and 1 SCTModel(s)
# Top 10 variable features:
#  FTL, CCL2, IGKC, GNLY, IGLC2, CCL3, TIMP1, IGHM, PPBP, CCL4

# Q & A -----------------------------------------------
## Q1: Are the same assays available for the "stim" samples?
# Yes.
splitSeurat$stim@assays
# $RNA
# Assay (v5) data with 14244 features for 14841 cells
# Top 10 variable features:
#  HBB, CCL3, CXCL10, HBA2, CCL4L2, CCL4, TXN, PPBP, CCL2, CCL7
# Layers:
#  counts, data, scale.data

# $SCT
# SCTAssay data with 14244 features for 14841 cells, and 1 SCTModel(s)
# Top 10 variable features:
#  FTL, CCL2, IGKC, GNLY, IGLC2, CCL3, TIMP1, IGHM, PPBP, CCL4

## Q2: Are the top features different ctrl vs stim?
# They appear to be the same

# Save -----------------------------------------------
saveRDS(splitSeurat, file = file.path(rDataDir, "splitSeurat.rds"))
