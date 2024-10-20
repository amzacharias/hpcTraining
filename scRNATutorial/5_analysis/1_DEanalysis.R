#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Begin investigating whether cell types differ ctrl vs stim
# Author: Amanda Zacharias
# Date: 2024-10-20
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========
inRDataDir <- file.path(projDir, "4_cluster", "rDataDir")

## Output ===========
rDataDir <- file.path(projDir, "5_analysis", "rDataDir")
plotsDir <- file.path(projDir, "5_analysis", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
seuratIntegrated <- readRDS(file.path(inRDataDir, "finalSeuratIntegrated.rds"))

# Add cell type annotation to metadata -----------------------------------------------
seuratIntegrated$celltype <- Idents(seuratIntegrated)

# Examine number of cells per celltype -----------------------------------------------
## Calculate n cells per celltype ===========
nCells <- FetchData(seuratIntegrated, vars = c("celltype", "sample")) %>%
  dplyr::count(celltype, sample)
## Barplot ===========
nCellsPerCelltypeBarplot <- nCells %>%
  ggplot(aes(x = celltype, y = n, fill = sample)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_text(aes(label = n), vjust = -0.2, position = position_dodge(1)) +
  theme_bw()
ggsave(
  plot = nCellsPerCelltypeBarplot,
  filename = "nCellsPerCelltypeBarplot.pdf", path = plotsDir,
  height = 90, width = 270, units = "mm"
)

# Pseudobulk DE analysis -----------------------------------------------
## Just B-cells ===========
seuratBCells <- subset(seuratIntegrated, subset = (celltype == "B cells"))
## Wilcox test ===========
bMarkers <- FindMarkers(
  seuratBCells,
  ident.1 = "ctrl",
  ident.2 = "stim",
  grouping.var = "sample",
  only.pos = FALSE,
  logfc.threshold = 0.25
)
# Error in WhichCells.Seurat(object = object, idents = ident.1) : 
#   Cannot find the following identities in the object: ctrl

# If it had worked, you could then make a volcano plot from the results.
# The tutorial recommends the `EnhancedVolcano` package.

# Description of pseudobulk analysis in more detail: https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md

# Other options -----------------------------------------------
# Sub-clustering of cells ex. t-cells: https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/09_merged_SC_marker_identification.md#:~:text=cells%20as%20described-,here,-Trajectory%20analysis%2C%20or
#   you could simply run FindNeighbors and FindClusters() on your subset() of cells.
#   you might also start with the raw counts and re-run SCTransform to get most variable genes again.
#   recall that your # of cells will decrease, which should be considered ex. for setting number of K neighbors
# Trajectory analysis or lineage tracing...
