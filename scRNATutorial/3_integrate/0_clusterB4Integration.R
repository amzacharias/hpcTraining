#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Clustering before integration
# Author: Amanda Zacharias
# Date: 2024-10-18
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06_integration.md
# Look at clustering of cells before sample integration b.c. it may be
# interesting to see differences between conditions beforehand.
# If we normalized ctrl and stim samples together, we would get condition-specific clustering.
# See below:
# This isn't ideal because we want cells of the same type to group together,
# facilitating downstream analyses.

# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========
cycleDataPath <- file.path(projDir, "0_data", "cycle.rda")
inRDataDir <- file.path(projDir, "2_normalize", "rDataDir")

## Output ===========
rDataDir <- file.path(projDir, "3_integrate", "rDataDir")
plotsDir <- file.path(projDir, "3_integrate", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
seuratPhase <- readRDS(file.path(inRDataDir, "seuratPhase.rds"))

# Generate UMAP -----------------------------------------------
seuratPhase <- RunUMAP(seuratPhase, dims = 1:40, reduction = "pca")
umapPlot <- DimPlot(seuratPhase, combine = TRUE)
# adding the combine=TRUE, so output is a patchwork ggplot object

ggsave(
  plot = umapPlot,
  filename = "b4integration_umapPlot.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

# The plot didn't really come out as expcted. We should see that
# cells from stim are clearly not clustering with cells from ctrl.