#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Clustering Quality Control; Are the clusters true?
# Author: Amanda Zacharias
# Date: 2024-10-18
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/08_SC_clustering_quality_control.md
#
# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(purrr) # 1.0.2; used for `map()`
library(patchwork) # 1.3.0

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========

## Output ===========
rDataDir <- file.path(projDir, "4_cluster", "rDataDir")
plotsDir <- file.path(projDir, "4_cluster", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
seuratIntegrated <- readRDS(file.path(rDataDir, "seuratIntegratedHbc.rds"))

# Explore QC Metrics -----------------------------------------------
## Segregate clusters by sample ===========
# Get identity and sample information then count number of cells per cluster per sample
nCells <- Seurat::FetchData(seuratIntegrated, vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)
head(nCells)
#   ident sample    n
# 1     0   ctrl 2017
# 2     0   stim 2203
# 3     1   ctrl 2398
# 4     1   stim 1320
# 5     2   ctrl 1840
# 6     2   stim 1809

# Barplot of number of cells per cluster by sample
clustersBySampleBarplot <- nCells %>%
  ggplot(aes(x = ident, y = n, fill = sample)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_text(aes(label = n), vjust = -0.2, position = position_dodge(1), size = 2) +
  scale_x_discrete(name = "Cluster identity") +
  scale_y_continuous(name = "Number of cells") +
  ggtitle("Number of cells per cluster") +
  theme_bw()
ggsave(
  plot = clustersBySampleBarplot,
  filename = "clustersBySampleBarplot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)

# UMAP cells per cluster per sample
clustersBySampleUmapPlot <- DimPlot(
  seuratIntegrated, label = TRUE, split.by = "sample"
) + NoLegend()
ggsave(
  plot = clustersBySampleUmapPlot,
  filename = "clustersBySampleUmapPlot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)

# Proportion of cells from a sample in each cluster
propCellsPerClustPlot <- seuratIntegrated@meta.data %>%
  ggplot() +
  geom_bar(aes(x = integrated_snn_res.0.8, fill = sample), position = position_fill()) +
  scale_x_discrete(name = "Cluster identity") +
  scale_y_continuous(name = "Cell / cluster ratio") +
  theme_bw()
ggsave(
  plot = propCellsPerClustPlot,
  filename = "propCellsPerClustPlot.pdf", path = plotsDir,
  width = 90, height = 90, units = "mm"
)

# We expect to see an even distribution of conditions across clusters,
# though depending on the experiment, some clusters might be condition specific.
# No clear clustering by sample.

# Segregate clusters by cell cycle phase -----------------------------------------------
# We didn't regress out cell cycle phase during SCTransform. Let's see if that
#   needs to be revisited.
clustersByPhaseUmapPlot <- DimPlot(
  seuratIntegrated, label = TRUE, split.by = "Phase"
) + NoLegend()
ggsave(
  plot = clustersByPhaseUmapPlot,
  filename = "clustersByPhaseUmapPlot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)
# No clear clustering by phase.

# Segregate clusters by var. uninteresting variation -----------------------------------------------
# Metrics of interest in meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

# Plotting
varMetricsPlot <- FeaturePlot(
  seuratIntegrated,
  reduction = "umap",
  features = metrics,
  pt.size = 0.4,
  order = TRUE, # plot cells in order of expression, high expression = above
  min.cutoff = "q10", # lower 10 quanitle of lowly expressed genes will be grayed out
  label = TRUE
)
ggsave(
  plot = varMetricsPlot,
  filename = "varMetricsPlot.pdf", path = plotsDir,
  width = 180, height = 270, units = "mm"
)

# Metrics rel. even across clusters except nGene is high in top-left clusters.

## Boxplot nGene per cluster ===========
nGenePerClusterBarplot <- seuratIntegrated@meta.data %>%
  ggplot(aes(x = integrated_snn_res.0.8, y = nGene, fill = integrated_snn_res.0.8)) +
  geom_boxplot() +
  scale_x_discrete("Cluster identity") +
  scale_y_continuous("Number of genes") +
  theme_bw() +
  Seurat::NoLegend()
ggsave(
  plot = nGenePerClusterBarplot,
  filename = "nGenePerClusterBarplot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)
# Boxplot shows a similar trend

# Explore PCs driving clustering -----------------------------------------------
# Hope defined PCs will be good at differentiating clusters.
## For each cell, extract UMAP coordinates and PC scores ===========
columnsOfInt <- c(paste0("PC_", 1:16), "ident", "UMAP_1", "UMAP_2")
pcData <- FetchData(seuratIntegrated, vars = columnsOfInt)

# A synonymous method to extract the reduction metadata:
seuratIntegrated@reductions$pca@cell.embeddings[1:10, 1:16]
seuratIntegrated@integrated_snn_res.0.8
seuratIntegrated@reductions$umap@cell.embeddings[1:10, 1:2]

# NOTE: Newer Seurat objects will have "umap_1" NOT "UMAP_1"

## UMAP plot for each PC ===========
umapLabels <- pcData %>%
  group_by(ident) %>%
  summarise(x = mean(UMAP_1), y = mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
umapPcPlots <- map(
  paste0("PC_", 1:16), function(pc) {
    pcData %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(colour = .data[[pc]]), alpha = 0.7) +
      scale_color_gradient(guide = "none", low = "grey90", high = "blue") +
      geom_text(data = umapLabels, aes(label = ident, x, y)) +
      ggtitle(pc)
  }
) %>%
  patchwork::wrap_plots()
ggsave(
  plot = umapPcPlots,
  filename = "umapPcPlots.pdf", path = plotsDir,
  width = 360, height = 360, units = "mm"
)

# We can see what we hoped to see!
# Example, PC_2 seems to drive clusters 8 and 12.
# We can relate this to the top variable genes positive for clusters 8 and 12.

# Explore known gene markers -----------------------------------------------
# Help determine if identity and resolution are appropriate.
# CD14+ monocytes	--> CD14, LYZ
# FCGR3A+ monocytes	--> FCGR3A, MS4A7
# Conventional dendritic cells --> FCER1A, CST3
# Plasmacytoid dendritic cells --> IL3RA, GZMB, SERPINF1, ITM2C
# B cells	--> CD79A, MS4A1
# T cells	--> CD3D
# CD4+ T cells --> CD3D, IL7R, CCR7
# CD8+ T cells --> CD3D, CD8A
# NK cells --> GNLY, NKG7
# Megakaryocytes --> PPBP
# Erythrocytes --> HBB, HBA2

## RNA counts as default assay ===========
# Note: SCTransform normalization was only top 3000 variable genes, so our genes of
# interest may not be present. So, switching to RNA counts as assay.
DefaultAssay(seuratIntegrated) <- "RNA"

# Normalize for visualization purposes
seuratIntegrated <- NormalizeData(seuratIntegrated, verbose = FALSE)
# In the assay, count slot = raw, data count = normalized

## CD14+ monocyte markers ===========
# Clusters must express both CD14 and LYZ to be considered true CD14+ monocytes
cd14MonocytePlot <- FeaturePlot(
  seuratIntegrated,
  reduction = "umap",
  features = c("CD14", "LYZ"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
ggsave(
  plot = cd14MonocytePlot,
  filename = "cd14MonocytePlot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)
# Clusters 1 & 3 have both markers. 10 is only high in 14, and CD14 only high in 10.
# ... etc. all the cell marker combinations

## Dotplot ===========
# Average expression of markers across clusters
# Can't use a gene twice.
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
dotPlot <- DotPlot(seuratIntegrated, markers, assay = "RNA")
ggsave(
  plot = dotPlot,
  filename = "dotPlot.pdf", path = plotsDir,
  width = 360, height = 90, units = "mm"
)

# Q & A -----------------------------------------------
# What clusters do you think correspond to each cell type?
# ...

# If a clusters appears to have multiple cell-types, consider increasing resolution.
# If that doesn't work, considering including more PCs.

# Remaining Q's:
# 1) Many clusters are T-cell characteristic. Can we subset further?
# 2) Are clusters with similar markers biologically different enough?
# 3) Can we identify other markers to give us more clustering/identity confidence?

