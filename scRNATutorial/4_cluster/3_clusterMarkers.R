#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Identify markers for individual clusters
# Author: Amanda Zacharias
# Date: 2024-10-20
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/09_merged_SC_marker_identification.md
# Goals: determine gene markers for clusters, identify cell types for each cluster using markers,
#         determine if re-clustering is necessary.
# Challenges: Over-interpretation of results, combining different types of marker identification
# Recommedations: Think of results has hypotheses that need verification. Inflated p-values happen
#   b.c. each cell is considered a replicate. Top markers are most trustworthy.
#   Identify markers conserved between conditions, identify markers DE between clusters.
#
# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(ggplot2) # 3.5.1
library(dplyr) # 1.1.4
library(purrr) # 1.0.2; used for `map()`

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
## Input ===========
annoPath <- file.path(projDir, "0_data", "annotation.csv")

## Output ===========
rDataDir <- file.path(projDir, "4_cluster", "rDataDir")
plotsDir <- file.path(projDir, "4_cluster", "plots")
dir.create(rDataDir, showWarnings = FALSE)
dir.create(plotsDir, showWarnings = FALSE)

# Load data -----------------------------------------------
seuratIntegrated <- readRDS(file.path(rDataDir, "seuratIntegratedHbc.rds"))
annotations <- read.csv(annoPath)

# Identify markers for each cluster -----------------------------------------------
# Recommended for when there is one group/condition.
# Wilcoxon Rank Sum Test for differential expression analysis.
# FindAllMarkers() important params:
#   logfc.threshold: min logfc, default = 0.25.
#     could miss markers expressed in fraction of cells in a cluster but not in other clusters.
#     could return metabolic/ribosomal genes, which aren't useful for distinguishing cell identities
#   min.diff.pct: min % diff between % of cells expressing the gene in cluster for all other clusters combined.
#     could miss markers expressed in all cells, but highly up in a specific cell type
#   min.pct: only test genes in min fraction of cells. i.e. don't test lowly expressed genes default: 0.1.
#     could get false negatives if threshold is too high
# Option to get both positive and negative markers. Hbc recommends just getting positive.
# markersWilcoxon <- FindAllMarkers(seuratIntegrated, only.pos = TRUE, logfc.threshold = 0.25) # nolint
# Not actually running this code b.c. we have both ctrl and stim conditions.

# Identification of conserved markers in all conditions -----------------------------------------------
# Better if we have multiple conditions.
# Seprates cells by sample condition, performs de testing for each cluster against all other clusters.
# P-values for each condition are combined across groups using `MetaDE` R package.
## Default Assay != SCTransform ===========
DefaultAssay(seuratIntegrated) <- "RNA"
# RNA has all genes, SCTransform only as top 30000. Also, With SCT, DE analysis would be performed on
# pearson residuals for NB regression, which would be hard to interpret.

## Try one cluster ===========
cluster0_conservedMarkers <- FindConservedMarkers(
  seuratIntegrated,
  ident.1 = 0,
  grouping.var = "sample",
  only.pos = TRUE,
  logfc.threshold = 0.25
)
# Install presto for faster implementation
head(cluster0_conservedMarkers, n = 3)
#      stim_p_val stim_avg_log2FC stim_pct.1 stim_pct.2 stim_p_val_adj
# CCR7          0        1.279817      0.927      0.430              0
# SELL          0        1.402343      0.832      0.370              0
# LDHB          0        1.444768      0.732      0.304              0
#         ctrl_p_val ctrl_avg_log2FC ctrl_pct.1 ctrl_pct.2 ctrl_p_val_adj
# CCR7  0.000000e+00        1.318506      0.839      0.368   0.000000e+00
# SELL  0.000000e+00        1.856410      0.610      0.186   0.000000e+00
# LDHB 2.056415e-293        1.207480      0.734      0.356  2.892347e-289
#           max_pval minimump_p_val
# CCR7  0.000000e+00              0
# SELL  0.000000e+00              0
# LDHB 2.056415e-293              0
#
# positive logfc = gene is more highly expressed in the cluster of interest.
# pct.1 = % of cells where the gene is detected in the cluster
# pct.2 = % of cells where the gene is detected on avg. in other clusters
# max_pval = largest p-value calculated by each condition
# minimump_p_val = combined p-value.
# Recall: p-values will be inflated because each cell is a "replicate"!
# Recommendation: Look for markers w/ lrg diffs between pct1. and pct.2, and larger fold changes.

## Adding gene annotations ===========
cluster0_ann_markers <- cluster0_conservedMarkers %>%
  rownames_to_column(var = "gene") %>%
  left_join(
    y = unique(annotations[, c("gene_name", "description")]),
    by = c("gene" = "gene_name")
  )

## Try cluster 10 ===========
cluster10_conservedMarkers <- FindConservedMarkers(
  seuratIntegrated,
  ident.1 = 10,
  grouping.var = "sample",
  only.pos = TRUE,
  logfc.threshold = 0.25
) %>%
  # Add annotations
  tibble::rownames_to_column(var = "gene") %>%
  left_join(
    y = unique(annotations[, c("gene_name", "description")]),
    by = c("gene" = "gene_name")
  )
# See if FCGR3A and MS4A7 are top genes for cluster 10
cluster10_conservedMarkers %>%
  subset((gene == "FCGR3A") | (gene == "MS4A7"))
# 1 FCGR3A          0        3.971008      0.988      0.101              0
# 3  MS4A7          0        3.242592      0.995      0.159              0
#   ctrl_p_val ctrl_avg_log2FC ctrl_pct.1 ctrl_pct.2 ctrl_p_val_adj max_pval
# 1          0        3.365744      0.980      0.141              0        0
# 3          0        3.471399      0.962      0.122              0        0
#   minimump_p_val
# 1              0
# 3              0
#                                                           description
# 1 Fc fragment of IgG receptor IIIa [Source:HGNC Symbol;Acc:HGNC:3619]
# 3  membrane spanning 4-domains A7 [Source:HGNC Symbol;Acc:HGNC:13378]

# The genes that we expected to be good markers do seem to be good!

## Running on multiple samples ===========
GetConserved <- function(clusterNumber) {
  #' Identify conserved clusters, add gene annotations, and add cluster number
  #'
  #' @param clusterNumber The number of cluster of interest. Integer.
  #' @return Returns a dataframe
  #' @example
  tmp_conservedMarkers <- FindConservedMarkers(
    seuratIntegrated, ident.1 = clusterNumber, grouping.var = "sample", only.pos = TRUE
  ) %>%
    tibble::rownames_to_column(var = "gene") %>%
    left_join(
      y = unique(annotations[, c("gene_name", "description")]),
      by = c("gene" = "gene_name")
    ) %>%
    cbind("cluster_id" = clusterNumber, .)
  return(tmp_conservedMarkers)
}
# Execute
conservedMarkers <- map_dfr(c(4, 0, 6, 2), GetConserved) # only doing 4 clusters for tutorial
# Some clusters may not have enough cells for ctrl or stim,
# in which case you need the FindAllMarkers() function instead.

# Evalulating marker genes -----------------------------------------------
## Get top 10 markers per cluster ===========
top10 <- conservedMarkers %>%
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
  top_n(n = 10, wt = avg_fc)
View(top10)

# Both cluster 0 and 6 have CCR7 and SELL, markers of T-cells.
# However, cluster 2 has CREM, which is a marker of activation.
# So it might make sense to keep them separate instead of merging.

# Cluster 4 has many heat shock and DNA damage genes. Cells are likely
# stressed / dying. Quality metrics like mitoRatio and nuMi don't support this.
# So, this could contain reactive T-cells ?

## Visualize marker genes ===========
# Feature plot
cluster4_markerGenesPlot <- FeaturePlot(
  seuratIntegrated,
  features = c("HSPH1", "HSPE1", "DNAJB1"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE,
  repel = TRUE
)
ggsave(
  plot = cluster4_markerGenesPlot,
  filename = "cluster4_markerGenesPlot.pdf", path = plotsDir,
  width = 180, height = 180, units = "mm"
)

# Violin plot
cluster4_violinPlot <- VlnPlot(
  object = seuratIntegrated,
  features = c("HSPH1", "HSPE1", "DNAJB1")
)
ggsave(
  plot = cluster4_violinPlot,
  filename = "cluster4_violinPlot.pdf", path = plotsDir,
  width = 180, height = 90, units = "mm"
)

# Refining marker genes for each cluster -----------------------------------------------
# The list of markers may not be sufficient to separate clusters.
# Use FindMarkers() to differentiate two specific clusters; e.g., multiple t-cell like clusters.
## Determine markers for CD4+ T cells ===========
cd4Tcells <- FindMarkers(seuratIntegrated, ident.1 = 2, ident.2 = c(0, 4, 6)) %>%
  # Add annotations
  tibble::rownames_to_column("gene") %>%
  left_join(
     y = unique(annotations[, c("gene_name", "description")]),
     by = c("gene" = "gene_name")
  ) %>%
  # Reorder columns
  dplyr::select(c(1, 3:5, 2, 6:7)) %>%
  # Sort
  arrange(p_val_adj)
# Inspect
View(cd4Tcells)

# CREM has positive fold change as marker of activation
# SELL and CCR7 negative fold change marker of naive / memory cells
# Based on results, hbc authors estimated cell type labels for the different clusters

## Reassign identity ===========
seuratIntegrated <- RenameIdents(
  object = seuratIntegrated,
  "0" = "Naive or memory CD4+ T cells",
  "1" = "CD14+ monocytes",
  "2" = "Activated T cells",
  "3" = "CD14+ monocytes",
  "4" = "Stressed cells / Unknown",
  "5" = "CD8+ T cells",
  "6" = "Naive or memory CD4+ T cells",
  "7" = "B cells",
  "8" = "NK cells",
  "9" = "CD8+ T cells",
  "10" = "FCGR3A+ monocytes",
  "11" = "B cells",
  "12" = "NK cells",
  "13" = "B cells",
  "14" = "Conventional dendritic cells",
  "15" = "Megakaryocytes",
  "16" = "Plasmacytoid dendritic cells"
)
## Plot UMAP ===========
umapPlot_renamedIdents <- DimPlot(
  object = seuratIntegrated,
  reduction = "umap",
  label = TRUE,
  label.size = 3,
  repel = TRUE
)
ggsave(
  plot = umapPlot_renamedIdents,
  filename = "umapPlot_renamedIdents.pdf", path = plotsDir,
  width = 270, height = 180, units = "mm"
)

## Plot w/o funky cells ===========
# Remove cluster 4 with markers for damage and stress
seuratSubsetLabeled <- subset(seuratIntegrated, idents = "Stressed cells / Unknown", invert = TRUE)
umapPlot_renamedIdentsSubset <- DimPlot(
  object = seuratSubsetLabeled,
  reduction = "umap",
  label = TRUE,
  label.size = 3,
  repel = TRUE
)
ggsave(
  plot = umapPlot_renamedIdentsSubset,
  filename = "umapPlot_renamedIdentsSubset.pdf", path = plotsDir,
  width = 270, height = 180, units = "mm"
)

# Save -----------------------------------------------
saveRDS(seuratIntegrated, file = file.path(rDataDir, "finalSeuratIntegrated.rds"))
sink(file.path(projDir, "sessionInfo_scRNAseq_20241020.txt"))
sessionInfo()
sink()
