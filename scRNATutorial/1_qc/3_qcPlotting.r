#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Generate quality control plots before filtering
# Author: Amanda Zacharias
# Date: 2024-09-30
# Email: 16amz1@queensu.ca
#-----------------------------------------------
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.4.0
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md

# Run this script before and after filtering to generate QC plots.
# For the sake of reducing code, pathws at top are simply changed with the "prefix" option.

# **We will not check for doublets.**
# Doublets: one sample generated from two cells. Incorrectly suggest
#   intermediate populations or transitory states.
# Why not check?
#   1) Traditional high thresholds for UMIs or genes seem logical but != accurate.
#   2) Existing tools falsely remove true signals. Scrublet is a popular tool to be explored.
# Once clusters are identified, see if a cell-marker applies to multiple cell types.

# Options -----------------------------------------------
# Is this plotting for before and after QC filtering?
prefix <- "aft" # b4 OR aft

# Packages -----------------------------------------------
library(Seurat) # 5.1.0
library(dplyr) # 1.1.4
library(ggplot2) # 3.5.1

# Source -----------------------------------------------

# Pathways -----------------------------------------------
projDir <- file.path("scRNATutorial")
# Input ===========
rDataDir <- file.path(projDir, "1_qc", "rDataDir")

# Output ===========
plotsDir <- file.path(projDir, "1_qc", "plots")
plotsB4FiltDir <- file.path(plotsDir, prefix)
dir.create(plotsDir, showWarnings = FALSE)
dir.create(plotsB4FiltDir, showWarnings = FALSE)

# Load data -----------------------------------------------
if (prefix == "b4") {
  mergedSeurat <- readRDS(file = file.path(rDataDir, "mergedSeuratQcMetrics.rds"))
  metadata <- readRDS(file = file.path(rDataDir, "metadata.rds"))
} else if (prefix == "aft") {
  mergedSeurat <- readRDS(file = file.path(rDataDir, "filteredSeurat.rds"))
  metadata <- readRDS(file = file.path(rDataDir, "filteredMetadata.rds"))
}
metadata <- metadata %>%
  mutate(sample = factor(sample, levels = c("ctrl", "stim")))

# Visualizations -----------------------------------------------
# Cell counts ===========
# Different platforms have different capture rates, so be thoughtful about your thresholds.
# It's also possible to get way more "captured" cells than expected because of dying cells
#   & other technical issues.
# For this experiment, we expect 12,000 to 13,000 cells.
cellCountsPlot <- metadata %>%
  ggplot(aes(x = sample, fill = sample)) +
  geom_bar() +
  ggtitle("n cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(
  plot = cellCountsPlot,
  filename = "cellCountsPlot.pdf", path = plotsB4FiltDir,
  width = 90, height = 90, units = "mm"
)

# UMI (transcripts) counts per cell ===========
# Generally should at least 500.
# If counts are 500 - 1000, they're usable but more depth is ideal.
umiCountsPlot <- metadata %>%
  ggplot(aes(x = nUMI, colour = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 500) +
  scale_x_log10() +
  theme_bw() +
  ylab("cell density")
ggsave(
  plot = umiCountsPlot,
  filename = "umiCountsPlot.pdf", path = plotsB4FiltDir,
  width = 90, height = 90, units = "mm"
)

# Genes per Cell ===========
# Similar pattern expected for gene counts per cell
# Expect a **single large peak of encapsulated cells**.
# A small peak to left or bimodal peak could mean...
#   1) set of cells failed
#   2) biologically dif cells, e.g., inactive or less complex cells
#   3) cells are smaller
geneCountsPlot <- metadata %>%
  ggplot(aes(x = nGene, colour = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 300) +
  scale_x_log10() +
  theme_bw() +
  ylab("cell density")
ggsave(
  plot = geneCountsPlot,
  filename = "geneCountsPlot.pdf", path = plotsB4FiltDir,
  width = 90, height = 90, units = "mm"
)

# Complexity ===========
# Novely score we calculated: nGene / UMI
# If high UMI but low nGene --> the same genes were sequenced many times
# Low complexity/novelty could rep a single cell type (e.g. red blood cells), or
#   an artefact / contamination.
# A novelty score > 0.80 is generally good.
complexityPlot <- metadata %>%
  ggplot(aes(x = log10GenesPerUMI, colour = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 0.8) +
  theme_bw()
ggsave(
  plot = complexityPlot,
  filename = "complexityPlot.pdf", path = plotsB4FiltDir,
  width = 90, height = 90, units = "mm"
)

# Mitochondrial counts ratio ===========
# Is there mitochondrial contamination from dead/dying cells
# Cells with ratio > 0.2 are generally poor quality unless
#   you expect otherwise.
mitoPlot <- metadata %>%
  ggplot(aes(x = mitoRatio, colour = sample, fill = sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  geom_vline(xintercept = 0.2) +
  theme_bw()
ggsave(
  plot = mitoPlot,
  filename = "mitoPlot.pdf", path = plotsB4FiltDir,
  width = 90, height = 90, units = "mm"
)

# Reads per cell ===========
# Could be worth exploring, but not in this tutorial.
#   Because the workflow used would need to save this information to assess.
# Hope to se all samples with peaks in ~same location between 10,000 and 100,000.

# Joint filtering effects ===========
# Considering metrics in isolation can lead to misinterpretation.
#   Ex. cells with high mito counts may be involved in respiratory processes & not be dying.
# Therefore, thresholds for indiv metrics should be as permissive as possible.
#   Consider joint effects where possible.
# Number of UMIs and genes per cell are often considered jointly.
# Here: # UMIs per cell vs # genes per cell coloured by fraction of mito reads
countsXMitoPlot <- metadata %>%
  ggplot(aes(x = nUMI, y = nGene, colour = mitoRatio)) +
  geom_point() +
  stat_smooth(method = lm) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  scale_x_log10() +
  scale_y_log10() +
  scale_x_continuous(name = "# UMIs per cell") +
  scale_y_continuous(name = "# genes per cell") +
  scale_colour_gradient(name = "mito ratio", low = "grey90", high = "black") +
  facet_wrap(~sample) +
  theme_bw()
ggsave(
  plot = countsXMitoPlot,
  filename = "countsXMitoPlot.pdf", path = plotsB4FiltDir,
  width = 180, height = 90, units = "mm"
)

# Upper right quadrant = generally high quality
# Lower left quadrant = generally low quality
# Bottom right quadrant = high UMIs but low genes, so dying or low complexity
# Also evalulate the slope of the line.

# High mitoRatio cells are in extra low count cells w/ few genes, which could
# indicate dying cells rather than respiratory/biologically interesting cells.
