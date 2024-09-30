# Title: Final Q and A for QC of scRNA tutorial
# Author: Amanda Zacharias
# Date: 2024-09-30
# Email: 16amz1@queensu.ca

Relevant link: [hpctraining scRNA QC tutorial](https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md)

1. Report the number of cells left for each sample, and comment on whether the number of cells removed is high or low. Can you give reasons why this number is still not ~12K (which is how many cells were loaded for the experiment)?

There are ~15,000 cells per sample, very similar to before filtering. Only ~1.7K cells were removed overall, so this is not surprising. Perhaps this number exceeds the number of loaded cells because there are still "junk" cells remaining, and more fine-tunign of thresholds may be necessary.

2. After filtering for nGene per cell, you should still observe a small shoulder to the right of the main peak. What might this shoulder represent?

The cells in the right should have a relatively high number of detected genes. These genes may have failed, be dying, dormant, less complex, or larger in size. Given that the shoulder is quite small, it's likely not a major issue.

3. When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant of the plot? What can you say about these cells that have been removed?

There are arguably no cells in the bottom right quadrant. Whatever cells appear to overlap with that region are likely due to the point size being large. The cells that were removed likely did not pass the filtering thresholds due to a low number of genes per number of UMIs ratio.