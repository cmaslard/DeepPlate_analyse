library(dplyr)
library(tidyverse)
library(DESeq2)

get_genes <- function(deseq_contrast, condition_stress, condition_control, direction,lfc_cutoff=1,padj_cutoff_i=0.05) {
  if (deseq_contrast == "climat_condition") {
    dds = dds_climat_condition
  } else {
    dds = dds_condition
  }
  
  results(dds, contrast = c(deseq_contrast, condition_stress, condition_control)) %>%
    as.data.frame() %>%
    filter(padj < padj_cutoff_i) %>%
    filter(abs(log2FoldChange) > lfc_cutoff) %>%
    rownames_to_column("gene") %>%
    filter(
      if (direction == "Up") log2FoldChange > 0
      else if (direction == "Down") log2FoldChange < 0
      else TRUE  # For both directions
    ) %>%
    pull(gene)
}