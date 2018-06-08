library(tidyverse)

select_data <- function(compartment, minEHK) {
  rna_seq_df <- read_tsv("./data/rna_seq_tubes.tsv") %>%
    filter(eisenberg_score >= 7) %>%
    filter(tube_comp == "tumor")
  
  responder_df <- read_tsv("./data/responder_status.txt")
  
  rna_seq_df <- rna_seq_df %>%
    left_join(responder_df, by = "sample") %>%
    filter(`responder status` %in% c("NR", "R"))
  
  counts <- read_tsv("./data/counts.tsv")
  counts_cols <- intersect(colnames(counts), rna_seq_df$tube_name)
  counts <- counts %>% dplyr::select(gene_id, counts_cols)
  counts_matrix <- as.matrix(counts[, 2:ncol(counts)])
  rownames(counts_matrix) <- counts$gene_id
  
  responder <- rna_seq_df %>% filter(tube_name %in% counts_cols) %>% filter(`responder status` == "R") %>% dplyr::select(tube_name) %>% unlist
  nonresponder <- rna_seq_df %>% filter(tube_name %in% counts_cols) %>% filter(`responder status` == "NR") %>% dplyr::select(tube_name) %>% unlist
  
  list(
    minEHK = minEHK,
    responder = responder,
    nonresponder = nonresponder,
    counts = counts_matrix,
    compartment = compartment,
    date = date()
  )
}