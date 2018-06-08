library(tidyverse)
setwd("~/Projects/ipi_melanoma/responder_de_pipeline")
source("./scripts/select-rna-samples.R")
source("./scripts/differential-expression.R")
source("./scripts/plots.R")
source("./scripts/report.R")

run <- function(minEHK) {
  compartments <- read_tsv("./data/rna_seq_tubes.tsv") %>%
    filter(eisenberg_score >= minEHK) %>%
    dplyr::select(tube_comp) %>%
    unique %>%
    unlist
  
  for (comp in compartments) {
    data <- select_data(comp, minEHK)
    
    if (length(data$responder) >= 2 & length(data$nonresponder) >= 2) {
      de_responder_samples(data) %>%
        de_plots %>%
        report
    } else print(paste("not enough samples for EHK", minEHK, "and compartment", comp))
  }
}