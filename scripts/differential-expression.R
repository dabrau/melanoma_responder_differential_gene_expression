library(tidyverse)
library(edgeR)
library(limma)
library(methods)
library(statmod)

DE <- function(counts, nonresponder, responder) {
  reordered_counts <- counts[, c(nonresponder, responder)]
  groups <- factor(
    c(
      rep("nonresponder", length(nonresponder)),
      rep("responder", length(responder))
    ),
    levels = c("nonresponder", "responder")
  )
  design <- model.matrix(~groups)
  colnames(design) <- levels(groups)
  
  dge <- DGEList(counts = reordered_counts)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  
  list(
    tt = topTable(fit, number = Inf, coef=ncol(design)),
    fit = fit
  )
}

de_responder_samples <- function(data) {
  de <- DE(data$counts, data$nonresponder, data$responder)
  
  hugo_names <- read_tsv("./data/ensembl-hugo-mapping.txt")
  de$tt <- de$tt %>%
    rownames_to_column("ensembl") %>%
    left_join(hugo_names, by = "ensembl") %>%
    mutate(label = coalesce(hugo, ensembl))
  
  data$de <- de
  data
}