library(tidyverse)
library(ggplot2)
library(gplots)
library(edgeR)

p_val_hist <- function(pVals) {
  function() hist(pVals, main = "Histogram of P-Values", xlab = "p-value")
}

volcano_plot <- function(tt) {
  ggplot(tt, aes(x = logFC, y = -log10(P.Value), label=label)) +
    geom_point(data = subset(tt, P.Value > 0.05), 
               color = "blue", alpha = 0.5) +
    geom_point(data = subset(tt, (P.Value <= 0.05 & P.Value > 0.01)),
               color = "red", alpha = 0.5) +
    geom_text(aes(label = ifelse(P.Value <= 0.01, as.character(label),'')),
              hjust = 0, vjust=  0) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme_linedraw() + 
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) +
    xlab("Fold change (log2)") +
    ylab("-log10(P-Value)") +
    ggtitle("Volcano Plot log2FC vs -log10(P-Value)")
}

tmm_norm <- function(counts) {
  # remove rows with all 0
  rows <- apply(counts, 1, function(row) sum(row != 0) != 0)
  filtered_counts <- counts[rows, ]
  
  dge <- DGEList(counts = filtered_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <-
    cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = TRUE,
      prior.count = 0.25
    )
  logCPM
}

hclust.ward = function(d) hclust(d,method="ward.D2")

heatmap <- function(counts, colors) {
  heatmap.2(
    counts,
    trace = "none",
    breaks = (-40:40 / 20),
    hclust = hclust.ward,
    col = colorRampPalette(c("green", "black", "red"))(80),
    ColSideColors = colors,
    srtCol = 45,
    margins = c(9, 8)
  )
  
  legend("topleft",
         legend = c("Responder", "Non Responder"),
         col = c("orange", "purple"),
         lty= 1,
         lwd = 10
  )
}

scale_counts <- function(counts) {
  t(
    scale(
      t(counts)
    )
  )
}

de_heatmaps <- function(data) {
  heatmaps <- list()
  
  norm <- tmm_norm(data$counts)
  tt <- data$de$tt %>% arrange(P.Value)
  norm <- norm[tt$ensembl, c(data$nonresponder, data$responder)]
  rownames(norm) <- tt$label
  
  color_labels <- c(rep("purple", length(data$nonresponder)), rep("orange", length(data$responder)))
  
  if (sum(tt$P.Value <= 0.05) > 300) {
    heatmaps$top300 <- function() heatmap(scale_counts(norm[1:300, ]), color_labels)
  }
  
  sig05Labels <- tt %>% filter(P.Value <= 0.05) %>% dplyr::select(label) %>% unlist
  scaled_sig05 <- scale_counts(norm[sig05Labels, ])
  heatmaps$scaled_sig05_matrix <- scaled_sig05
  heatmaps$sig05 <- function() heatmap(scaled_sig05, color_labels)
  
  sig01Labels <- tt %>% filter(P.Value <= 0.01) %>% dplyr::select(label) %>% unlist
  heatmaps$sig01 <- function() heatmap(scale_counts(norm[sig01Labels, ]), color_labels)
  
  data$heatmaps <- heatmaps
  data
}

de_plots <- function(data) {
  data <- de_heatmaps(data)
  
  data$P.Val_hist <- p_val_hist(data$de$tt$P.Value)
  data$volcano <- volcano_plot(data$de$tt)
  
  data
}