library(tidyverse)

de_summary <- function(data) {
  data.frame(
    info = c(
      "date",
      "compartment",
      "min EHK",
      "n",
      "n-responders",
      "n-nonresponders",
      "less than p-value 0.05",
      "less than p-value 0.01"
    ),
    value = c(
      data$date,
      data$compartment,
      as.character(data$minEHK),
      as.character(length(data$responder) + length(data$nonresponder)),
      as.character(length(data$responder)),
      as.character(length(data$nonresponder)),
      as.character(sum(data$de$tt$P.Value <= 0.05)),
      as.character(sum(data$de$tt$P.Value <= 0.01))
    )
  )
}

de_samples <- function(data) {
  data.frame(
    samples = c(
      data$responder,
      data$nonresponder
    ),
    responder_status = c(
      rep("responder", length(data$responder)),
      rep("nonresponder", length(data$nonresponder))
    )
  )
}

matrix_to_tsv <- function(mat, rowname, dir) {
  as.data.frame(mat) %>%
    rownames_to_column(rowname) %>%
    write_tsv(dir)
}

save_plot <- function(report_plot, dir, filename) {
  filepath = paste0(dir, "/", filename)
  if (first(class(report_plot)) == "gg") {
    ggsave(filepath, plot = report_plot, device = "png")
    return ()
  }
  
  png(filename = filepath)
  if (class(report_plot) == "function") {
    report_plot()
    dev.off()
  } else if (class(plot) == "histogram") {
    plot(report_plot)
    dev.off()
  }
}

report <- function(data) {
  output_dir <- paste0(getwd(), "/outputs")
  output_folder <- paste0(
    data$compartment,
    "-EHK",
    data$minEHK,
    "-",
    str_replace_all(
      str_replace_all(data$date, " ", "_"),
      ":", "-"
    ))
  
  output <- paste0(output_dir, "/", output_folder)
  dir.create(output)
  
  data %>% de_summary %>% write_tsv(paste0(output, "/summary.tsv"))
  data %>% de_samples %>% write_tsv(paste0(output, "/samples.tsv"))
  data$counts %>% matrix_to_tsv("gene_id", paste0(output, "/counts.tsv"))
  
  plots_output <- paste0(output, "/plots")
  dir.create(plots_output)
  save_plot(data$P.Val_hist, plots_output, "pval_hist.png")
  save_plot(data$volcano, plots_output, "volcano.png")
  
  for (name in names(data$heatmaps)) {
    item <- data$heatmaps[[name]]
    if (class(item) == "matrix") {
      item %>% matrix_to_tsv("gene_id", paste0(plots_output, "/heatmap_", name, ".tsv"))
    } else if (class(item) == "function") {
      save_plot(item, plots_output, paste0("heatmap_", name, ".png"))
    }
  }
}