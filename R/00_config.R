## Shared configuration for Chile Archaea analysis scripts.
## Set CHILE_ARCHAEA_ROOT to override the repository root when running scripts elsewhere.

suppressPackageStartupMessages({
  library(tidyverse)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = FALSE))
} else {
  normalizePath(getwd(), mustWork = FALSE)
}

default_root <- if (file.exists(file.path(SCRIPT_DIR, "R", "00_config.R"))) {
  SCRIPT_DIR
} else {
  normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
}
PROJECT_ROOT <- Sys.getenv("CHILE_ARCHAEA_ROOT", unset = default_root)
SCRIPT_DIR <- file.path(PROJECT_ROOT, "R")
DATA_DIR <- file.path(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
REPORT_DIR <- file.path(PROJECT_ROOT, "results")
PAPER_FIG_DIR <- file.path(REPORT_DIR, "figures")
TABLE_DIR <- file.path(REPORT_DIR, "tables")
ICLOUD_DATA_DIR <- DATA_DIR
ICLOUD_REPORT_DIR <- REPORT_DIR

SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
SITE_LABELS <- c("AZ", "SG", "LC", "NA")
DNA_LEVELS <- c("iDNA", "eDNA")
DNA_COLORS <- c("iDNA" = "dodgerblue2", "eDNA" = "darkorange")

DEPTH_LEVELS_RAW <- c(
  "12_0_5", "11_5_10", "10_10_20", "09_20_40", "08_40_60",
  "07_60_80", "06_80_100", "05_100_120", "04_120_140",
  "03_140_160", "02_160_180", "01_180_200"
)
DEPTH_LABELS <- c(
  "0-5", "5-10", "10-20", "20-40", "40-60",
  "60-80", "80-100", "100-120", "120-140",
  "140-160", "160-180", "180-200"
)
DEPTH_60_LABELS <- c("0-5", "5-10", "10-20", "20-40", "40-60")

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_tsv_rownames <- function(path) {
  read.table(path, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
}

theme_clean <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.45),
      axis.text = element_text(colour = "grey20"),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
}

format_site_depth <- function(data) {
  data %>%
    mutate(
      site = factor(site, levels = SITE_LEVELS, labels = SITE_LABELS),
      dna_type = factor(dna_type, levels = DNA_LEVELS),
      depths_cm_re = factor(depths_cm_re, levels = DEPTH_LEVELS_RAW, labels = DEPTH_LABELS),
      depth_order = factor(depth_order, levels = as.character(seq_along(DEPTH_LABELS)), labels = paste(DEPTH_LABELS, "cm"))
    )
}
