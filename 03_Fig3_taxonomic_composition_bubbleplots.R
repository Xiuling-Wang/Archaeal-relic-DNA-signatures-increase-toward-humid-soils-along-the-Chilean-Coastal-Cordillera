## Fig. 3 / supplementary taxonomy - top archaeal taxa bubble plots.

source("00_config.R")
suppressPackageStartupMessages({
  library(RColorBrewer)
})

asv_path <- file.path(DATA_DIR, "file_A.txt")
env_path <- file.path(DATA_DIR, "env_2024.csv")
out_dir <- ensure_dir(file.path(REPORT_DIR, "bubble plot"))

tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
tax_cols <- c("Kingdom", tax_ranks)

asv_raw <- read.delim(asv_path, header = TRUE, row.names = 1, check.names = FALSE)
env_raw <- read.csv(env_path, row.names = 1, header = TRUE, check.names = FALSE)

asv <- asv_raw %>%
  filter(Kingdom == "Archaea", !is.na(Phylum), Phylum != "")

tax <- asv[, tax_cols, drop = FALSE]
asv_counts <- asv[, setdiff(colnames(asv), tax_cols), drop = FALSE]

env <- env_raw %>%
  filter(row.names(.) %in% colnames(asv_counts), depth < 60) %>%
  mutate(group_col = paste(site_dna, depth_label, sep = "_"))

asv_counts <- asv_counts[, row.names(env), drop = FALSE]
asv_counts <- asv_counts[rowSums(asv_counts) > 0, colSums(asv_counts) > 200, drop = FALSE]
tax <- tax[rownames(asv_counts), , drop = FALSE]
env <- env[row.names(env) %in% colnames(asv_counts), , drop = FALSE]

group_info <- tibble(sample_id = row.names(env), group = env$group_col)
sample_levels <- expand_grid(
  site = SITE_LEVELS,
  dna_type = DNA_LEVELS,
  depth = DEPTH_60_LABELS
) %>%
  mutate(sample_id = paste(site, dna_type, str_replace_all(depth, "-", "_"), sep = "_")) %>%
  pull(sample_id)

relative_abundance <- asv_counts %>%
  rownames_to_column("asv_id") %>%
  pivot_longer(-asv_id, names_to = "sample_id", values_to = "reads") %>%
  left_join(group_info, by = "sample_id") %>%
  group_by(asv_id, group) %>%
  summarise(reads = mean(reads), .groups = "drop") %>%
  group_by(group) %>%
  mutate(relative_abundance = 100 * reads / sum(reads)) %>%
  ungroup() %>%
  select(asv_id, sample_id = group, relative_abundance) %>%
  left_join(tax %>% rownames_to_column("asv_id"), by = "asv_id")

make_rank_plot <- function(rank, top_n = 12) {
  rank_data <- relative_abundance %>%
    filter(!is.na(.data[[rank]]), .data[[rank]] != paste0(str_sub(rank, 1, 1), "__NA")) %>%
    group_by(across(all_of(c("Phylum", rank))), sample_id) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
    group_by(across(all_of(rank))) %>%
    mutate(mean_abundance = mean(relative_abundance)) %>%
    ungroup()

  top_taxa <- rank_data %>%
    distinct(.data[[rank]], mean_abundance) %>%
    slice_max(mean_abundance, n = top_n, with_ties = FALSE) %>%
    pull(.data[[rank]])

  plot_data <- rank_data %>%
    filter(.data[[rank]] %in% top_taxa) %>%
    filter(relative_abundance > 0) %>%
    mutate(
      sample_id = factor(sample_id, levels = sample_levels),
      sample_label = str_replace(as.character(sample_id), "^NB", "NA"),
      taxon = fct_reorder(.data[[rank]], mean_abundance)
    )

  pal <- brewer.pal(max(3, min(8, n_distinct(plot_data$Phylum))), "Set2")

  ggplot(plot_data, aes(x = sample_id, y = taxon)) +
    geom_point(aes(size = relative_abundance, fill = Phylum),
               shape = 21, alpha = 0.95, colour = "grey25", stroke = 0.18) +
    scale_x_discrete(labels = function(x) str_replace(x, "^NB", "NA"), drop = FALSE) +
    scale_fill_manual(values = setNames(pal, sort(unique(plot_data$Phylum)))) +
    scale_size_continuous(range = c(0.5, 8), breaks = c(1, 10, 30, 60, 100), limits = c(0, 100)) +
    labs(x = "Sample", y = NULL, title = paste("Top", rank), size = "Relative abundance (%)") +
    theme_clean(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

plots <- setNames(lapply(tax_ranks, make_rank_plot), tax_ranks)

iwalk(plots, function(plot, rank) {
  ggsave(file.path(out_dir, paste0("Fig3_top_", tolower(rank), "_bubbleplot.pdf")),
         plot, width = 14, height = 4.8, units = "in", device = cairo_pdf)
})
