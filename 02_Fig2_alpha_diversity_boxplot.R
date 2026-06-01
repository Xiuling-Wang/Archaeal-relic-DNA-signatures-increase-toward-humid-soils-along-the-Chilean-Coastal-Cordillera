## Fig. 2 - alpha diversity boxplots for archaeal communities.

source("00_config.R")
suppressPackageStartupMessages({
  library(ggpubr)
  library(patchwork)
  library(vegan)
})

alpha_candidates <- c(
  file.path(REPORT_DIR, "submission_statistics", "alpha_diversity_60cm.csv"),
  file.path(ICLOUD_DATA_DIR, "arc_60cm_202rare", "alpha_diversity.txt")
)
alpha_path <- alpha_candidates[file.exists(alpha_candidates)][1]
env_path <- file.path(DATA_DIR, "env_2024.csv")
asv_path <- file.path(DATA_DIR, "arc_60cm_202rare", "asv_202.txt")
out_dir <- ensure_dir(file.path(REPORT_DIR, "alpha diversity"))

env <- read.csv(env_path, header = TRUE, check.names = FALSE)

if (!is.na(alpha_path)) {
  alpha <- if (str_detect(alpha_path, "\\.csv$")) {
    read.csv(alpha_path, header = TRUE, check.names = FALSE)
  } else {
    read.table(alpha_path, header = TRUE, check.names = FALSE)
  }
} else {
  asv <- read.table(asv_path, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
  asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
  asv_t <- t(asv) %>% as.data.frame()
  alpha <- tibble(
    sample_id = rownames(asv_t),
    Shannon = diversity(asv_t, index = "shannon"),
    Simpson = diversity(asv_t, index = "simpson"),
    Chao1 = as.numeric(estimateR(asv_t)["S.chao1", ])
  )
}

alpha_input <- alpha %>%
  select(sample_id, Shannon, Chao1, Simpson, any_of(c("site", "dna_type", "depth_order", "depth")))
for (col in c("site", "dna_type", "depth_order", "depth")) {
  if (!col %in% names(alpha_input)) alpha_input[[col]] <- NA
}

env_meta <- env %>%
  select(sample_id, site_env = site, dna_type_env = dna_type,
         depth_order_env = depth_order, depth_env = depth)

alpha_data <- alpha_input %>%
  left_join(env_meta, by = "sample_id") %>%
  mutate(
    site = coalesce(as.character(site), as.character(site_env)),
    dna_type = coalesce(as.character(dna_type), as.character(dna_type_env)),
    depth_order = coalesce(as.character(depth_order), as.character(depth_order_env)),
    depth = coalesce(as.numeric(depth), as.numeric(depth_env))
  ) %>%
  filter(depth < 60) %>%
  select(sample_id, site, dna_type, depth_order, Shannon, Chao1, Simpson) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS, labels = SITE_LABELS),
    dna_type = factor(dna_type, levels = DNA_LEVELS),
    depth_order = factor(depth_order, levels = as.character(seq_along(DEPTH_LABELS)),
                         labels = paste(DEPTH_LABELS, "cm"))
  ) %>%
  drop_na(Shannon, Chao1, Simpson)

make_alpha_panel <- function(metric, y_label, tag = NULL) {
  ggplot(alpha_data, aes(x = site, y = .data[[metric]], fill = dna_type)) +
    geom_boxplot(
      width = 0.56,
      alpha = 0.78,
      linewidth = 0.32,
      outlier.shape = NA,
      position = position_dodge(width = 0.68)
    ) +
    geom_point(
      aes(colour = dna_type),
      size = 1.6,
      alpha = 0.72,
      stroke = 0,
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.68, seed = 1)
    ) +
    stat_compare_means(
      aes(group = dna_type),
      method = "anova",
      label = "p.signif",
      size = 3.4,
      hide.ns = FALSE
    ) +
    scale_fill_manual(values = DNA_COLORS, name = "DNA pool") +
    scale_colour_manual(values = DNA_COLORS, guide = "none") +
    labs(x = "Sites", y = y_label, tag = tag) +
    theme_clean(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

p_shannon <- make_alpha_panel("Shannon", "Shannon (H)", "A")
p_chao1 <- make_alpha_panel("Chao1", "Chao1", "B")
p_simpson <- make_alpha_panel("Simpson", "Simpson", "C")

fig2 <- (p_shannon | p_chao1) / p_simpson +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Fig2_alpha_diversity_boxplot.pdf"), fig2,
       width = 9.5, height = 7.2, units = "in", device = cairo_pdf)
ggsave(file.path(PAPER_FIG_DIR, "02_Fig2_redrawn.pdf"), fig2,
       width = 9.5, height = 7.2, units = "in", device = cairo_pdf)
