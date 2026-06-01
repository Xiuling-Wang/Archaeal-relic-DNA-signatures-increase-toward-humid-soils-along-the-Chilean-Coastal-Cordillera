library(tidyverse)
library(vegan)
library(rdacca.hp)
library(cowplot)
library(grid)
library(png)

base_dir <- "/Users/xwang/🐣Academia"
data_dir <- file.path(base_dir, "data/02 🌼Chile Archaea")
report_dir <- file.path(base_dir, "report/02 🌼Chile Archaea/Network")
fig_dir <- file.path(base_dir, "🌸Paper🌸/02 🌼Chile Archaea重点关注/Figures/FigS")
tmp_dir <- "/tmp/archaea_figs_env_format"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

env_vars <- c("pH", "Conductivity", "moisture", "CN", "Feo", "Alo", "Mno", "Sio", "NH4", "NO3", "Po", "Pi")

pretty_labels <- c(
  depth = "Depth",
  pH = "pH",
  Conductivity = "Conductivity",
  moisture = "Moisture",
  N = "N",
  C = "C",
  CN = "C:N",
  Feo = "Fe[ox]",
  Fed = "Fe[d]",
  FeoFed = "Fe[ox]:Fe[d]",
  Fep = "Fe[p]",
  Alo = "Al[ox]",
  Alp = "Al[p]",
  Sio = "Si[ox]",
  NH4 = "NH[4]^'+'",
  NO3 = "NO[3]^'-'",
  Pt = "P[t]",
  Pi = "P[i]",
  Po = "P[o]",
  Mnd = "Mn[d]",
  Mno = "Mn[ox]",
  Mnp = "Mn[p]"
)

parse_labeller <- function(x) {
  out <- unname(pretty_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  parse(text = out)
}

sig_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

format_r <- function(x) {
  out <- sprintf("%.2f", x)
  out <- sub("0+$", "", out)
  sub("\\.$", "", out)
}

theme_pub <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
      axis.text = element_text(color = "#4A4A4A"),
      axis.title = element_text(color = "black"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

# Fig. S7: hierarchical partitioning ---------------------------------------
asv_raw <- read.table(file.path(data_dir, "arc_60cm_202rare/asv_202.txt"),
                      header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
env <- read.csv(file.path(data_dir, "env_2024.csv"), header = TRUE, check.names = FALSE)
env_60 <- env %>%
  filter(depth < 60, sample_id %in% colnames(asv_raw)) %>%
  arrange(match(sample_id, colnames(asv_raw)))
asv <- asv_raw[, env_60$sample_id, drop = FALSE]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0, drop = FALSE]
env_60 <- env_60 %>% filter(sample_id %in% colnames(asv))
asv <- asv[, env_60$sample_id, drop = FALSE]

make_hp_df <- function(pool) {
  env_pool <- env_60 %>%
    filter(dna_type == pool) %>%
    select(all_of(env_vars), sample_id) %>%
    na.omit()
  asv_pool <- asv[, env_pool$sample_id, drop = FALSE]
  asv_pool <- asv_pool[rowSums(asv_pool) > 0, , drop = FALSE]
  hp <- rdacca.hp(
    vegdist(t(asv_pool)),
    env_pool %>% select(all_of(env_vars)),
    method = "dbRDA",
    type = "adjR2"
  )
  hp$Hier.part %>%
    as.data.frame() %>%
    rownames_to_column("Variable") %>%
    transmute(Variable, Individual_percent = `I.perc(%)`, pool)
}

hp_df <- bind_rows(make_hp_df("iDNA"), make_hp_df("eDNA")) %>%
  group_by(pool) %>%
  arrange(desc(Individual_percent), .by_group = TRUE) %>%
  mutate(Variable = factor(Variable, levels = unique(Variable))) %>%
  ungroup()

make_hp_panel <- function(pool, label) {
  ggplot(hp_df %>% filter(pool == !!pool),
         aes(x = reorder(Variable, -Individual_percent), y = Individual_percent)) +
    geom_col(fill = "#5B5B5B", width = 0.9) +
    scale_x_discrete(labels = parse_labeller) +
    labs(x = "Variables", y = expression("% Individual effect to " * R^2 * " (%)"),
         tag = label) +
    theme_pub(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.tag = element_text(size = 24, face = "bold")
    )
}

p_s7 <- plot_grid(
  make_hp_panel("iDNA", "A"),
  make_hp_panel("eDNA", "B"),
  ncol = 1,
  align = "v"
)

ggsave(file.path(fig_dir, "07_FigS7_env_hp_formatted.pdf"),
       p_s7, width = 8.6, height = 13.1, units = "in", device = cairo_pdf)

# Fig. S9B: eDNA module-environment heatmap --------------------------------
r_mat <- read.csv(file.path(report_dir, "eDNA_Module_trait_r.csv"), row.names = 1, check.names = FALSE)
p_mat <- read.csv(file.path(report_dir, "eDNA_Module_trait_p-value.csv"), row.names = 1, check.names = FALSE)

module_order <- c("MEred", "MEblack", "MEpink", "MEmagenta", "MEbrown",
                  "MEturquoise", "MEyellow", "MEblue", "MEgreen")
trait_order <- c("depth", "pH", "Conductivity", "moisture", "N", "C", "CN",
                 "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                 "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")

r_mat <- r_mat[module_order[module_order %in% rownames(r_mat)],
               trait_order[trait_order %in% colnames(r_mat)], drop = FALSE]
p_mat <- p_mat[rownames(r_mat), colnames(r_mat), drop = FALSE]

heat_df <- as.data.frame(as.table(as.matrix(r_mat))) %>%
  rename(Module = Var1, Trait = Var2, r = Freq) %>%
  mutate(
    p = as.vector(as.matrix(p_mat)),
    Module = factor(Module, levels = rev(rownames(r_mat))),
    Trait = factor(Trait, levels = colnames(r_mat)),
    label = paste0(format_r(r), sig_stars(p))
  )

p_heat <- ggplot(heat_df, aes(x = Trait, y = Module, fill = r)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.25) +
  geom_text(aes(label = label), size = 3.35) +
  scale_fill_gradientn(
    colours = c("blue1", "skyblue", "white", "pink", "red"),
    values = scales::rescale(c(-0.9, -0.45, 0, 0.45, 0.9)),
    limits = c(-0.9, 0.9),
    breaks = seq(-0.8, 0.8, by = 0.4),
    name = NULL,
    guide = guide_colorbar(
      barheight = unit(58, "mm"),
      barwidth = unit(6, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  scale_x_discrete(labels = parse_labeller, position = "bottom") +
  scale_y_discrete(labels = function(x) sub("^ME", "ME", x)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    legend.position = "right",
    legend.text = element_text(color = "black", size = 10),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 3, 3, 3)
  )

network_png <- file.path(tmp_dir, "eDNA_network.png")
if (!file.exists(network_png)) {
  stop("Missing network PNG: run pdftoppm on eDNA_co-ocurrence network.pdf to create ", network_png)
}
network_img <- png::readPNG(network_png)
network_grob <- rasterGrob(network_img, interpolate = TRUE)
p_network <- ggdraw() +
  draw_grob(network_grob) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             size = 28, fontface = "bold")

p_s9 <- plot_grid(
  p_network,
  ggdraw(p_heat) + draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
                              size = 28, fontface = "bold"),
  ncol = 1,
  rel_heights = c(1.08, 0.88)
)

ggsave(file.path(fig_dir, "09_FigS9_eDNA_co-ocurrence_network_formatted_v2.pdf"),
       p_s9, width = 15.2, height = 14.4, units = "in", device = cairo_pdf)

# Fig. S9 original-style label overlay -------------------------------------
# This version preserves the original assembled Fig. S9 exactly and only redraws
# the heatmap x-axis labels with plotmath formatting.
original_s9_pdf <- file.path(fig_dir, "09_FigS9_eDNA_co-ocurrence network.pdf")
original_s9_png <- file.path(tmp_dir, "original_s9_300.png")
if (!file.exists(original_s9_png)) {
  original_s9_tmp_pdf <- file.path(tmp_dir, "original_s9.pdf")
  file.copy(original_s9_pdf, original_s9_tmp_pdf, overwrite = TRUE)
  system2("pdftoppm", c("-png", "-singlefile", "-r", "300",
                        shQuote(original_s9_tmp_pdf),
                        shQuote(file.path(tmp_dir, "original_s9_300"))))
}

original_s9_img <- png::readPNG(original_s9_png)
original_s9_grob <- rasterGrob(original_s9_img, interpolate = TRUE)

overlay_trait_order <- c("depth", "pH", "Conductivity", "moisture", "N", "C", "CN",
                         "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                         "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")
overlay_labels <- parse_labeller(overlay_trait_order)
left_edge <- 113 / 1613
right_edge <- 1516 / 1613
cell_width <- (right_edge - left_edge) / length(overlay_trait_order)
label_x <- left_edge + ((seq_along(overlay_trait_order) - 0.5) * cell_width)

p_s9_original_style <- ggdraw() +
  draw_grob(original_s9_grob) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)),
            x = 0, y = 0, width = 1, height = 0.058)

for (i in seq_along(overlay_trait_order)) {
  p_s9_original_style <- p_s9_original_style +
    draw_label(overlay_labels[i], x = label_x[i], y = 0.030,
               hjust = 0.5, vjust = 0.5, size = 13, color = "black")
}

ggsave(file.path(fig_dir, "09_FigS9_eDNA_co-ocurrence_network_original_style_labels.pdf"),
       p_s9_original_style, width = 17.9167, height = 14.4057,
       units = "in", device = cairo_pdf)

# Fig. S9 v2: redraw eDNA heatmap with site-module grouping ----------------
edna_site_cor <- tribble(
  ~Module,        ~AZ,    ~SG,    ~LC,    ~NB,
  "MEyellow",    -0.31,  -0.31,  -0.30,   0.92,
  "MEturquoise",  0.89,  -0.04,  -0.42,  -0.43,
  "MEred",       -0.38,   0.01,   0.72,  -0.35,
  "MEpink",      -0.21,   0.55,  -0.14,  -0.19,
  "MEmagenta",   -0.24,   0.74,  -0.23,  -0.27,
  "MEgreen",     -0.31,  -0.25,  -0.03,   0.59,
  "MEbrown",     -0.30,   0.79,  -0.18,  -0.31,
  "MEblue",      -0.28,  -0.28,  -0.26,   0.83,
  "MEblack",     -0.23,  -0.21,   0.71,  -0.27
)

edna_module_order <- c("MEred", "MEblack", "MEpink", "MEmagenta", "MEbrown",
                       "MEturquoise", "MEyellow", "MEblue", "MEgreen")
edna_module_labels <- c("red(7)", "black(1)", "pink(6)", "magenta(5)",
                        "brown(3)", "turquoise(8)", "yellow(9)",
                        "blue(2)", "green(4)")

# Fig. S8: iDNA module-site correlation heatmap ----------------------------
idna_site_heatmap_order <- c("MEyellow", "MEturquoise", "MEred", "MEpink",
                             "MEmagenta", "MEgreen", "MEbrown", "MEblue",
                             "MEblack")
idna_site_module_labels <- c("yellow (9)", "turquoise (8)", "red (7)",
                             "pink (6)", "magenta (5)", "green (4)",
                             "brown (3)", "blue (2)", "black (1)")
idna_site_label_mat <- tribble(
  ~Module,        ~AZ,        ~SG,        ~LC,        ~NA_site,
  "MEyellow",    "-0.37",    "-0.23",    "0.89***", "-0.29",
  "MEturquoise", "-0.27",    "-0.3",     "-0.31",   "0.87***",
  "MEred",       "-0.28",    "-0.24",    "-0.17",   "0.69***",
  "MEpink",      "-0.25",    "0.75***",  "-0.2",    "-0.3",
  "MEmagenta",   "0.82***",  "-0.15",    "-0.34",   "-0.33",
  "MEgreen",     "-0.26",    "0.83***",  "-0.26",   "-0.31",
  "MEbrown",     "-0.2",     "0.59**",   "-0.18",   "-0.21",
  "MEblue",      "-0.17",    "-0.28",    "-0.25",   "0.7***",
  "MEblack",     "0.58**",   "-0.23",    "-0.2",    "-0.15"
)

idna_site_heat_df <- idna_site_label_mat %>%
  pivot_longer(c(AZ, SG, LC, NA_site), names_to = "Site", values_to = "label") %>%
  mutate(
    Module = factor(Module, levels = rev(idna_site_heatmap_order)),
    Site = recode(Site, NA_site = "NA"),
    Site = factor(Site, levels = c("AZ", "SG", "LC", "NA")),
    module_label = factor(
      idna_site_module_labels[match(as.character(Module), idna_site_heatmap_order)],
      levels = rev(idna_site_module_labels)
    ),
    r = as.numeric(gsub("\\*+", "", label))
  )

idna_site_label_mat_for_export <- idna_site_label_mat
names(idna_site_label_mat_for_export)[names(idna_site_label_mat_for_export) == "NA_site"] <- "NA"
write.csv(idna_site_label_mat_for_export,
          file.path(report_dir, "iDNA_sites_modules_heatmap_labels_for_FigS8.csv"),
          row.names = FALSE)

p_idna_site_modules <- ggplot(idna_site_heat_df,
                              aes(x = Site, y = module_label, fill = r)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.35, width = 1, height = 1) +
  geom_text(aes(label = label), size = 3.8, color = "black") +
  scale_fill_gradientn(
    colours = c("blue1", "skyblue", "white", "pink", "red"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation",
    guide = guide_colorbar(
      barheight = unit(47, "mm"),
      barwidth = unit(6, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 12, face = "bold",
                               margin = margin(t = 6)),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    axis.ticks = element_blank(),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(file.path(fig_dir, "08_FigS8_sites_iDNAmodules_heatmap_redrawn.pdf"),
       p_idna_site_modules, width = 7.2, height = 3.2,
       units = "in", device = cairo_pdf)

# Fig. S10: eDNA module-site correlation heatmap ---------------------------
# Kept separate from the iDNA module-site heatmap because module colours and
# numbers are network-specific and should not be interpreted across pools.
edna_site_heatmap_order <- c("MEyellow", "MEturquoise", "MEred", "MEpink",
                             "MEmagenta", "MEgreen", "MEbrown", "MEblue",
                             "MEblack")
edna_site_module_labels <- c("yellow (9)", "turquoise (8)", "red (7)",
                             "pink (6)", "magenta (5)", "green (4)",
                             "brown (3)", "blue (2)", "black (1)")
edna_site_label_mat <- tribble(
  ~Module,        ~AZ,        ~SG,       ~LC,        ~NA_site,
  "MEyellow",    "-0.31",    "-0.31",   "-0.3",     "0.92***",
  "MEturquoise", "0.89***",  "-0.04",   "-0.42",    "-0.43",
  "MEred",       "-0.38",    "0.01",    "0.72***",  "-0.35",
  "MEpink",      "-0.21",    "0.55*",   "-0.14",    "-0.19",
  "MEmagenta",   "-0.24",    "0.74***", "-0.23",    "-0.27",
  "MEgreen",     "-0.31",    "-0.25",   "-0.03",    "0.59**",
  "MEbrown",     "-0.3",     "0.79***", "-0.18",    "-0.31",
  "MEblue",      "-0.28",    "-0.28",   "-0.26",    "0.83***",
  "MEblack",     "-0.23",    "-0.21",   "0.71***",  "-0.27"
)

edna_site_heat_df <- edna_site_label_mat %>%
  pivot_longer(c(AZ, SG, LC, NA_site), names_to = "Site", values_to = "label") %>%
  mutate(
    Module = factor(Module, levels = rev(edna_site_heatmap_order)),
    Site = recode(Site, NA_site = "NA"),
    Site = factor(Site, levels = c("AZ", "SG", "LC", "NA")),
    module_label = factor(
      edna_site_module_labels[match(as.character(Module), edna_site_heatmap_order)],
      levels = rev(edna_site_module_labels)
    ),
    r = as.numeric(gsub("\\*+", "", label))
  )

edna_site_label_mat_for_export <- edna_site_label_mat
names(edna_site_label_mat_for_export)[names(edna_site_label_mat_for_export) == "NA_site"] <- "NA"
write.csv(edna_site_label_mat_for_export,
          file.path(report_dir, "eDNA_sites_modules_heatmap_labels_for_FigS10.csv"),
          row.names = FALSE)

p_edna_site_modules <- ggplot(edna_site_heat_df,
                              aes(x = Site, y = module_label, fill = r)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.35, width = 1, height = 1) +
  geom_text(aes(label = label), size = 3.8, color = "black") +
  scale_fill_gradientn(
    colours = c("blue1", "skyblue", "white", "pink", "red"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation",
    guide = guide_colorbar(
      barheight = unit(47, "mm"),
      barwidth = unit(6, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 12, face = "bold",
                               margin = margin(t = 6)),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    axis.ticks = element_blank(),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(file.path(fig_dir, "10_FigS10_sites_eDNAmodules_heatmap.pdf"),
       p_edna_site_modules, width = 7.2, height = 3.2,
       units = "in", device = cairo_pdf)

edna_site_assignments <- edna_site_cor %>%
  pivot_longer(c(AZ, SG, LC, NB), names_to = "site_code", values_to = "r") %>%
  group_by(Module) %>%
  slice_max(r, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    site_label = recode(site_code, NB = "NA"),
    Module = factor(Module, levels = edna_module_order)
  ) %>%
  arrange(Module)

write.csv(edna_site_assignments,
          file.path(report_dir, "eDNA_module_site_assignments_for_FigS9.csv"),
          row.names = FALSE)

edna_groups <- edna_site_assignments %>%
  mutate(row_id = row_number()) %>%
  group_by(site_label) %>%
  summarise(first_row = min(row_id), last_row = max(row_id), .groups = "drop") %>%
  arrange(first_row)

original_s9_img <- png::readPNG(original_s9_png)
original_s9_top <- original_s9_img[1:2500, , , drop = FALSE]
p_s9_top <- ggdraw() +
  draw_grob(rasterGrob(original_s9_top, interpolate = TRUE)) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)),
            x = 0, y = 0, width = 0.09, height = 0.13)

edna_data_trait_order <- c("depth", "pH", "Conductivity", "moisture", "N", "C", "CN",
                           "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                           "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")
edna_display_trait_order <- c("depth", "pH", "Cond.", "moisture", "N", "C", "CN",
                              "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                              "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")

edna_r_mat <- read.csv(file.path(report_dir, "eDNA_Module_trait_r.csv"),
                       row.names = 1, check.names = FALSE)
edna_p_mat <- read.csv(file.path(report_dir, "eDNA_Module_trait_p-value.csv"),
                       row.names = 1, check.names = FALSE)
edna_r_mat <- edna_r_mat[edna_module_order, edna_data_trait_order, drop = FALSE]
edna_p_mat <- edna_p_mat[edna_module_order, colnames(edna_r_mat), drop = FALSE]
colnames(edna_r_mat)[colnames(edna_r_mat) == "Conductivity"] <- "Cond."
colnames(edna_p_mat)[colnames(edna_p_mat) == "Conductivity"] <- "Cond."

edna_heat_df <- as.data.frame(as.table(as.matrix(edna_r_mat))) %>%
  rename(Module = Var1, Trait = Var2, r = Freq) %>%
  mutate(
    p = as.vector(as.matrix(edna_p_mat)),
    row_id = match(Module, edna_module_order),
    y = length(edna_module_order) - row_id + 1,
    x = match(Trait, edna_display_trait_order),
    label = paste0(format_r(r), sig_stars(p))
  )

edna_module_label_df <- tibble(
  Module = edna_module_order,
  module_label = edna_module_labels,
  row_id = seq_along(edna_module_order),
  y = length(edna_module_order) - row_id + 1
)

edna_group_df <- edna_groups %>%
  mutate(
    y_top = length(edna_module_order) - first_row + 1 + 0.5,
    y_bottom = length(edna_module_order) - last_row + 1 - 0.5,
    y_mid = (y_top + y_bottom) / 2
  )

edna_label_strings <- c(
  depth = "Depth", pH = "pH", `Cond.` = "Cond.", moisture = "Moist.",
  N = "N", C = "C", CN = "C:N", Feo = "Fe[ox]", Fed = "Fe[d]",
  FeoFed = "Fe[ox]:Fe[d]", Fep = "Fe[p]", Alo = "Al[ox]",
  Alp = "Al[p]", Sio = "Si[ox]", NH4 = "NH[4]^'+'",
  NO3 = "NO[3]^'-'", Pt = "P[t]", Pi = "P[i]", Po = "P[o]",
  Mnd = "Mn[d]", Mno = "Mn[ox]", Mnp = "Mn[p]"
)
edna_axis_labels <- parse(text = unname(edna_label_strings[edna_display_trait_order]))

p_s9_heat_v2 <- ggplot(edna_heat_df, aes(x = x, y = y, fill = r)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.25, width = 1, height = 1) +
  geom_text(aes(label = label), size = 4.35, color = "black") +
  geom_text(data = edna_module_label_df,
            aes(x = 0.45, y = y, label = module_label),
            inherit.aes = FALSE, hjust = 1, size = 5.0,
            fontface = "bold", color = "black") +
  geom_segment(data = edna_group_df,
               aes(x = -1.70, xend = -1.70, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = edna_group_df,
               aes(x = -1.70, xend = -1.40, y = y_top, yend = y_top),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = edna_group_df,
               aes(x = -1.70, xend = -1.40, y = y_bottom, yend = y_bottom),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = edna_group_df,
            aes(x = -2.18, y = y_mid, label = site_label),
            inherit.aes = FALSE, hjust = 0.5, size = 6.5,
            fontface = "bold", color = "black") +
  scale_fill_gradientn(
    colours = c("blue1", "skyblue", "white", "pink", "red"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation",
    guide = guide_colorbar(
      barheight = unit(76, "mm"),
      barwidth = unit(7, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(edna_display_trait_order),
    labels = edna_axis_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(edna_module_order),
    labels = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(-2.55, length(edna_display_trait_order) + 0.5),
                  ylim = c(0.5, length(edna_module_order) + 0.5),
                  clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 13.5, margin = margin(t = 8)),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(color = "black", size = 15, face = "bold"),
    legend.text = element_text(color = "black", size = 13),
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 14)
  )

p_s9_v2 <- plot_grid(
  p_s9_top,
  ggdraw(p_s9_heat_v2) +
    draw_label("B", x = 0.018, y = 0.98, hjust = 0, vjust = 1,
               size = 32, fontface = "bold"),
  ncol = 1,
  rel_heights = c(0.98, 0.58)
)

ggsave(file.path(fig_dir, "09_FigS9_eDNA_co-ocurrence_network_redrawn_heatmap_v2.pdf"),
       p_s9_v2, width = 17.92, height = 14.4,
       units = "in", device = cairo_pdf)

# Fig. 6 original-style label overlay --------------------------------------
# Preserve the original main Fig. 6 and redraw only the site group labels and
# environmental labels. Site groups are assigned from the module-site
# correlation result generated by the WGCNA analysis.
fig_main_dir <- file.path(base_dir, "🌸Paper🌸/02 🌼Chile Archaea重点关注/Figures")
original_fig6_pdf <- file.path(fig_main_dir, "06_Fig6_iDNA_environment_modules_heatmap.pdf")
original_fig6_png <- file.path(tmp_dir, "original_fig6_300.png")
if (!file.exists(original_fig6_png)) {
  original_fig6_tmp_pdf <- file.path(tmp_dir, "original_fig6.pdf")
  file.copy(original_fig6_pdf, original_fig6_tmp_pdf, overwrite = TRUE)
  system2("pdftoppm", c("-png", "-singlefile", "-r", "300",
                        shQuote(original_fig6_tmp_pdf),
                        shQuote(file.path(tmp_dir, "original_fig6_300"))))
}

fig6_site_cor <- tribble(
  ~Module,       ~AZ,    ~SG,    ~LC,    ~NB,
  "MEyellow",   -0.37,  -0.23,   0.89,  -0.29,
  "MEturquoise",-0.27,  -0.30,  -0.31,   0.87,
  "MEred",      -0.28,  -0.24,  -0.17,   0.69,
  "MEpink",     -0.25,   0.75,  -0.20,  -0.30,
  "MEmagenta",   0.82,  -0.15,  -0.34,  -0.33,
  "MEgreen",    -0.26,   0.83,  -0.26,  -0.31,
  "MEbrown",    -0.20,   0.59,  -0.18,  -0.21,
  "MEblue",     -0.17,  -0.28,  -0.25,   0.70,
  "MEblack",     0.58,  -0.23,  -0.20,  -0.15
)

fig6_module_order <- c("MEmagenta", "MEblack", "MEpink", "MEgreen", "MEbrown",
                       "MEturquoise", "MEblue", "MEred", "MEyellow")
fig6_site_assignments <- fig6_site_cor %>%
  pivot_longer(c(AZ, SG, LC, NB), names_to = "site_code", values_to = "r") %>%
  group_by(Module) %>%
  slice_max(r, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    site_label = recode(site_code, NB = "NA"),
    Module = factor(Module, levels = fig6_module_order)
  ) %>%
  arrange(Module)

write.csv(fig6_site_assignments,
          file.path(report_dir, "iDNA_module_site_assignments_for_Fig6.csv"),
          row.names = FALSE)

fig6_groups <- fig6_site_assignments %>%
  mutate(row_id = row_number()) %>%
  group_by(site_label) %>%
  summarise(first_row = min(row_id), last_row = max(row_id), .groups = "drop") %>%
  mutate(site_label = factor(site_label, levels = c("AZ", "SG", "NA", "LC"))) %>%
  arrange(site_label)

fig6_img <- png::readPNG(original_fig6_png)
fig6_grob <- rasterGrob(fig6_img, interpolate = TRUE)

fig6_trait_order <- c("depth", "pH", "Cond.", "moisture", "N", "C", "CN",
                      "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                      "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")
fig6_labels <- c(
  depth = "Depth", pH = "pH", `Cond.` = "Cond.", moisture = "Moist.",
  N = "N", C = "C", CN = "C:N", Feo = "Fe[ox]", Fed = "Fe[d]",
  FeoFed = "Fe[ox]:Fe[d]", Fep = "Fe[p]", Alo = "Al[ox]",
  Alp = "Al[p]", Sio = "Si[ox]", NH4 = "NH[4]^'+'",
  NO3 = "NO[3]^'-'", Pt = "P[t]", Pi = "P[i]", Po = "P[o]",
  Mnd = "Mn[d]", Mno = "Mn[ox]", Mnp = "Mn[p]"
)
fig6_parsed_labels <- parse(text = unname(fig6_labels[fig6_trait_order]))

fig6_left_edge <- 0.092
fig6_right_edge <- 0.926
fig6_cell_width <- (fig6_right_edge - fig6_left_edge) / length(fig6_trait_order)
fig6_label_x <- fig6_left_edge + ((seq_along(fig6_trait_order) - 0.5) * fig6_cell_width)
fig6_top_row_y <- 0.302
fig6_row_step <- 0.0276

p_fig6_original_style <- ggdraw() +
  draw_grob(fig6_grob) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)),
            x = 0, y = 0, width = 1, height = 0.061) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)),
            x = 0, y = 0.055, width = 0.092, height = 0.270)

for (i in seq_along(fig6_trait_order)) {
  p_fig6_original_style <- p_fig6_original_style +
    draw_label(fig6_parsed_labels[i], x = fig6_label_x[i], y = 0.033,
               hjust = 0.5, vjust = 0.5, size = 14, color = "black")
}

fig6_module_labels <- c("magenta(5)", "black(1)", "pink(6)", "green(4)",
                        "brown(3)", "turquoise(8)", "blue(2)", "red(7)",
                        "yellow(9)")
for (i in seq_along(fig6_module_labels)) {
  p_fig6_original_style <- p_fig6_original_style +
    draw_label(fig6_module_labels[i],
               x = 0.088,
               y = fig6_top_row_y - (i - 1) * fig6_row_step,
               hjust = 1, vjust = 0.5, size = 11.5,
               fontface = "bold", color = "black")
}

for (i in seq_len(nrow(fig6_groups))) {
  y_top <- fig6_top_row_y - (fig6_groups$first_row[i] - 1) * fig6_row_step + 0.014
  y_bottom <- fig6_top_row_y - (fig6_groups$last_row[i] - 1) * fig6_row_step - 0.014
  y_mid <- (y_top + y_bottom) / 2
  p_fig6_original_style <- p_fig6_original_style +
    draw_line(x = c(0.031, 0.031), y = c(y_bottom, y_top), size = 1.1, color = "black") +
    draw_line(x = c(0.031, 0.044), y = c(y_top, y_top), size = 1.1, color = "black") +
    draw_line(x = c(0.031, 0.044), y = c(y_bottom, y_bottom), size = 1.1, color = "black") +
    draw_label(as.character(fig6_groups$site_label[i]), x = 0.013, y = y_mid,
               hjust = 0.5, vjust = 0.5, size = 16, fontface = "bold",
               color = "black")
}

ggsave(file.path(fig_main_dir, "06_Fig6_iDNA_environment_modules_heatmap_original_style_labels.pdf"),
       p_fig6_original_style, width = 18.3922, height = 17.8022,
       units = "in", device = cairo_pdf)

# Fig. 6 v2: redraw the heatmap so all text is readable --------------------
fig6_img_top <- fig6_img[1:3600, , , drop = FALSE]
fig6_top_grob <- rasterGrob(fig6_img_top, interpolate = TRUE)
p_fig6_top <- ggdraw() + draw_grob(fig6_top_grob)

fig6_data_trait_order <- c("depth", "pH", "Conductivity", "moisture", "N", "C", "CN",
                           "Feo", "Fed", "FeoFed", "Fep", "Alo", "Alp", "Sio",
                           "NH4", "NO3", "Pt", "Pi", "Po", "Mnd", "Mno", "Mnp")

# The current Module_trait_r.csv is identical to eDNA_Module_trait_r.csv, so it
# cannot be used for the iDNA main figure. These labels are recovered from the
# original iDNA_environment_modules_heatmap.pdf and keep the original stars.
fig6_idna_labels <- tribble(
  ~Module, ~depth, ~pH, ~Conductivity, ~moisture, ~N, ~C, ~CN, ~Feo, ~Fed, ~FeoFed, ~Fep, ~Alo, ~Alp, ~Sio, ~NH4, ~NO3, ~Pt, ~Pi, ~Po, ~Mnd, ~Mno, ~Mnp,
  "MEmagenta", "-0.3", "0.6**", "0.04", "-0.35", "-0.37", "-0.37", "-0.22", "-0.26", "-0.1", "-0.29", "-0.48*", "-0.36", "-0.41", "-0.61**", "0.09", "0.75***", "-0.05", "0.41", "-0.12", "-0.58**", "-0.47*", "-0.42",
  "MEblack", "0.42", "0.35", "-0.07", "-0.16", "-0.21", "-0.19", "0.53*", "-0.22", "0.35", "-0.24", "-0.17", "-0.19", "-0.15", "-0.45*", "-0.17", "-0.13", "-0.17", "-0.05", "-0.15", "-0.4", "-0.26", "-0.24",
  "MEpink", "-0.36", "0.4", "-0.05", "-0.31", "-0.33", "-0.32", "-0.29", "-0.04", "-0.32", "0.01", "0.08", "-0.27", "-0.21", "0.02", "-0.19", "-0.28", "-0.35", "-0.51*", "-0.26", "-0.25", "-0.34", "-0.32",
  "MEgreen", "-0.03", "0.37", "-0.02", "-0.3", "-0.35", "-0.34", "-0.32", "-0.03", "-0.09", "-0.06", "0.23", "-0.28", "-0.17", "-0.01", "-0.16", "-0.28", "-0.08", "-0.57**", "0.02", "-0.32", "-0.37", "-0.36",
  "MEbrown", "0.39", "0.24", "0.02", "-0.21", "-0.24", "-0.23", "-0.22", "-0.11", "-0.09", "-0.11", "0.08", "-0.2", "-0.14", "0.01", "-0.02", "-0.23", "0.05", "-0.4", "0.11", "-0.24", "-0.26", "-0.25",
  "MEturquoise", "-0.01", "-0.69***", "-0.19", "0.94***", "0.76***", "0.67**", "0.11", "0.82***", "0.09", "0.76***", "0.72***", "0.9***", "0.84***", "0.7***", "0.17", "-0.01", "-0.47*", "0.48*", "-0.53*", "0.73***", "0.8***", "0.72***",
  "MEblue", "0.33", "-0.53*", "-0.28", "0.67**", "0.37", "0.33", "0.11", "0.67**", "0.11", "0.66**", "0.58**", "0.77***", "0.68***", "0.67**", "-0.16", "-0.08", "-0.3", "0.31", "-0.34", "0.4", "0.45*", "0.25",
  "MEred", "-0.36", "-0.58**", "0.02", "0.67**", "0.84***", "0.76***", "0.11", "0.52*", "-0.12", "0.52*", "0.56*", "0.59**", "0.67**", "0.44", "0.85***", "-0.04", "-0.32", "0.51*", "-0.39", "0.73***", "0.79***", "0.88***",
  "MEyellow", "0.15", "-0.26", "-0.1", "-0.26", "-0.17", "-0.18", "-0.16", "-0.46*", "0.18", "-0.42", "-0.45*", "-0.27", "-0.31", "-0.04", "-0.33", "-0.31", "0.68***", "-0.39", "0.72***", "0.19", "-0.03", "-0.03"
)

fig6_heat_df <- fig6_idna_labels %>%
  pivot_longer(all_of(fig6_data_trait_order), names_to = "Trait", values_to = "label") %>%
  mutate(
    row_id = match(Module, fig6_module_order),
    y = length(fig6_module_order) - row_id + 1,
    Trait = if_else(Trait == "Conductivity", "Cond.", Trait),
    x = match(Trait, fig6_trait_order),
    r = as.numeric(gsub("\\*+", "", label))
  )

fig6_r_recovered <- fig6_heat_df %>%
  select(Module, Trait, r) %>%
  pivot_wider(names_from = Trait, values_from = r) %>%
  arrange(match(Module, fig6_module_order))
fig6_label_recovered <- fig6_heat_df %>%
  select(Module, Trait, label) %>%
  pivot_wider(names_from = Trait, values_from = label) %>%
  arrange(match(Module, fig6_module_order))
write.csv(fig6_r_recovered,
          file.path(report_dir, "iDNA_Module_trait_r_recovered_from_original_Fig6.csv"),
          row.names = FALSE)
write.csv(fig6_label_recovered,
          file.path(report_dir, "iDNA_Module_trait_labels_recovered_from_original_Fig6.csv"),
          row.names = FALSE)

fig6_module_label_df <- tibble(
  Module = fig6_module_order,
  module_label = fig6_module_labels,
  row_id = seq_along(fig6_module_order),
  y = length(fig6_module_order) - row_id + 1
)

fig6_group_df <- fig6_groups %>%
  mutate(
    y_top = length(fig6_module_order) - first_row + 1 + 0.5,
    y_bottom = length(fig6_module_order) - last_row + 1 - 0.5,
    y_mid = (y_top + y_bottom) / 2
  )

fig6_axis_df <- tibble(
  Trait = fig6_trait_order,
  x = seq_along(fig6_trait_order)
)

p_fig6_heat_v2 <- ggplot(fig6_heat_df, aes(x = x, y = y, fill = r)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.25, width = 1, height = 1) +
  geom_text(aes(label = label), size = 4.35, color = "black") +
  geom_text(data = fig6_module_label_df,
            aes(x = 0.45, y = y, label = module_label),
            inherit.aes = FALSE, hjust = 1, size = 5.0,
            fontface = "bold", color = "black") +
  geom_segment(data = fig6_group_df,
               aes(x = -1.70, xend = -1.70, y = y_bottom, yend = y_top),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = fig6_group_df,
               aes(x = -1.70, xend = -1.40, y = y_top, yend = y_top),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = fig6_group_df,
               aes(x = -1.70, xend = -1.40, y = y_bottom, yend = y_bottom),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = fig6_group_df,
            aes(x = -2.18, y = y_mid, label = site_label),
            inherit.aes = FALSE, hjust = 0.5, size = 6.5,
            fontface = "bold", color = "black") +
  scale_fill_gradientn(
    colours = c("blue1", "skyblue", "white", "pink", "red"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation",
    guide = guide_colorbar(
      barheight = unit(76, "mm"),
      barwidth = unit(7, "mm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  scale_x_continuous(
    breaks = fig6_axis_df$x,
    labels = fig6_parsed_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(fig6_module_order),
    labels = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(-2.55, length(fig6_trait_order) + 0.5),
                  ylim = c(0.5, length(fig6_module_order) + 0.5),
                  clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, margin = margin(t = 8)),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(color = "black", size = 15, face = "bold"),
    legend.text = element_text(color = "black", size = 13),
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 14)
  )

p_fig6_v2 <- plot_grid(
  p_fig6_top,
  ggdraw(p_fig6_heat_v2),
  ncol = 1,
  rel_heights = c(0.96, 0.48)
)

ggsave(file.path(fig_main_dir, "06_Fig6_iDNA_environment_modules_heatmap_redrawn_heatmap_v3_correct_iDNA.pdf"),
       p_fig6_v2, width = 18.4, height = 17.2,
       units = "in", device = cairo_pdf)
