## Supplementary rarefaction curve — archaeal ASVs
## Fig. S10: shows 202-read threshold is appropriate
## Output: FigS/10_FigS10_rarefaction.pdf + .png

SCRIPT_BASE <- "/Users/xwang/🐣Academia/script/02 🌼Chile Archaea"
source(file.path(SCRIPT_BASE, "00_config.R"), chdir = FALSE, local = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# ── 1. Paths ─────────────────────────────────────────────────────────────────
unrarefied_path <- file.path(DATA_DIR, "arc_unrarefild",
                             "ASV_Arc_200cm_delete_less200.txt")
rarefied_path   <- file.path(DATA_DIR, "arc_60cm_202rare", "asv_202.txt")
fig_dir         <- file.path(PAPER_FIG_DIR, "FigS")
ensure_dir(fig_dir)

RARE_DEPTH <- 202   # rarefaction threshold

# ── 2. Load tables ────────────────────────────────────────────────────────────
raw  <- read_tsv_rownames(unrarefied_path)   # rows = ASVs, cols = samples
rare <- read_tsv_rownames(rarefied_path)     # rows = ASVs, cols = 95 retained

# Keep only the 95 retained samples (0–60 cm), using their names
kept_samples <- colnames(rare)
# Strip any leading/trailing whitespace or BOM artifacts
kept_samples <- trimws(kept_samples)
colnames(raw) <- trimws(colnames(raw))

# Intersect (some unrarefied names may include extra suffix)
use_cols <- intersect(colnames(raw), kept_samples)
cat(sprintf("Retained samples found in unrarefied table: %d / %d\n",
            length(use_cols), length(kept_samples)))

otu <- t(raw[, use_cols])   # vegan wants samples as rows

# ── 3. Sample metadata ────────────────────────────────────────────────────────
# Naming rule:
#   leading 'e' → eDNA ; leading 'i' OR no leading letter → iDNA
#   site codes: PA/PB/PC → AZ | SA/SB/SC → SG | LA/LB/LC → LC | NA/NAC/NB → NA

parse_sample <- function(nm) {
  nm_clean <- sub("_lib.*", "", nm)          # strip _lib suffix
  # DNA type
  dna <- if (grepl("^e", nm_clean)) "eDNA" else "iDNA"
  core <- sub("^[ei]", "", nm_clean)         # remove leading e/i
  # Site
  site <- case_when(
    grepl("^P",  core) ~ "AZ",
    grepl("^S",  core) ~ "SG",
    grepl("^L",  core) ~ "LC",
    grepl("^N",  core) ~ "NA",
    TRUE               ~ "Unknown"
  )
  list(dna_type = dna, site = site)
}

meta <- tibble(sample = use_cols) %>%
  mutate(
    parsed   = map(sample, parse_sample),
    dna_type = map_chr(parsed, "dna_type"),
    site     = map_chr(parsed, "site"),
    site     = factor(site, levels = SITE_LABELS),
    dna_type = factor(dna_type, levels = DNA_LEVELS)
  ) %>%
  select(-parsed)

# ── 4. Compute rarefaction curves (vegan) ─────────────────────────────────────
set.seed(42)
rc <- rarecurve(otu, step = 10, sample = RARE_DEPTH, tidy = TRUE)
# tidy = TRUE returns a data.frame: Site, Sample, Species

# Rename columns for clarity
rc <- rc %>%
  rename(sample = Site, n_reads = Sample, n_asvs = Species)

# Attach metadata
rc <- rc %>% left_join(meta, by = "sample")

# Sample read totals (for labelling excluded samples)
read_totals <- rowSums(otu)

# ── 5. Plot ───────────────────────────────────────────────────────────────────
SITE_COLORS <- c(
  "AZ" = "#D62728",
  "SG" = "#FF7F0E",
  "LC" = "#2CA02C",
  "NA" = "#1F77B4"
)

XMAX <- 1500   # zoom window; captures all samples near threshold

p <- ggplot(rc %>% dplyr::filter(n_reads <= XMAX),
            aes(x = n_reads, y = n_asvs,
                group = sample,
                colour = site,
                linetype = dna_type)) +
  geom_line(alpha = 0.60, linewidth = 0.5) +
  geom_vline(xintercept = RARE_DEPTH,
             colour = "grey30", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = RARE_DEPTH + 15, y = 1,
           label = paste0(RARE_DEPTH, " reads\n(rarefaction depth)"),
           hjust = 0, vjust = 0, size = 2.5, colour = "grey30") +
  scale_colour_manual(values = SITE_COLORS, name = "Site") +
  scale_linetype_manual(values = c("iDNA" = "solid", "eDNA" = "dashed"),
                        name   = "DNA pool") +
  scale_x_continuous(limits = c(0, XMAX),
                     labels = scales::comma) +
  labs(
    x     = "Sequencing depth (reads)",
    y     = "Observed ASVs",
    caption = paste0("All 95 retained samples shown (0–60 cm, n ≥ 202 reads).",
                     " Curves truncated at 1,500 reads for clarity.")
  ) +
  theme_clean(base_size = 9) +
  theme(
    legend.position  = "right",
    legend.key.width = unit(1.2, "lines"),
    plot.caption     = element_text(size = 6.5, colour = "grey40")
  )

# ── 6. Save ───────────────────────────────────────────────────────────────────
out_pdf <- file.path(fig_dir, "10_FigS10_rarefaction.pdf")
out_png <- file.path(fig_dir, "10_FigS10_rarefaction.png")

ggsave(out_pdf, p, width = 12, height = 8, units = "cm")
ggsave(out_png, p, width = 12, height = 8, units = "cm",
       dpi = 300, bg = "white")

cat("Figure saved:\n  ", out_pdf, "\n  ", out_png, "\n")

# ── 7. Quick summary table (for Methods text) ─────────────────────────────────
totals <- tibble(sample = rownames(otu), reads = rowSums(otu)) %>%
  left_join(meta, by = "sample")

cat("\n── Read depth summary by site/pool ──\n")
totals %>%
  group_by(site, dna_type) %>%
  summarise(
    n        = n(),
    min_reads = min(reads),
    median   = median(reads),
    max_reads = max(reads),
    .groups  = "drop"
  ) %>%
  print(n = Inf)

cat(sprintf("\nOverall retained samples: %d\n", nrow(totals)))
cat(sprintf("Min reads: %d | Median reads: %.0f | Max reads: %d\n",
            min(totals$reads), median(totals$reads), max(totals$reads)))
