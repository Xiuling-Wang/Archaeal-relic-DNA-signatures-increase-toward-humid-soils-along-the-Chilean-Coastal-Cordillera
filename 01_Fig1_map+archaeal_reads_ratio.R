## ============================================================
## New Figure 1 = Panel A (map) + Panel B (archaeal read ratio)
## ============================================================
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr)
library(grid)

# Ensure a UTF-8 text locale so accented site names (e.g. "Pan de Azúcar")
# render correctly. Under a "C" locale the multi-byte "ú" splits into ".."
# on both cairo_pdf and png devices. Harmless if the locale is already UTF-8.
try(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)

output_dir <- "/Users/xwang/🐣Academia/🌸Paper🌸/02 🌼Chile Archaea重点关注/Figures"

# -----------------------------------------------------------
# 0. Load data from archived workspace
# -----------------------------------------------------------
load('/tmp/archaea_env.RData')

# df02 = individual samples; df01 = site x dna_type x depth means
# Limit to 0-60 cm (already filtered in the workspace)
df02 <- df02 %>%
  filter(!is.na(arcpercent)) %>%
  mutate(
    site  = factor(site, levels = c("AZ","SG","LC","NA"),
                   labels = c("AZ\nPan de Azúcar",
                              "SG\nSanta Gracia",
                              "LC\nLa Campana",
                              "NA\nNahuelbuta")),
    dna_type = factor(dna_type, levels = c("iDNA","eDNA"))
  )

df01 <- df01 %>%
  filter(!is.na(mean_value)) %>%
  mutate(
    site  = factor(site, levels = c("AZ","SG","LC","NA"),
                   labels = c("AZ\nPan de Azúcar",
                              "SG\nSanta Gracia",
                              "LC\nLa Campana",
                              "NA\nNahuelbuta")),
    dna_type = factor(dna_type, levels = c("iDNA","eDNA"))
  )

# depth factor: shallow -> deep (Y axis goes 0-5 at top, 40-60 at bottom)
depth_levels <- c("0-5","5-10","10-20","20-40","40-60")
df02$depths_cm_re <- factor(df02$depths_cm_re, levels = rev(depth_levels))
df01$depths_cm_re <- factor(df01$depths_cm_re, levels = rev(depth_levels))

DNA_COLS  <- c("iDNA" = "dodgerblue2", "eDNA" = "darkorange")
DNA_LTYPE <- c("iDNA" = "solid",   "eDNA" = "dashed")

# -----------------------------------------------------------
# Panel B - Archaeal read ratio, unified x-axis 0-16 %
# -----------------------------------------------------------
pb <- ggplot(df02,
             aes(x = arcpercent, y = depths_cm_re, colour = dna_type)) +

  # individual replicates
  geom_point(size = 1.9, alpha = 0.78,
             position = position_jitter(width = 0, height = 0.055, seed = 1)) +

  # depth-mean lines
  geom_path(data  = df01,
            aes(x = mean_value, y = depths_cm_re,
                colour = dna_type, linetype = dna_type,
                group  = dna_type),
            linewidth = 0.75,
            lineend = "round") +

  facet_wrap(~ site, nrow = 1, strip.position = "bottom") +

  scale_x_continuous(
    breaks = c(0, 4, 8, 12, 16),
    expand = expansion(mult = c(0.01, 0.02)),
    position = "top"
  ) +
  coord_cartesian(xlim = c(0, 16), clip = "on") +

  scale_colour_manual(
    values = DNA_COLS,
    labels = c("iDNA (intracellular)", "eDNA (extracellular)"),
    name   = "DNA pool"
  ) +
  scale_linetype_manual(
    values = DNA_LTYPE,
    labels = c("iDNA (intracellular)", "eDNA (extracellular)"),
    name   = "DNA pool"
  ) +

  scale_y_discrete(
    labels = rev(depth_levels)   # show human-readable labels
  ) +

  labs(
    x = "Archaeal reads (% of total prokaryotic reads)",
    y = "Soil depth (cm)"
  ) +

  theme_bw(base_size = 10) +
  theme(
    strip.placement     = "outside",
    strip.background    = element_rect(fill = "grey94", colour = "grey70", linewidth = 0.35),
    strip.text          = element_text(face = "bold", size = 8.5,
                                       lineheight = 1.1),
    panel.border        = element_rect(colour = "grey35", fill = NA, linewidth = 0.45),
    panel.grid.major    = element_line(colour = "grey91", linewidth = 0.35),
    panel.grid.minor    = element_blank(),
    panel.spacing       = unit(0.7, "lines"),
    axis.ticks          = element_line(colour = "grey35", linewidth = 0.35),
    axis.text.x         = element_text(size = 8.2, colour = "grey25"),
    axis.text.y         = element_text(size = 8.2, colour = "grey25"),
    axis.title.x        = element_text(size = 9.2, margin = margin(b = 4)),
    axis.title.y        = element_text(size = 9.2, margin = margin(r = 5)),
    legend.position     = "bottom",
    legend.direction    = "horizontal",
    legend.title        = element_text(size = 9, face = "bold"),
    legend.text         = element_text(size = 8.2),
    legend.box.spacing  = unit(1, "pt"),
    legend.key.width    = unit(1.15, "cm"),
    legend.margin       = margin(t = -2, r = 0, b = 0, l = 0),
    plot.margin         = margin(t = 3, r = 5, b = 2, l = 2)
  ) +
  guides(
    colour   = guide_legend(title.position = "left"),
    linetype = guide_legend(title.position = "left")
  )

# -----------------------------------------------------------
# Panel A - Map
# -----------------------------------------------------------
world  <- ne_countries(scale = "medium", returnclass = "sf")
chile  <- ne_countries(country = "Chile", scale = "medium", returnclass = "sf")

# Site data — "Pan de Azúcar" with proper Spanish accent (cairo_pdf + ragg embed the glyph)
sites <- data.frame(
  code    = c("AZ","SG","LC","NA"),
  label   = c("AZ\nPan de Azúcar", "SG\nSanta Gracia",
              "LC\nLa Campana",    "NA\nNahuelbuta"),
  lat     = c(-26.30, -29.75, -33.00, -37.82),
  lon     = c(-70.46, -71.00, -71.03, -73.00),
  colour  = c("#D73027","#FC8D59","#1A9850","#4575B4")
)
sites_sf <- st_as_sf(sites, coords = c("lon","lat"), crs = 4326)

# Main (zoomed) map
p_main <- ggplot() +
  geom_sf(data = world,
          fill = "grey88", colour = "grey62", linewidth = 0.18) +
  geom_sf(data = chile,
          fill = "grey78", colour = "grey45", linewidth = 0.32) +

  # Climate gradient arrow (north to south, in Pacific Ocean margin)
  annotate("segment",
           x = -77, xend = -77, y = -25.8, yend = -38.8,
           arrow = arrow(length = unit(0.12,"cm"), ends = "last",
                         type = "closed"),
           colour = "grey30", linewidth = 0.55) +

  # "Climate gradient" header above arrow
  annotate("text", x = -77, y = -24.5,
           label = "Climate\ngradient", hjust = 0.5, size = 2.3,
           colour = "grey20", fontface = "italic", lineheight = 0.9) +

  # Dotted line connecting sites
  geom_sf(data = sites_sf, colour = "grey30", size = 0, show.legend = FALSE) +
  annotate("path",
           x   = sites$lon,
           y   = sites$lat,
           colour = "grey35", linetype = "dotted", linewidth = 0.65) +

  # Site points
  geom_point(data = sites,
             aes(x = lon, y = lat, fill = code),
             colour = "white", shape = 21, size = 4.2, stroke = 0.65,
             show.legend = FALSE) +
  scale_fill_manual(values = setNames(sites$colour, sites$code)) +

  # Site labels (right of dots)
  geom_text(data = sites,
            aes(x = lon + 0.8, y = lat, label = label),
            hjust = 0, vjust = 0.5, size = 2.9,
            lineheight = 0.9, fontface = "plain") +

  # Climate zone labels: LEFT-aligned starting just RIGHT of arrow
  # so they never extend past the left xlim boundary
  annotate("text",
           x = -76.5,
           y = c(-26.0, -29.75, -33.0, -37.82),
           label = c("Hyperarid", "Semi-arid", "Mediterranean", "Humid-\ntemperate"),
           hjust = 0, vjust = 0.5, size = 2.1,
           colour = "grey25", lineheight = 0.85) +

  # Minimal cartographic cues for the location panel.
  annotate("segment",
           x = -63.4, xend = -63.4, y = -25.4, yend = -24.15,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           colour = "grey25", linewidth = 0.45) +
  annotate("text", x = -63.4, y = -23.75,
           label = "N", size = 2.4, fontface = "bold", colour = "grey20") +
  annotate("segment",
           x = -78.25, xend = -74.75, y = -40.65, yend = -40.65,
           colour = "grey20", linewidth = 0.55) +
  annotate("segment",
           x = c(-78.25, -74.75), xend = c(-78.25, -74.75),
           y = -40.85, yend = -40.45,
           colour = "grey20", linewidth = 0.55) +
  annotate("text", x = -76.5, y = -40.15,
           label = "300 km", size = 2.0, colour = "grey20") +

  coord_sf(
    xlim   = c(-79, -62),
    ylim   = c(-41.5, -23.0),
    expand = FALSE
  ) +
  scale_x_continuous(breaks = c(-76, -72, -68, -64),
                     labels = function(x) parse(text = paste0(abs(x), "*degree*W"))) +
  scale_y_continuous(breaks = seq(-40, -24, by = 4),
                     labels = function(y) parse(text = paste0(abs(y), "*degree*S"))) +

  theme_bw(base_size = 10) +
  theme(
    axis.text        = element_text(size = 7.6, colour = "grey25"),
    axis.ticks       = element_line(colour = "grey35", linewidth = 0.35),
    axis.title       = element_blank(),
    panel.border     = element_rect(colour = "grey25", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.28),
    panel.background = element_rect(fill = "#D6EAF8"),
    plot.margin      = margin(1, 2, 2, 1)
  )

# Inset: South America - placed in UPPER-RIGHT corner (east of sites, no overlap)
p_inset <- ggplot() +
  geom_sf(data = world, fill = "grey85", colour = "grey55", linewidth = 0.15) +
  geom_sf(data = chile, fill = "grey60", colour = "grey40", linewidth = 0.2) +
  # Study area box
  annotate("rect",
           xmin = -77.5, xmax = -64.5, ymin = -41.0, ymax = -23.5,
           fill = NA, colour = "#E74C3C", linewidth = 0.75) +
  coord_sf(xlim = c(-83, -33), ylim = c(-58, 13), expand = FALSE) +
  theme_void() +
  theme(
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.45),
    panel.background = element_rect(fill = "#D6EAF8")
  )

# Combine map panels: inset in bottom-right corner (below LC/south of all site labels)
pa <- p_main +
  annotation_custom(
    grob   = ggplotGrob(p_inset),
    xmin   = -66.5, xmax = -62,
    ymin   = -41.5, ymax = -34.5
  )

# -----------------------------------------------------------
# Combine A + B with patchwork
# -----------------------------------------------------------
fig1 <- pa + pb +
  plot_layout(ncol = 2, widths = c(1.05, 2.05)) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(plot.tag = element_text(face = "bold", size = 11.5))
  )

# -----------------------------------------------------------
# Save — wider for more strip space in panel B
# -----------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(output_dir, "01_Fig1_redrawn.pdf")
out_png <- file.path(output_dir, "01_Fig1_redrawn.png")
ggsave(out, fig1,
       width = 22.5, height = 11.2, units = "cm",
       device = cairo_pdf)
tmp_png <- file.path(tempdir(), "01_Fig1_redrawn.png")
ggsave(tmp_png, fig1,
       width = 22.5, height = 11.2, units = "cm",
       dpi = 600, bg = "white")
file.copy(tmp_png, out_png, overwrite = TRUE)
cat("Saved to", out, "\n")
