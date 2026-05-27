# Script Inventory

This repository contains the active Chile Archaea analysis scripts copied from the local working script folder and cleaned for public release.

## Active scripts

| Order | Script | Main output or purpose |
|---:|---|---|
| 00 | `R/00_config.R` | Shared repository-relative paths, site/depth levels, colors, and helper functions. |
| 01 | `R/01_Fig1_map+archaeal_reads_ratio.R` | Fig. 1 map and archaeal read ratio by depth and DNA pool. |
| 02 | `R/02_Fig2_alpha_diversity_boxplot.R` | Fig. 2 alpha diversity boxplots. |
| 03 | `R/03_Fig3_taxonomic_composition_bubbleplots.R` | Fig. 3 and supplementary taxonomic bubble plots from phylum to genus. |
| 04 | `R/04_Fig4_dbRDA+environment_partitioning.R` | Fig. 4 dbRDA ordinations and environmental partitioning. |
| 05 | `R/05_Fig5_specialist+generalist_bubbleplots.R` | Fig. 5 specialist/generalist indicator ASV bubble plots. |
| 06 | `R/06_Fig6_WGCNA_network+module_heatmaps.R` | Fig. 6 WGCNA networks, module tables, and module-trait heatmaps. |
| 07 | `R/07_FigS1_iDNA+eDNA_Venn.R` | Supplementary iDNA/eDNA Venn diagrams. |
| 08 | `R/08_FigS2_top_class_barplot.R` | Supplementary top-class composition barplots. |
| 09 | `R/09_FigS5_top50_ASV_bubbleplot.R` | Supplementary top 50 ASV bubble plot. |
| 10 | `R/10_FigS2-S3_top_class+order_heatmap.R` | Supplementary top class/order heatmaps. |
| 11 | `R/11_NMDS+PERMANOVA.R` | NMDS ordination and PERMANOVA checks. |
| 12 | `R/12_LEfSe_input_prepare.R` | LEfSe input table preparation. |
| 13 | `R/13_submission_statistics.R` | Submission/statistical summary tables and text outputs. |
| 14 | `R/14_rarefaction_curve.R` | Supplementary rarefaction curve for retained archaeal samples. |

## Public-release cleanup

- Local absolute paths were replaced with repository-relative path variables defined in `R/00_config.R`.
- R session files, generated plots, and local binary outputs are excluded by `.gitignore`.
- Legacy superseded and non-archaea scripts are not included in this public repository snapshot.
- All active scripts pass `parse(file = ...)`; full execution requires the processed input tables listed in each script.
