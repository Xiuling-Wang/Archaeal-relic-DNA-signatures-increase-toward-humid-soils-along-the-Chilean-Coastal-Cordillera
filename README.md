# Chile Archaea Scripts

This folder is organized as the active analysis/code index for the Chile Archaea manuscript.

## Active scripts

| Order | Script | Purpose |
|---:|---|---|
| 00 | `00_config.R` | Shared paths, colours, site/depth levels, and small helper functions. |
| 01 | `01_Fig1_map+archaeal_reads_ratio.R` | Fig. 1 map plus archaeal reads ratio by depth and DNA pool. |
| 02 | `02_Fig2_alpha_diversity_boxplot.R` | Fig. 2 alpha diversity boxplots. Replaces the older alpha script and avoids `introdataviz`. |
| 03 | `03_Fig3_taxonomic_composition_bubbleplots.R` | Fig. 3 / supplementary top taxa bubble plots from phylum to genus. Replaces repeated rank-specific blocks. |
| 04 | `04_Fig4_dbRDA+environment_partitioning.R` | Fig. 4 dbRDA ordinations and environmental partitioning. |
| 05 | `05_Fig5_specialist+generalist_bubbleplots.R` | Fig. 5 specialist/generalist indicator ASV bubble plots. |
| 06 | `06_Fig6_WGCNA_network+module_heatmaps.R` | Fig. 6 WGCNA networks, module tables, and module-trait heatmaps. |
| 07 | `07_FigS1_iDNA+eDNA_Venn.R` | Supplementary iDNA/eDNA Venn diagrams. |
| 08 | `08_FigS2_top_class_barplot.R` | Supplementary top-class composition barplots. |
| 09 | `09_FigS5_top50_ASV_bubbleplot.R` | Supplementary top 50 ASV bubble plot. |
| 10 | `10_FigS2-S3_top_class+order_heatmap.R` | Supplementary top class/order heatmaps. |
| 11 | `11_NMDS+PERMANOVA.R` | NMDS ordination and PERMANOVA checks. |
| 12 | `12_LEfSe_input_prepare.R` | LEfSe input table preparation. |
| 13 | `13_submission_statistics.R` | Submission/statistical summary tables and text outputs. |

## Legacy folders

| Folder | Contents |
|---|---|
| `99_legacy_superseded/` | Older scripts replaced by cleaner active scripts, kept only for traceability. |
| `99_legacy_non_archaea/` | Scripts that point to the bacteria project or older non-archaea workflows. |

## Notes

- Active scripts were checked with `parse(file = ...)`.
- `02_Fig2_alpha_diversity_boxplot.R` and `03_Fig3_taxonomic_composition_bubbleplots.R` were run successfully after cleanup.
- Some advanced scripts still require specialist packages (`WGCNA`, `CorLevelPlot`, `ggraph`, `tidygraph`, `VennDiagram`) when their analyses are rerun.
