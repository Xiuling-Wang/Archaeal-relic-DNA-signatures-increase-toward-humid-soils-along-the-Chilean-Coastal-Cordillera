from pathlib import Path


REQUIRED_PATHS = [
    "README.md",
    "docs/variables.md",
    "docs/session-info.md",
    "docs/script-inventory.md",
    "R/00_config.R",
    "R/01_Fig1_map+archaeal_reads_ratio.R",
    "R/02_Fig2_alpha_diversity_boxplot.R",
    "R/03_Fig3_taxonomic_composition_bubbleplots.R",
    "R/04_Fig4_dbRDA+environment_partitioning.R",
    "R/05_Fig5_specialist+generalist_bubbleplots.R",
    "R/06_Fig6_WGCNA_network+module_heatmaps.R",
    "R/07_FigS1_iDNA+eDNA_Venn.R",
    "R/08_FigS2_top_class_barplot.R",
    "R/09_FigS5_top50_ASV_bubbleplot.R",
    "R/10_FigS2-S3_top_class+order_heatmap.R",
    "R/11_NMDS+PERMANOVA.R",
    "R/12_LEfSe_input_prepare.R",
    "R/13_submission_statistics.R",
    "R/14_rarefaction_curve.R",
]


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    missing = [path for path in REQUIRED_PATHS if not (root / path).exists()]
    if missing:
        raise SystemExit("Missing required files:\n" + "\n".join(missing))
    print("Repository scaffold looks complete.")


if __name__ == "__main__":
    main()
