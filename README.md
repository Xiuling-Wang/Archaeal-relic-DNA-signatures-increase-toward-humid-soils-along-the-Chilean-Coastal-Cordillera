# Chile archaea DNA-pool analyses

Code repository for the manuscript:

**DNA-pool-resolved archaeal communities along a Chilean climate and soil-depth gradient**

Authors: Xiuling Wang and Dirk Wagner

This repository is intended to support the archaeal community analyses reported in the manuscript. The study separates intracellular DNA (iDNA; cell-associated, putatively active fraction) and extracellular DNA (eDNA; extracellular or relic fraction) across four sites along the Chilean Coastal Cordillera climate gradient.

## Data availability

Demultiplexed 16S rRNA gene sequences are deposited in the European Nucleotide Archive:

- ENA project: [PRJEB73502](https://www.ebi.ac.uk/ena/browser/view/PRJEB73502)

The processed archaeal ASV table, sample metadata, supplementary tables, and analysis/visualization code for the archaeal analyses will be added here before manuscript submission. The current `R/` folder contains the cleaned active analysis scripts, with local absolute paths replaced by repository-relative configuration in `R/00_config.R`.

## Repository structure

```text
R/                 R scripts for statistical analyses and figure generation
scripts/           Helper scripts for data checks or format conversion
metadata/          Sample metadata and variable descriptions
data/raw/          Local raw/intermediate data, not tracked by Git
data/processed/    Lightweight processed tables intended for reproducibility
results/figures/   Generated figure outputs, not tracked unless selected
results/tables/    Generated result tables
docs/              Notes on workflow, variables, and reproducibility
```

## Active script workflow

1. Configure paths and shared plotting helpers with `R/00_config.R`.
2. Generate Fig. 1 map and archaeal-read ratio summaries with `R/01_Fig1_map+archaeal_reads_ratio.R`.
3. Generate alpha diversity figures with `R/02_Fig2_alpha_diversity_boxplot.R`.
4. Generate taxonomic composition figures and supplementary plots with `R/03_*`, `R/08_*`, `R/09_*`, and `R/10_*`.
5. Run dbRDA, NMDS, PERMANOVA, and environmental partitioning with `R/04_*` and `R/11_*`.
6. Identify generalist and specialist ASVs with `R/05_Fig5_specialist+generalist_bubbleplots.R`.
7. Build WGCNA co-occurrence modules and module-trait heatmaps with `R/06_Fig6_WGCNA_network+module_heatmaps.R`.
8. Generate supplementary Venn, LEfSe input, submission statistics, and rarefaction outputs with `R/07_*`, `R/12_*`, `R/13_*`, and `R/14_*`.

## Reproducibility notes

- Large raw sequencing outputs and local intermediate files should remain outside Git.
- Processed tables required to reproduce manuscript figures should be placed in `data/processed/` or `metadata/`.
- Scripts should use relative paths from the repository root whenever possible.
- To run scripts outside the repository root, set `CHILE_ARCHAEA_ROOT=/path/to/chile-archaea-dna-pools`.
- Record package versions in `docs/session-info.md` before submission.
