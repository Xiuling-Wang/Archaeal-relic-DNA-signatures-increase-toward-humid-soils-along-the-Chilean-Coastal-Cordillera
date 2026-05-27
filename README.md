# Chile archaea DNA-pool analyses

Code repository for the manuscript:

**DNA-pool-resolved archaeal communities along a Chilean climate and soil-depth gradient**

Authors: Xiuling Wang and Dirk Wagner

This repository is intended to support the archaeal community analyses reported in the manuscript. The study separates intracellular DNA (iDNA; cell-associated, putatively active fraction) and extracellular DNA (eDNA; extracellular or relic fraction) across four sites along the Chilean Coastal Cordillera climate gradient.

## Data availability

Demultiplexed 16S rRNA gene sequences are deposited in the European Nucleotide Archive:

- ENA project: [PRJEB73502](https://www.ebi.ac.uk/ena/browser/view/PRJEB73502)

The processed archaeal ASV table, sample metadata, supplementary tables, and analysis/visualization code for the archaeal analyses will be added here before manuscript submission.

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

## Planned workflow

1. Prepare archaeal ASV table and metadata.
2. Filter archaeal ASVs and remove negative-control-enriched taxa.
3. Generate alpha diversity summaries and rarefaction sensitivity checks.
4. Run beta diversity, PERMANOVA, PERMDISP, and dbRDA analyses.
5. Identify generalist and specialist ASVs.
6. Build WGCNA co-occurrence modules for iDNA and eDNA pools.
7. Generate publication figures and supplementary tables.

## Reproducibility notes

- Large raw sequencing outputs and local intermediate files should remain outside Git.
- Processed tables required to reproduce manuscript figures should be placed in `data/processed/` or `metadata/`.
- Scripts should use relative paths from the repository root whenever possible.
- Record package versions in `docs/session-info.md` before submission.
