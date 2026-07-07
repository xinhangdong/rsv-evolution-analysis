# rsv-evolution-analysis

Genomic epidemiology and evolutionary surveillance of Respiratory Syncytial Virus (RSV).

**Author:** Xinhang (Gary) Dong — The Webb Schools
**Mentor:** Prof. Sean Wang, Scripps Research
**Status:** Phase 1 (setup and skills), started July 2026

## Question
RSV vaccines and the monoclonal antibody nirsevimab arrived in 2023 and, for the first time,
put immune pressure on the virus. Is RSV beginning to show signs of immune-escape or
resistance evolution? This project builds an open, reproducible surveillance pipeline over
public RSV whole-genome sequences from NCBI to find out.

## Approach (planned)
1. Collect RSV-A / RSV-B whole genomes and metadata from NCBI.
2. Subtype and lineage typing (Nextclade); multiple sequence alignment (MAFFT).
3. Phylogenetics: maximum-likelihood and time-scaled trees; SNP-distance clustering.
4. Gene-level analysis of the G and F proteins, including vaccine / antibody escape mutations.

## Repository layout
- `data/` — raw and processed sequence data (large files are not committed)
- `scripts/` — analysis scripts
- `notebooks/` — Jupyter notebooks
- `results/` — figures and tables

## Reproducing the environment
```
conda env create -f environment.yml
conda activate rsv
```

*This README is a work in progress. I will expand each section in my own words as the project develops.*
