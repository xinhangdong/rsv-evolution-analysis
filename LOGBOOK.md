# Project Logbook — rsv-evolution-analysis

A running record of work, commands, decisions, and results.
Kept continuously from day one (ISEF requires a continuous project logbook).

---

## 2026-07-07 — Day 1: environment setup
- Installed Miniforge (conda) on macOS, Apple Silicon.
- Created conda environment `rsv` (Python 3.11) with Biopython, pandas, numpy, matplotlib, JupyterLab.
- Configured git; set GitHub username to `xinhangdong`.
- Created the project folder structure and this repository; first commit.
- Wrote `scripts/fetch_one_genome.py`; fetched RSV-B (NC_001781): 15,225 bp, molecule RNA, 10 genes (NS1, NS2, N, P, M, SH, G, F, M2, L).
- Pushed the repository to GitHub (public).
- **Step 1 (environment setup) complete.** Next: continue Phase 1 skills; work toward the August milestone — batch-download RSV genomes from NCBI and build metadata.csv (>= 2000 sequences).

## 2026-07-08 — RSV-A vs RSV-B reference genomes
- Switched the script to NC_038235 (RSV-A) and compared it to RSV-B (NC_001781): both have the same 10 genes in the same order and nearly the same length (RSV-A 15,222 bp vs RSV-B 15,225 bp).
- Genes sit at slightly different positions, and the G gene length differs (RSV-A 923 bp vs RSV-B 926 bp) while the N gene is identical (1124-2327) in both — same genome layout, different sequence.

<!-- Add a new dated entry every working session: what I did, commands run, decisions made, results, and what is next. -->
