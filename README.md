# PacBio Structural Variant Pipeline

---

## 1. Overview
This repository provides an end-to-end R-based pipeline (`structural_variant_finder.R`) to:
1. Parse long-read structural-variant (SV) calls from PacBio data (tumor vs. normal).  
2. Annotate each SV with gene, functional impact, and clinical evidence (CIViC, OncoKB, COSMIC).  
3. Compute a combined pathogenicity score (ACMG/AMP–style thresholds) and classify each SV as Pathogenic, Likely Pathogenic, VUS, or Benign.  
4. Produce a final Excel report with clinical-significance flags and summary plots (via `ggplot2`).

This approach was designed for cancer‐susceptibility studies—e.g., comparing tumor tissue against matched normal to identify clinically actionable SVs.

---

## 2. Biological & Analytical Context
- **Why PacBio?**  
  Long‐read (PacBio) sequencing can capture insertions, deletions, and complex rearrangements (translocations, inversions, tandem duplications) that short‐read methods often miss.  
- **Clinical Significance:**  
  By integrating CIViC, OncoKB, and COSMIC data, this pipeline assigns each SV a “Clinical Evidence Score.” Combined with functional-impact tags (Oncogene, Tumor Suppressor, Cancer Gene, or Uncertain), the pipeline applies simple ACMG/AMP–style thresholds to classify variants for further curation or reporting.  
- **Intended Use Case:**  
  Applied to tumor/normal pairs in cancer samples, with the goal of flagging novel or known SVs driving oncogenesis. Ideal for research labs or clinical‐research settings that need a reproducible, R-based workflow (no Conda environment required).

---

## 3. Prerequisites

1. **Operating System**  
   - Linux or macOS (tested on Ubuntu 20.04 and macOS Big Sur)

2. **Software**  
   - **R** (≥ 4.0)  
   - **R packages** (install with `install.packages()` or `BiocManager::install()`):  
     ```
     readr
     dplyr
     stringr
     openxlsx
     clusterProfiler
     org.Hs.eg.db
     httr
     jsonlite
     ggplot2
     ```
   - **Java** (≥ 8) — required by some annotation tools (if you choose to extend the pipeline).  
   - **(Optional) SLURM / HPC** — for large sample sets, you can submit `structural_variant_finder.R` jobs in parallel on a cluster.

3. **Input Files**  
   - Tumor SV calls exported as a tab‐delimited file (`.tsv` or `.txt`) with columns:  
     ```
     Chr, Start, End, SV_Type, Gene_name, [other fields]
     ```
   - Matched normal SV calls in the same format.  
   - Local COSMIC TSV (to annotate known cancer SVs). Place it under `data/` (see Folder Structure).  
   - **API keys** for querying CIViC and OncoKB (see Usage).

4. **Quick Start**  
For full methods, scoring criteria, advanced usage notes (e.g. API‐skipping modes, SLURM examples, change logs), and example outputs, see [docs/README.pdf](docs/README.pdf).
