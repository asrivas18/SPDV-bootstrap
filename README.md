# SPDV-Bootstrap: Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation
This repository contains the code, data, and materials to reproduce the results in the manuscript:

**Title**: Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation
**Authors**: Sudesh K. Srivastav, Apurv Srivastav  
**Journal**: Journal of Statistical Computation and Simulation (JSCS)  
**Submitted**: [Month Year]  
**GitHub**: https://github.com/asrivas18/SPDV-bootstrap 
**Zenodo Archive**: DOI forthcoming upon acceptance (full replication kit with large simulation summaries)

## Overview

This work introduces a scalable framework for nonparametric studentized (bootstrap-t) inference on the population variance using the algebraic identity between the Sample Pairwise Difference Variance (SPDV) and the classical unbiased sample variance s². The identity enables O(n) per-resample computation instead of O(n²), making high-replication (B ≥ 50,000) bootstrap feasible for n ≈ 10⁶ on standard hardware.

Key results:
- Theoretical: Bootstrap consistency under finite fourth moments.
- Empirical: Studentized intervals substantially outperform classical chi-square in non-normal settings (Table 3).
- Computational: Up to 5 orders of magnitude speedup (Tables 1 & 3).

## Repository Contents

- `simulation.R` / `coverage_study.R`: Main script to reproduce Table 3 (coverage probabilities) and Figure 2.
- `bootstrap_functions.R` / `spdv_boot.R`: Core functions for SPDV, studentized bootstrap, influence-function SE.
- `coverage_results.rds` / `coverage_results.csv`: Raw Monte Carlo results (5,000 × 4 × 3 × 5,000 = 300 million bootstrap runs).
- `Figure2_Coverage_color.pdf`: Color version of coverage plot.
- `Figure2_Coverage_CSDA.pdf`: Grayscale/Elsevier-style version for submission.
- `report.Rmd`: R Markdown source for full simulation report and reproducibility.
- `sessionInfo.txt`: R session information (R 4.4.1, packages, etc.).
- `data/bmi_subset.rds` or `.csv`: Preprocessed BRFSS BMI subset for real-data application.
- `LICENSE`: MIT License (permissive for academic reuse).

## System Requirements

- R ≥ 4.4.0
- Packages: ggplot2, dplyr, tidyr, foreach, doParallel, doRNG, kableExtra
- Hardware: 16-core CPU recommended for parallel simulation (43 min wall-clock on AMD EPYC 7543P, 128GB RAM).
- RAM: ≥ 16 GB (for large B and n)

## Installation

```bash
# Install required packages (run in R console)
install.packages(c("ggplot2", "dplyr", "tidyr", "foreach", "doParallel", "doRNG", "kableExtra"))

#########################################################

How to Reproduce Results

Reproduce Table 3 (Coverage) and Figure 2:R# Set working directory to repo root
setwd("path/to/SPDV-bootstrap")

# Run the full simulation (M=5000, B=5000 – ~43 minutes parallel)
source("simulation.R")

# Or quick test mode (M=500, B=500 – ~2 minutes)
source("simulation_quick.R")Output:
Console: Table 3 printed
Files: coverage_results.rds, coverage_results.csv
Plots: Figure2_Coverage_color.pdf, Figure2_Coverage_CSDA.pdf

Reproduce BMI Application (Section 6):Rsource("applications/bmi_application.R")
Loads bmi_subset.rds
Runs B=50,000 studentized bootstrap
Outputs CI: [21.88, 23.59] for variance

View session info:RsessionInfo()Compare to sessionInfo.txt

Reproducibility Notes

Random seed: 12345 (via doRNG L'Ecuyer-CMRG)
Parallel backend: doParallel + foreach
All results are deterministic given the seed and R version.
Skipped cases in studentized bootstrap (unstable SE) are logged in simulation output.

Citation
If you use this code, please cite the manuscript:
Srivastav, S. K., & Srivastav, A. (submitted). Scalable Studentized Bootstrap Inference for Variance via Pairwise Difference Representation. Computational Statistics & Data Analysis.
License
MIT License – see LICENSE file.
Contact
Sudesh K. Srivastav
ssrivas@tulane.edu

Department of Biostatistics and Data Science, Tulane University

