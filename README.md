# SPDV-Bootstrap: Scalable Studentized Bootstrap Variance Inference  
via Linear-Time Pairwise Variance Representation

This repository contains all code, data, and materials to reproduce the results in the manuscript:

**Title**: Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation  
**Authors**: Sudesh K. Srivastav, Apurv Srivastav  
**Journal**: Journal of Statistical Computation and Simulation (JSCS)  
**Status**: Submitted  
**GitHub**: https://github.com/asrivas18/SPDV-bootstrap  
**Zenodo Archive**: DOI forthcoming upon acceptance

## Overview

This work presents a scalable nonparametric framework for studentized (bootstrap-t) inference on the population variance. It exploits the classical algebraic identity between the Sample Pairwise Difference Variance (SPDV) and the unbiased sample variance $s^2$, enabling exact $O(n)$ computation per bootstrap resample instead of $O(n^2)$. This reduces total complexity to $O(Bn)$, making high-replication ($B \ge 50{,}000$) inference feasible for $n \approx 10^6$ on standard parallel hardware.

Key results:
- Theoretical: Bootstrap consistency under finite fourth moments (Theorem 1).
- Empirical: Studentized intervals substantially outperform classical chi-square in non-normal settings (Table 2, Figure 1).
- Computational: Up to 5–6 orders of magnitude speedup vs. naive pairwise implementation (Tables 1 & 2).

## Repository Contents
.
├── README.md                          # This file
├── LICENSE                            # MIT License
├── sessionInfo.txt                    # R session & package versions
├── SPDV_manuscript.tex                # Full LaTeX source
├── SPDV_manuscript.pdf                # Compiled PDF
├── R/
│   ├── spdv_functions.R               # Core functions: spdv_fast(), studentized_bootstrap()
│   ├── run_simulation_study.R         # Reproduce Table 3, Figure 2, Table 5 (sensitivity to B)
│   ├── run_real_data_examples.R       # BMI & CPS real-data examples (Section 6)
│   ├── process_cps_2022.R             # Preprocess CPS income subset
│   └── process_brfss_2022.R           # Preprocess BRFSS BMI subset
├── data/
│   ├── README-data.md
│   ├── brfss_bmi_subset.rds           # Louisiana BRFSS BMI (~5,432 obs)
│   └── cps_income_subset_60k.rds      # CPS household income subset (~60,000 obs)
├── results/
│   ├── coverage_results.rds           # Raw Monte Carlo output
│   └── runtime_comparison.csv         # Benchmark timings
└── figures/
├── coverage_plot_JSCS.pdf         # Main coverage figure (Figure 1)
└── sensitivity_B_plot.pdf         # Coverage vs. B (Figure 2)
text## System Requirements

- R ≥ 4.4.0 (tested on 4.4.1)
- Key packages: ggplot2, dplyr, tidyr, foreach, doParallel, doRNG, kableExtra, haven, moments
- Hardware: 16-core CPU recommended (full simulation ~43 min parallel); ≥16 GB RAM

## Installation

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "foreach", "doParallel", "doRNG", "kableExtra", "haven", "moments"))
How to Reproduce Results
1. Simulation Study (Tables 2, 3; Figures 1, 2)
Rsource("R/run_simulation_study.R")
Outputs:

Console: Table 2 (coverage), Table 3 (sensitivity to B)
Files: results/coverage_results.rds, figures/coverage_plot_JSCS.pdf, figures/sensitivity_B_plot.pdf

2. Real-Data Examples (Section 6)
Rsource("R/run_real_data_examples.R")
Outputs:

BMI (Louisiana BRFSS 2022): variance ≈22.71, 95% CI [21.88, 23.59], root pairwise [6.61, 6.88]
CPS income (2022 ASEC): variance ≈2.1×10⁹ USD², 95% studentized CI (heavy-tailed)

See data/README-data.md for subset details.
3. Preprocess CPS data from scratch (optional)
Rsource("R/process_cps_2022.R")
Requires raw hhpub22.sas7bdat in raw_data/. Produces data/cps_income_subset_60k.rds/csv.
4. View session information
RsessionInfo()
# Compare to sessionInfo.txt in repo root
Reproducibility Notes

Random seeds: Fixed via doRNG (L'Ecuyer-CMRG) in all scripts.
Parallel: doParallel + foreach (auto-detects cores).
All results deterministic given seed, R version, and packages.
Unstable studentized pivots (<5% in small n heavy tails) retained — contributes to slight conservativeness.

Citation
If you use this code or data, please cite:
Srivastav, S. K., & Srivastav, A. (submitted). Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation. Journal of Statistical Computation and Simulation.
bibtex@article{srivastav2025spdv,
  title={Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation},
  author={Srivastav, Sudesh K. and Srivastav, Apurv},
  journal={Journal of Statistical Computation and Simulation},
  year={submitted},
  note={Code available at https://github.com/asrivas18/SPDV-bootstrap}
}
License
MIT License — see LICENSE file.
Contact
Sudesh K. Srivastav
ssrivas@tulane.edu
Department of Biostatistics and Data Science, Tulane University