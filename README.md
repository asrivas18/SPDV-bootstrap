# SPDV-Bootstrap: Scalable studentized bootstrap variance inference under skewness and heavy tails

This repository contains code, data, and supplementary materials to reproduce the main results from the manuscript:

**Title:** Scalable studentized bootstrap variance inference under skewness and heavy tails  
**Authors:** Sudesh K. Srivastav, Apurv Srivastav  
**Target journal:** Journal of Applied Statistics (JAS)  
**GitHub:** [https://github.com/asrivas18/SPDV-bootstrap](https://github.com/asrivas18/SPDV-bootstrap)  
**Zenodo archive:** DOI to be added upon acceptance

## Overview

This project develops a scalable nonparametric framework for studentized (bootstrap-t) inference for the population variance. It exploits the classical algebraic identity between the sample pairwise difference variance (SPDV) and the unbiased sample variance (s^2), allowing exact (O(n)) variance computation per bootstrap resample instead of direct (O(n^2)) pairwise evaluation. As a result, the total computational cost becomes (O(Bn)), making high-replication bootstrap inference feasible for large sample sizes and modern parallel hardware.

### Key results

* **Theory:** Nonparametric bootstrap consistency is established under a finite fourth-moment condition using the equivalent non-degenerate U-statistic representation of the sample variance.
* **Empirical performance:** In simulation settings with skewness and heavy tails, studentized intervals substantially outperform classical chi-square intervals and improve markedly over normal-approximation and percentile bootstrap intervals.
* **Scalability:** The identity-based implementation delivers several orders of magnitude speedup relative to direct pairwise evaluation and enables high-replication studentized inference at modern data scales.

## Repository contents

```text
.
├── README.md
├── LICENSE
├── sessionInfo.txt
├── ms\_article.tex
├── ms\_article.pdf
├── supplement.tex
├── supplement.pdf
├── R/
│   ├── spdv\_functions.R
│   ├── run\_simulation\_study.R
│   ├── run\_real\_data\_examples.R
│   ├── process\_cps\_2022.R
│   └── process\_brfss\_2022.R
├── data/
│   ├── README-data.md
│   ├── brfss\_bmi\_subset.rds
│   └── cps\_income\_subset\_60k.rds
├── results/
│   ├── coverage\_results.rds
│   ├── runtime\_comparison.csv
│   └── real\_data\_results.rds
└── figures/
    ├── coverage\_plot.pdf
    └── sensitivity\_B\_plot.pdf
```

### File guide

* `ms\_article.tex` / `ms\_article.pdf`: Main manuscript.
* `supplement.tex` / `supplement.pdf`: Online supplement with additional computational benchmarks.
* `R/spdv\_functions.R`: Core computational functions for SPDV-based variance estimation and studentized bootstrap inference.
* `R/run\_simulation\_study.R`: Reproduces the simulation study, including empirical coverage and sensitivity to the number of bootstrap replications.
* `R/run\_real\_data\_examples.R`: Reproduces the BMI and CPS real-data examples.
* `R/process\_cps\_2022.R`, `R/process\_brfss\_2022.R`: Data preprocessing scripts.
* `data/`: Preprocessed subsets used in the manuscript examples.
* `results/`: Saved output objects and benchmark summaries.
* `figures/`: Plots used in the manuscript.

## System requirements

* R >= 4.4.0 (tested on R 4.4.1)
* Recommended packages: `ggplot2`, `dplyr`, `tidyr`, `foreach`, `doParallel`, `doRNG`, `kableExtra`, `haven`, `moments`
* Recommended hardware: multi-core CPU (16 cores used for reported benchmarks), at least 16 GB RAM

## Installation

Install required packages in R:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "foreach", "doParallel",
  "doRNG", "kableExtra", "haven", "moments"
))
```

## How to reproduce results

### 1\. Simulation study

Reproduces the main coverage results and sensitivity-to-(B) analysis.

```r
source("R/run\_simulation\_study.R")
```

Expected outputs:

* Empirical coverage summaries for the competing confidence intervals.
* Sensitivity results for the number of bootstrap replications.
* Saved files such as:

  * `results/coverage\_results.rds`
  * `figures/coverage\_plot.pdf`
  * `figures/sensitivity\_B\_plot.pdf`

### 2\. Real-data examples

Reproduces the BRFSS BMI and CPS household income analyses from the manuscript.

```r
source("R/run\_real\_data\_examples.R")
```

Expected outputs:

* **BMI (Louisiana BRFSS 2022, (n = 5{,}117))**: variance approximately 48.03, with studentized interval approximately \[45.49, 50.90].
* **CPS income (2022 CPS, (n = 60{,}000))**: variance approximately (16.89 \\times 10^9), with studentized interval approximately \[16.04, 17.85] in units of (10^9).

See `data/README-data.md` for details on the included subsets.

### 3\. Preprocess CPS data from scratch (optional)

```r
source("R/process\_cps\_2022.R")
```

This script requires the raw CPS data file in the expected `raw\_data/` location and produces the processed subset used in the manuscript.

### 4\. View session information

```r
sessionInfo()
```

Compare the output to `sessionInfo.txt` in the repository root.

## Simulation design summary

The simulation study considers sample sizes (n \\in {20, 50, 100}) under four distributions with unit variance:

* Normal
* (t\_5)
* (t\_3)
* Lognormal

All distributions are linearly rescaled to unit variance so that comparisons reflect differences in skewness and tail behavior rather than scale. The study compares chi-square, normal-approximation, percentile bootstrap, and studentized bootstrap-(t) intervals.

## Reproducibility notes

* Random-number generation is controlled using `doRNG` with L'Ecuyer-CMRG streams.
* Parallel execution uses `foreach` and `doParallel`.
* Results are deterministic given the same random seed, R version, and package versions.
* In a small fraction of heavy-tailed simulation replications, near-zero resample-specific standard errors can produce extreme studentized pivots; these cases are retained in the reported results.

## Data notes

* **BRFSS BMI subset:** Louisiana adults with complete BMI records after preprocessing.
* **CPS income subset:** Positive household income values from the 2022 CPS extract used for the manuscript illustration.
* Preprocessing scripts are included so that readers can reconstruct the analysis datasets from the original public sources when permitted.

## Citation

If this repository is useful in your work, please cite the associated manuscript:

```bibtex
@article{srivastav\_spdv\_bootstrap,
  title   = {Scalable studentized bootstrap variance inference under skewness and heavy tails},
  author  = {Srivastav, Sudesh K. and Srivastav, Apurv},
  journal = {Journal of Applied Statistics},
  year    = {submitted},
  note    = {Code repository: https://github.com/asrivas18/SPDV-bootstrap}
}
```

## License

This repository is released under the MIT License. See `LICENSE` for details.

## Contact

**Sudesh K. Srivastav**  
Department of Biostatistics and Data Science  
Tulane University  
Email: [ssrivas@tulane.edu](mailto:ssrivas@tulane.edu)

