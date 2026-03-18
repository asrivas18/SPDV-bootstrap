# Real Data Used in the SPDV-Bootstrap Manuscript

This folder contains **preprocessed, public-use subsets** of the two real datasets analyzed in Section 6 of the manuscript:

> *Scalable Studentized Bootstrap Variance Inference via Linear-Time Pairwise Variance Representation*

The subsets are small, anonymized, and analysis-ready. They are provided to enable full reproducibility of the real-data results without requiring users to download and process large raw files themselves.

## Files

### 1. `brfss_bmi_subset.rds`
- **Source**: 2022 Louisiana Behavioral Risk Factor Surveillance System (BRFSS) public-use data  
- **Provider**: Louisiana Department of Health (cited as LA_BRFSS_2022 in the manuscript)  
- **Contents**:  
  - Single column: `Total_HH_Income` (numeric, body mass index in kg/m²)  
  - Or a simple numeric vector of BMI values  
- **Processing steps**:  
  - Restricted to Louisiana respondents (`_STATE == 22`)  
  - Age ≥ 18  
  - Non-missing BMI values only (from `_BMI5 / 100`)  
  - No personal identifiers retained  
- **Size**: ≈ 5,117 observations (exact count as reported in the paper)  
- **Used in**: Section 6.1 (Public Health: BMI Heterogeneity)

### 2. `cps_income_subset_60k.rds` (and `.csv`)
- **Source**: 2022 Current Population Survey (CPS) Annual Social and Economic Supplement (ASEC) public-use microdata  
- **Provider**: U.S. Census Bureau, accessed via IPUMS CPS (cited as Flood et al., 2023)  
- **Contents**:  
  - Single column: `Total_HH_Income` (numeric, total household income in USD)  
- **Processing steps**:  
  - Extracted total household income (typically `HTOTVAL`)  
  - Filtered to positive, non-missing values  
  - Random subsample to exactly 60,000 observations (seed = 2022)  
  - No restricted variables or personal identifiers retained  
- **Size**: 60,000 observations (exact count as reported in the paper)  
- **Used in**: Section 6.2 (Economic Inequality: Household Income)

## How to Recreate These Subsets from Raw Public Data

Due to size and licensing constraints, only cleaned subsets are stored here. To recreate from original public sources:

1. **BRFSS BMI subset**  
   - Download the 2022 BRFSS public-use dataset from the Louisiana Department of Health or CDC website.  
   - Filter to Louisiana (`_STATE == 22`), age ≥ 18, non-missing BMI (`_BMI5 / 100`).  
   - Save as:  
     ```r
     brfss_bmi_subset <- data.frame(Total_HH_Income = bmi_vector)
     saveRDS(brfss_bmi_subset, "data/brfss_bmi_subset.rds")