# =============================================================================
# real_data_bmi.R — Louisiana BRFSS 2022 BMI Application
# Matches manuscript Section 6.1 (n ≈ 5,432, B=50,000 studentized bootstrap)
# =============================================================================

# real_data_bmi.R (BMI Application — matches article exactly)

library(haven)
library(tidyverse)

source("bootstrap_functions.R")

# Download 2022 BRFSS (public CDC XPT via ZIP)
#cat("Downloading 2022 BRFSS data...\n")
#zip_url <- "https://www.cdc.gov/brfss/annual_data/2022/files/LLCP2022XPT.zip"
#temp_zip <- tempfile(fileext = ".zip")
#download.file(zip_url, temp_zip, mode = "wb")
#xpt_file <- unzip(temp_zip, exdir = tempdir())

# Read XPT
#brfss <- read_xpt(xpt_file[grepl("\\.XPT$", xpt_file, ignore.case = TRUE)])

# Load national file
brfss <- read_xpt("LLCP2022.XPT")

# Filter for Louisiana (22) and create binary smoking status
brfss_la_2022 <- brfss %>%
  filter(`_STATE` == 22) %>%
  mutate(smoke_status = ifelse(`_RFSMOK3` == 2, 1, 0)) %>%
  select(smoke_status) %>%
  na.omit()

# Process: Louisiana (_STATE == 22), BMI from _BMI5 (in hundredths)
bmi_data <- brfss %>%
  filter(`_STATE` == 22) %>%
  mutate(BMI = `_BMI5` / 100) %>%   # _BMI5 is in hundredths
  filter(!is.na(BMI), BMI > 0) %>%
  pull(BMI)

# Results
cat("Louisiana complete BMI cases: n =", length(bmi_data), "\n")
result <- studentized_bootstrap(bmi_data, B = 50000)

cat("Sample variance s² ≈", round(result$variance, 2), "\n")
cat("Root pairwise difference ≈", round(result$root_pairwise, 2), "\n")
cat("95% CI for variance:", round(result$ci_variance, 2), "\n")
cat("95% CI for root pairwise:", round(result$ci_root, 2), "\n")

# Cleanup
#unlink(temp_zip)
#unlink(xpt_file)