# Step 1: Download 2022 CPS ASEC public-use microdata
# URL: https://www.census.gov/data/datasets/time-series/demo/cps/cps-asec.html
# File: https://www2.census.gov/programs-surveys/cps/datasets/2022/march/asecpub22sas.zip (SAS format)
# Or use IPUMS CPS extract (easier): https://cps.ipums.org/cps/

library(haven)
library(dplyr)

cps_full <- read_sas("path/to/asecpub22sas/hhld.sas7bdat")  # household file

cps_income <- cps_full %>%
  filter(YEAR == 2022) %>%
  transmute(
    SERIAL = as.character(SERIAL),
    HHINC  = as.numeric(HTOTVAL),     # Total household income (adjust var name if needed)
    YEAR   = as.integer(YEAR),
    WEIGHT = as.numeric(HWHHWGT)      # Household weight (optional)
  ) %>%
  filter(HHINC >= 0 & !is.na(HHINC))  # Remove negative/missing income

# Optional: take a 60,000-row random subset if full file too large
set.seed(2022)
cps_income <- cps_income %>% slice_sample(n = 60000)

saveRDS(cps_income, file = "data/cps_income_subset.rds")
cat("Saved:", nrow(cps_income), "rows\n")