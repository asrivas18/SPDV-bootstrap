# Step 1: Download full 2022 BRFSS public-use data (LLCP 2022)
# URL: https://www.cdc.gov/brfss/annual_data/annual_2022.html
# Direct data file: https://www.cdc.gov/brfss/annual_data/2022/LLCP2022XPT.zip (SAS XPT format)
# Unzip and convert to CSV or read directly with foreign/Haven

library(haven)   # for read_xpt
library(dplyr)

# 1. Define the path clearly
file_path <- "C:/Users/ssrivas/Box/Important_Folder/Research_Folders/Working_folder_articles/papers/New_papers/JSCS/Repository/Data/R-programs"

# 2. Read and filter by Louisiana FIPS code (22)
# Note: Variables starting with underscores must be in backticks
brfss_full <- read_xpt("LLCP2022.XPT") %>%
  filter(`_STATE` == 22)  

# 3. Clean BMI data
# Standard BRFSS cleaning: _BMI5 is the calculated variable. 
# It is scaled by 100 (e.g., 2400 = 24.0)
brfss_bmi <- brfss_full %>%
  filter(!is.na(`_BMI5`)) %>%           
  transmute(
    SEQNO  = as.character(SEQNO),
    BMI    = as.numeric(`_BMI5`) / 100,  # Convert to standard BMI scale
    STATE  = factor(`_STATE`, labels = "Louisiana"),
    # IDATE is usually MMDDYYYY string
    IDATE  = as.Date(IDATE, format = "%m%d%Y") 
  ) %>%
  filter(complete.cases(.))              

# 4. Save and Verify
saveRDS(brfss_bmi, file = "brfss_bmi_subset.rds")

cat("Successfully cleaned BRFSS subset.\n")
cat("Final Count:", nrow(brfss_bmi), "rows\n")
print(summary(brfss_bmi$BMI))
