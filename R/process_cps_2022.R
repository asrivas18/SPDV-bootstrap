# =============================================================================
# R script: Read and preprocess CPS household income dataset (cps_income.sas7bdat)
# Equivalent to your SAS program
# For manuscript: Scalable Studentized Bootstrap Variance Inference...
# Output:
#   - data/cps_income_clean.rds / .csv    (full cleaned data)
#   - data/cps_income_subset_60k.rds / .csv   (optional 60,000-row subsample)
# =============================================================================

# Load required packages
library(haven)      # for read_sas()
library(dplyr)      # for data manipulation
library(data.table) # for fast CSV writing (fwrite)

# -----------------------------------------------------------------------------
# 1. File paths – adjust these to your local locations
# -----------------------------------------------------------------------------
input_file <- "cps_income.sas7bdat"     # input SAS dataset
output_dir <- "data_folder"          # output directory

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
 
 # -----------------------------------------------------------------------------
 # 2. Read the dataset
 # -----------------------------------------------------------------------------
 cat("Reading CPS household file...\n")
Reading CPS household file...
 cps_raw <- read_sas(input_file)
 str(cps_raw)

# ----------------------------------------------------------------------------- 
# 3. Basic cleaning – keep only relevant variables and positive income
# -----------------------------------------------------------------------------
cps_clean <- cps_raw %>%
     transmute(
         Household_ID    = as.character(SERIAL),          # or H_SEQ if different
         Total_HH_Income = as.numeric(HHINCOME),           # Total household income
         HH_Weight       = as.numeric(HWTFINL),           # Household weight
         Survey_Year     = as.integer(YEAR)             # Year (if present)
     ) %>%
     filter(
         !is.na(Total_HH_Income),
         Total_HH_Income > 0
     ) %>%
     # Optional: keep only 2022 if multiple years exist
     filter(is.na(Survey_Year) | Survey_Year == 2022)
 
cat("Rows after cleaning:", nrow(cps_clean), "\n")



# -----------------------------------------------------------------------------
# 4. Summary statistics – for verification
# -----------------------------------------------------------------------------
cat("\nSummary of Total Household Income (after cleaning):\n")
print(summary(cps_clean$Total_HH_Income))

cat("\nNumber of observations:", nrow(cps_clean), "\n")

# -----------------------------------------------------------------------------
# 5. Save cleaned dataset (full version)
# -----------------------------------------------------------------------------
saveRDS(cps_clean, file.path(output_dir, "cps_income_clean.rds"))
fwrite(cps_clean, file.path(output_dir, "cps_income_clean.csv"))

cat("Saved:\n")
cat("  - cps_income_clean.rds\n")
cat("  - cps_income_clean.csv\n")

# -----------------------------------------------------------------------------
# 6. Optional: Create reproducible 60,000-row subsample (manuscript size)
# -----------------------------------------------------------------------------
set.seed(2022)  # Reproducible seed
target_n <- 60000

if (nrow(cps_clean) > target_n) {
  cps_subset <- cps_clean %>%
    slice_sample(n = target_n)
  
  cat("Subsampled to", nrow(cps_subset), "households.\n")
  
  # Save subsample
  saveRDS(cps_subset, file.path(output_dir, "cps_income_subset_60k.rds"))
  fwrite(cps_subset, file.path(output_dir, "cps_income_subset_60k.csv"))
  
  cat("Saved subsample:\n")
  cat("  - cps_income_subset_60k.rds\n")
  cat("  - cps_income_subset_60k.csv\n")
} else {
  cat("Dataset already smaller than 60,000 – no subsampling performed.\n")
}

# -----------------------------------------------------------------------------
# 7. Quick final check
# -----------------------------------------------------------------------------
cat("\nFinal verification – subsample summary:\n")
print(summary(cps_subset$Total_HH_Income))