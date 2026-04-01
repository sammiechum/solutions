# Title: "01_create_ae_summary_table.R"
# Author: "Sammie Chum"
# Date: "2026-03-31"

  
# Load libraries & data --------------------------------------------------------
# install.packages("pharmaverseadam")
# install.packages("gtsummary")
# install.packages("webshot2")

library(pharmaverseadam)
library(gtsummary)
library(dplyr)

adae <- pharmaverseadam::adae 
adsl <- pharmaverseadam::adsl


# Pre-processing ---------------------------------------------------------------
# Filter for Treatment-emergent AE records
teae <- adae %>%
  filter(TRTEMFL == "Y")


# Build Table ------------------------------------------------------------------
tbl <- teae %>%
  tbl_hierarchical(
    variables   = c(AESOC, AETERM),       
    by          = ACTARM,                
    id          = USUBJID,              
    denominator = adsl,                
    overall_row = TRUE,
    label       = list(
      AESOC ~ "Primary System Organ Class",
      AETERM ~ "Reporter Term for the Adverse Event Term"
    )
  ) %>%
  sort_hierarchical() # Sort by descending frequency by default

# Save to PDF
pdf_tbl <- tbl %>%
  as_gt() %>%                          # convert gtsummary object to gt
  gt::gtsave("data/teae_table.pdf")         # save directly as PDF
