# Title: "01_create_ds_domain.R"
# Author: "Sammie Chum"
# Date: "2026-03-31"

# Install and load required packages -------------------------------------------
# install.packages("sdtm.oak")
# install.packages("dplyr")
# install.packages("pharmaverseraw")
# install.packages("pharmaversesdtm")

library(sdtm.oak)
library(pharmaverseraw)
library(dplyr)
library(pharmaversesdtm)

# Load Data --------------------------------------------------------------------
dm <- pharmaversesdtm::dm

# Load DS raw input data
ds_raw <- pharmaverseraw::ds_raw

# Load study controlled terminology
study_ct <- read.csv("data/sdtm_ct.csv")


# Create oak_id_vars -----------------------------------------------------------
# Create oak id variables
ds_raw <- ds_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ds_raw"
  )

# Combine DSDTCOL (date) + DSTMCOL
ds_raw <- ds_raw %>%
  dplyr::mutate(
    DSDTCOL_DT = dplyr::if_else(
      !is.na(DSTMCOL) & DSTMCOL != "",
      paste0(DSDTCOL, " ", DSTMCOL),
      DSDTCOL
    )
  )

# Change to uppercase to match study controlled terminology
ds_raw <- ds_raw %>%
  mutate(IT.DSDECOD = toupper(IT.DSDECOD))

# Map Topic Variable -----------------------------------------------------------
ds <-
  # Derive topic variable
  # Map DSTERM using assign_no_ct, raw_var=IT.DSTERM, tgt_var=DSTERM
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSTERM",
    tgt_var = "DSTERM",
    id_vars = oak_id_vars()
  )

# Map Additional Variables  ----------------------------------------------------
ds <- ds %>%
  # Map DSDECOD using assign_ct, raw_var=IT.DSTERM, tgt_var=DSDECOD
  assign_ct(
    raw_dat = ds_raw,
    raw_var = "IT.DSDECOD",
    tgt_var = "DSDECOD",
    ct_spec = study_ct,
    ct_clst = "C66727", # From study controlled terminology
    id_vars = oak_id_vars()
  ) %>%
  # Map DSCAT: derived from DSDECOD and OTHERSP 
  # For DSDECOD, Randomized and Completed flags indicate Protocol Milestone,
  # if there is event in OTHERSP it also indicates Protocol Milestone (final visits)
  # Others (not null) in IT.DSDECOD set to disposition
  assign_no_ct(
    raw_dat = ds_raw %>%
      dplyr::mutate(
        DSCAT_RAW = dplyr::case_when(
          IT.DSDECOD %in% c("Randomized", "Completed")     ~ "PROTOCOL MILESTONE",
          !is.na(OTHERSP)                                   ~ "PROTOCOL MILESTONE",
          !is.na(IT.DSDECOD)                               ~ "DISPOSITION EVENT",
          TRUE                                              ~ NA_character_
        )
      ),
    raw_var = "DSCAT_RAW",
    tgt_var = "DSCAT",
    id_vars = oak_id_vars()
  ) %>%
  # Map VISIT using assign_no_ct, raw_var=INSTANCE, tgt_var=VISIT
  assign_no_ct(
    raw_dat = ds_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    id_vars = oak_id_vars()
  ) %>%
  # Map VISITNUM: derived from INSTANCE and assigned values accoording to 
  # INSTANCE name logically
  # Following SDTMIG v3.4 unplanned visits had VISITNUM set to 99
  assign_no_ct(
    raw_dat = ds_raw %>%
      dplyr::mutate(
        VISITNUM_RAW = dplyr::case_when(
          INSTANCE == "Screening 1"        ~   1,
          INSTANCE == "Baseline"           ~   2,
          INSTANCE == "Week 2"             ~   3,
          INSTANCE == "Week 4"             ~   4,
          INSTANCE == "Week 6"             ~   5,
          INSTANCE == "Week 8"             ~   6,
          INSTANCE == "Week 12"            ~   7,
          INSTANCE == "Week 16"            ~   8,
          INSTANCE == "Week 20"            ~   9,
          INSTANCE == "Week 24"            ~  10,
          INSTANCE == "Week 26"            ~  11,
          INSTANCE == "Retrieval"          ~  12,
          INSTANCE == "Ambul Ecg Removal"  ~  13,
          INSTANCE == "Unscheduled 1.1"    ~ 99,
          INSTANCE == "Unscheduled 4.1"    ~ 99,
          INSTANCE == "Unscheduled 5.1"    ~ 99,
          INSTANCE == "Unscheduled 6.1"    ~ 99,
          INSTANCE == "Unscheduled 8.2"    ~ 99,
          INSTANCE == "Unscheduled 13.1"   ~ 99,
          TRUE                             ~  NA_real_
        ) %>% as.character()
      ),
    raw_var = "VISITNUM_RAW",
    tgt_var = "VISITNUM",
    id_vars = oak_id_vars()
  ) %>%
  # Map DSDTC using assign_datetime, raw_var=IT.DSDTCOL
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "DSDTCOL",
    tgt_var = "DSDTC",
    raw_fmt = c("m-d-y")
  ) %>%
  # Map DSSTDTC using assign_datetime, raw_var=IT.DSSTDAT
  assign_datetime(
    raw_dat = ds_raw,
    raw_var = "IT.DSSTDAT",
    tgt_var = "DSSTDTC",
    raw_fmt = c("m-d-y")
  )


# Create SDTM Derived Variables  -----------------------------------------------
ds <- ds %>%
  # Add in STUDYID, DOMAIN, and USUBJID
  dplyr::mutate(
    STUDYID = ds_raw$STUDY,
    DOMAIN = "DS",
    USUBJID = paste0("01-", ds_raw$PATNUM)
  ) %>%
  # Add DSSEQ
  derive_seq(
    tgt_var = "DSSEQ",
    rec_vars = c("USUBJID", "DSTERM")
  ) %>%
  # Add Study Day
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "DSSTDTC",
    refdt = "RFXSTDTC",
    study_day_var = "DSSTDY"
  ) %>%
  # Include all specified variables
  select(
    "STUDYID", "DOMAIN", "USUBJID", "DSSEQ", "DSTERM", "DSDECOD", "DSCAT",
    "VISITNUM", "VISIT", "DSDTC", "DSSTDTC", "DSSTDY"
  )

# Save file as csv
write.csv(ds, "data/sdtm_dataset.csv")
