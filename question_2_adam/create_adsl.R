# Title: "create_adsl.R"
# Author: "Sammie Chum"
# Date: "2026-03-31"

# Install and load required packages -------------------------------------------
library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)

# Load Data --------------------------------------------------------------------
dm <- pharmaversesdtm::dm
vs <- pharmaversesdtm::vs
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae

dm <- convert_blanks_to_na(dm)
vs <- convert_blanks_to_na(vs)
ex <- convert_blanks_to_na(ex)
ds <- convert_blanks_to_na(ds)
ae <- convert_blanks_to_na(ae)

# Use DM domain as the basis for ADSL
adsl <- dm %>%
  select(-DOMAIN)


# TRTSDTM & TRTSTMF Derivation -------------------------------------------------
# Convert data from ex to datetime 
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc                = EXSTDTC,
    new_vars_prefix    = "EXST",
    time_imputation    = "00:00:00", # Impute missing time with 00:00:00
    # If only seconds are missing do not populate the imputation flag
    ignore_seconds_flag = TRUE 
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "00:00:00"
  )

# Derive treatment start and end time
# Only include observations where the patient received a valid dose 
# and datepart of Start Date/Time of Treatment [EX.EXSTDTC] is complete. 
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    # Filter for valid dose 
    filter_add  = (EXDOSE > 0 |
                     (EXDOSE == 0 &
                        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars    = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order       = exprs(EXSTDTM, EXSEQ),
    # Select patient's first exposure
    mode        = "first",                
    by_vars     = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )

adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))


# AGEGR9 & AGEGR9N Derivation -------------------------------------------------
# Group age into the following categories: <18, 18-50, >50
format_agegr9 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 50) ~ "18-50",
    var_input > 50 ~ ">50",
    TRUE ~ "Missing"
  )
}

# Provide numeric groupings (1, 2, 3) for the ages <18, 18-50, >50
format_agegr9n <- function(var_input) {
  case_when(
    var_input == "<18" ~ 1,
    var_input == "18-50" ~ 2,
    var_input == ">50" ~ 3,
    var_input == "Missing" ~ NA_real_
  )
}

# Apply first age grouping: <18, 18-50, >50
adsl <- adsl %>%
  mutate(
    AGEGR9 = format_agegr9(AGE)
  ) 

# Then apply numeric grouping: 1, 2, 3
adsl <- adsl %>%
  mutate(
    AGEGR9N = format_agegr9n(AGEGR9)
  )


# ITTFL Derivation -------------------------------------------------------------
# Set to "Y" if [DM.ARM] not equal to missing Else set to "N"
adsl <- adsl %>%
  mutate (
    ITTFL = if_else(!is.na(ARM), "Y", "N")
  )


# LSTAVLDT Derivation ----------------------------------------------------------
# Set to the last date patient has documented clinical data to show him/her
# alive, converted to numeric date, using a series of 4 possible dates
# Set to the max of these dates
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      # (1) Last complete date of vital assessment
      event(
        dataset_name = "vs",
        order        = exprs(VSDTC, VSSEQ),
        condition    = !is.na(VSDTC) &
          nchar(VSDTC) >= 10 &        # VS.VSDTC not missing
          !(is.na(VSSTRESN) & is.na(VSSTRESC)),     # Not both missing
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(VSDTC, highest_imputation = "M"),
          seq      = VSSEQ
        )
      ),
      # (2) Last complete onset date of AE 
      event(
        dataset_name = "ae",
        order        = exprs(AESTDTC, AESEQ),
        condition    = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq      = AESEQ
        )
      ),
      # (3) Last complete disposition date 
      event(
        dataset_name = "ds",
        order        = exprs(DSSTDTC, DSSEQ),
        condition    = !is.na(DSSTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(DSSTDTC, highest_imputation = "M"),
          seq      = DSSEQ
        )
      ),
      # (4) Last date of treatment administration
      event(
        dataset_name  = "adsl",
        condition     = !is.na(TRTEDTM),
        set_values_to = exprs(
          LSTALVDT = date(TRTEDTM),                                
          seq      = 0
        )
      )
    ),
    source_datasets = list(vs = vs, ae = ae, ds = ds, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order   = exprs(LSTALVDT, seq, event_nr),
    mode    = "last",                                            
    new_vars = exprs(LSTALVDT)
  )

write.csv(adsl, "data/adsl.csv")
