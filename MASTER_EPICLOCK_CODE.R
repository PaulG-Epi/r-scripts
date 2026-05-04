


########################################################
#######################################################

# Epiclock PAH manuscript master file

################################################################################

# Files for install


install.packages('sesame')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")

library(sesame)

## Run once after install
sesameDataCache()
#sesame_checkVersion()


install.packages('methylclock')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylclock")

library(methylclock)


install.packages("remotes")
remotes::install_github("danbelsky/DunedinPACE")

library(DunedinPACE)



install.packages(c("readxl", "dplyr", "gtsummary", "officer", "flextable", "patchwork", "lmerTest"))

library(readxl)
library(dplyr)
library(gtsummary)
library(officer)
library(flextable)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(readr)
library(broom)
library(purrr)
library(lme4)
library(lmerTest)
library(forcats)






################################################################################


## Load files for analyses



# - beta_noob                            (DNAm beta matrix)
# - demographics_0228                    (demographic table)
# - skin_PAH_exposure                    (skin PAH table)
# - gear_PAH_exposure                    (gear PAH table)
# - OKF_EPIC2_SampleInfoFile_SampleInfo_ (sample info with Link)
# - OKF_EPIC2_houseman_cbc               (cell proportions)
#



# names(demographics_0228)
# names(skin_PAH_exposure)
# names(gear_PAH_exposure)
# names(OKF_EPIC2_SampleInfoFile_SampleInfo_)





################################################################################

## MERGE DEMOGRAPHICS + SKIN PAH + LINK KEY

names(demographics_0228)
names(skin_PAH_exposure)



merged_demo <- demographics_0228 %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    skin_PAH_exposure %>% mutate(ID = as.character(ID)),
    by = "ID"
  ) %>%
  tidyr::drop_na()

merged_demo2 <- merged_demo %>%
  mutate(
    ID = as.character(ID),
    timepoint = case_when(
      Test == 1 ~ "V01",
      Test == 2 ~ "V02",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
      mutate(ID = as.character(ID)) %>%
      select(ID, VisitID, Link) %>%
      rename(timepoint = VisitID),
    by = c("ID", "timepoint")
  )

OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  count(ID, VisitID) %>%
  filter(n > 1)

OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  filter((ID == 9 & VisitID == "V02") |
           (ID == 29 & VisitID == "V01") |
           (ID == 37 & VisitID == "V01")) %>%
  select(ID, VisitID, Link)

link_key <- OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  mutate(ID = as.character(ID)) %>%
  select(ID, VisitID, Link) %>%
  arrange(ID, VisitID) %>%
  group_by(ID, VisitID) %>%
  slice(1) %>%
  ungroup()

merged_demo2 <- merged_demo %>%
  mutate(
    ID = as.character(ID),
    VisitID = case_when(
      Test == 1 ~ "V01",
      Test == 2 ~ "V02",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(link_key, by = c("ID", "VisitID"))


################################################################################

## CLOCK CALCULATION (beta_noob -> DNAmAge -> meta_use2)



#Check structure
str(beta_noob)
#need to remove suffixes from CpG site names
row.names(beta_noob) <- sub("_.*", "", row.names(beta_noob))
head(rownames(beta_noob))

colnames(beta_noob)[1:5]
merged_demo2$Link[1:5]

beta_ids <- colnames(beta_noob)

sum(is.na(merged_demo2$Link))                        # how many missing Links
sum(merged_demo2$Link %in% beta_ids, na.rm = TRUE)   # how many Links found in beta_noob
length(setdiff(na.omit(merged_demo2$Link), beta_ids))# Links in metadata not in beta
length(setdiff(beta_ids, na.omit(merged_demo2$Link)))# beta samples not in metadata

meta_use <- merged_demo2 %>%
  filter(!is.na(Link)) %>%
  filter(Link %in% colnames(beta_noob))

beta_use <- beta_noob[, meta_use$Link, drop = FALSE]  # reorders columns to match metadata

identical(colnames(beta_use), meta_use$Link)

meta_use %>%
  count(Link) %>%
  filter(n > 1)
meta_use %>% count(Link) %>% filter(n > 1)


####################################################
#DNAm methylation calculation
dnam_age <- DNAmAge(beta_use)

str(dnam_age)
head(dnam_age)

names(meta_use)

meta_use2 <- meta_use %>%
  left_join(
    dnam_age %>% rename(Link = id),
    by = "Link"
  )

#Calculate DIFFERENCE based age acceleration
meta_use2 <- meta_use2 %>%
  mutate(
    Horvath_AA_diff     = Horvath     - Age,
    Hannum_AA_diff      = Hannum      - Age,
    Levine_AA_diff      = Levine      - Age,
    skinHorvath_AA_diff = skinHorvath - Age,
    PedBE_AA_diff       = PedBE       - Age,
    Wu_AA_diff          = Wu          - Age,
    BLUP_AA_diff        = BLUP        - Age,
    EN_AA_diff          = EN          - Age
  )

summary(meta_use2$Horvath_AA_diff)
hist(meta_use2$Horvath_AA_diff, main = "Horvath age acceleration (difference)")

#Calculate residual based age-adjusted age acceleration
meta_use2 <- meta_use2 %>%
  mutate(
    Horvath_AA_resid     = resid(lm(Horvath     ~ Age, data = .)),
    Hannum_AA_resid      = resid(lm(Hannum      ~ Age, data = .)),
    Levine_AA_resid      = resid(lm(Levine      ~ Age, data = .)),
    skinHorvath_AA_resid = resid(lm(skinHorvath ~ Age, data = .)),
    PedBE_AA_resid       = resid(lm(PedBE       ~ Age, data = .)),
    Wu_AA_resid          = resid(lm(Wu          ~ Age, data = .)),
    BLUP_AA_resid        = resid(lm(BLUP        ~ Age, data = .)),
    EN_AA_resid          = resid(lm(EN          ~ Age, data = .))
  )

summary(meta_use2$Horvath_AA_resid)

names(meta_use)

#add duned in pace 
pace_list <- PACEProjector(beta_noob)
pace_vec  <- pace_list$DunedinPACE

pace_df <- tibble::tibble(
  Link = names(pace_vec),
  DunedinPACE = as.numeric(pace_vec)
)

meta_use2 <- meta_use2 %>%
  left_join(pace_df, by = "Link")

nrow(meta_use2)


#################################################################################################



################################################################################

## SKIN PAH REGRESSIONS (Unadjusted and adjusted for years served)

####
#Calculate BMI and classify FF_status


dat_skin <- meta_use2 %>%
  mutate(
    BMI = `Weight (kg)` / (`Height (m)`^2),
    `Years served` = as.numeric(`Years served`)
  )

#treat FF_status as categorical
dat_skin <- dat_skin %>% mutate(FF_status = as.factor(FF_status))

###############################################
#  name exposures
#-------------------------
pah_vars <- names(dat_skin)[grepl("\\(ng\\)|^PAH", names(dat_skin))]
extra_vars <- c("FF_status", "Years served", "BMI")
exposure_vars <- c(pah_vars, extra_vars)

##############################################
# name outcomes for loop analyses 

outcomes <- c(
  "Horvath","Hannum","Levine","BNN","skinHorvath","PedBE","Wu","TL","BLUP","EN",
  "Horvath_AA_diff","Hannum_AA_diff","Levine_AA_diff","skinHorvath_AA_diff",
  "PedBE_AA_diff","Wu_AA_diff","BLUP_AA_diff","EN_AA_diff",
  "Horvath_AA_resid","Hannum_AA_resid","Levine_AA_resid","skinHorvath_AA_resid",
  "PedBE_AA_resid","Wu_AA_resid","BLUP_AA_resid","EN_AA_resid", "DunedinPACE"
)

# Keep only outcomes that exist in the data
outcomes <- outcomes[outcomes %in% names(dat_skin)]

###################################################

# model set up

# test model
run_one_model <- function(df, outcome, exposure) {
  # skip if outcome or exposure has no variation / all missing
  y <- df[[outcome]]
  x <- df[[exposure]]
  
  if (all(is.na(y)) || all(is.na(x))) return(NULL)
  if (dplyr::n_distinct(y, na.rm = TRUE) < 2) return(NULL)
  if (dplyr::n_distinct(x, na.rm = TRUE) < 2) return(NULL)
  
  fml <- as.formula(paste0("`", outcome, "` ~ `", exposure, "`"))
  fit <- lm(fml, data = df)
  
  broom::tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      n = stats::nobs(fit)
    ) %>%
    select(outcome, exposure, term, estimate, std.error, statistic, p.value, n)
}

run_models_grid <- function(df, outcomes, exposures, timepoint_label) {
  grid <- tidyr::crossing(outcome = outcomes, exposure = exposures)
  
  res <- pmap_dfr(
    list(grid$outcome, grid$exposure),
    ~ run_one_model(df, ..1, ..2)
  )
  
  res %>%
    mutate(timepoint = timepoint_label) %>%
    select(timepoint, everything())
}

#######################################################
# Split pre/post and run everything


dat_pre_skin  <- dat_skin %>% filter(Test == 1)
dat_post_skin <- dat_skin %>% filter(Test == 2)

res_pre_skin  <- run_models_grid(dat_pre_skin,  outcomes, exposure_vars, "Pre")
res_post_skin <- run_models_grid(dat_post_skin, outcomes, exposure_vars, "Post")

results_all_skin <- bind_rows(res_pre_skin, res_post_skin)

# View
results_all_skin %>% arrange(timepoint, outcome, p.value)

results_all_skin
#######
# Export

write.csv(results_all_skin, "SKIN_assoc_unadjusted_all_outcomes_pre_post.csv", row.names = FALSE)


####### 
## Adjusting for years of service

run_one_model_yrs <- function(df, outcome, exposure) {
  
  y <- df[[outcome]]
  x <- df[[exposure]]
  yrs <- df[["Years served"]]
  
  if (all(is.na(y)) || all(is.na(x)) || all(is.na(yrs))) return(NULL)
  if (dplyr::n_distinct(y, na.rm = TRUE) < 2) return(NULL)
  if (dplyr::n_distinct(x, na.rm = TRUE) < 2) return(NULL)
  
  fml <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + `Years served`"))
  fit <- lm(fml, data = df)
  
  broom::tidy(fit) %>%
    filter(term == exposure | term == paste0("`", exposure, "`")) %>%  # <-- FIX
    mutate(
      outcome = outcome,
      exposure = exposure,
      covariates = "Years served",
      n = stats::nobs(fit)
    ) %>%
    select(outcome, exposure, covariates,
           estimate, std.error, statistic, p.value, n)
}

run_models_grid_yrs <- function(df, outcomes, exposures, timepoint_label) {
  
  grid <- tidyr::crossing(outcome = outcomes, exposure = exposures)
  
  res <- purrr::pmap_dfr(
    list(grid$outcome, grid$exposure),
    ~ run_one_model_yrs(df, ..1, ..2)
  )
  
  res %>%
    mutate(timepoint = timepoint_label) %>%
    select(timepoint, everything())
}

dat_pre_skin  <- dat_skin %>% filter(Test == 1)
dat_post_skin <- dat_skin %>% filter(Test == 2)

res_pre_yrs_skin  <- run_models_grid_yrs(dat_pre_skin,  outcomes, exposure_vars, "Pre")
res_post_yrs_skin <- run_models_grid_yrs(dat_post_skin, outcomes, exposure_vars, "Post")

results_yrs_skin <- dplyr::bind_rows(res_pre_yrs_skin, res_post_yrs_skin)

results_yrs_skin

###
# Export

write.csv(
  results_yrs_skin,
  "SKIN_assoc_adjusted_years_served_pre_post.csv",
  row.names = FALSE
)


##############################################
#Model assumptions

run_model_and_diagnostics <- function(df, outcome, exposure) {
  
  fml <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + `Years served`"))
  fit <- lm(fml, data = df)
  
  ### Linearity: residuals vs fitted ###
  p_linearity <- ggplot(data.frame(fitted = fitted(fit), resid = resid(fit)),
                        aes(x = fitted, y = resid)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste("Residuals vs Fitted:", outcome, "~", exposure, "+ Years served"),
      x = "Fitted values",
      y = "Residuals"
    ) +
    theme_bw()
  
  ### Normality: QQ plot ###
  p_qq <- ggplot(data.frame(sample = resid(fit)), aes(sample = sample)) +
    stat_qq() +
    stat_qq_line() +
    labs(
      title = paste("Normal Q-Q:", outcome, "~", exposure, "+ Years served"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_bw()
  
  list(
    fit = fit,
    residuals_vs_fitted_plot = p_linearity,
    qq_plot = p_qq
  )
}

dat_pre_skin <- dat_skin %>% filter(Test == 1)

diag_obj <- run_model_and_diagnostics(
  df = dat_pre_skin,
  outcome = "Hannum_AA_resid",
  exposure = "PAHsumskin...9"   
)

# Show plots
diag_obj$residuals_vs_fitted_plot
diag_obj$qq_plot

# diagnostics
par(mfrow = c(1,2))
plot(diag_obj$fit, which = 1)  # Residuals vs Fitted
plot(diag_obj$fit, which = 2)  # QQ plot
par(mfrow = c(1,1))




################################################################################
################################################################################

################################################################################
## GEAR PAH: MERGE DATA THROUGH ANALYSES
################################################################################

names(demographics_0228)
names(gear_PAH_exposure)

# Merge demographics + gear PAH exposure
merged_demo_gear <- demographics_0228 %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    gear_PAH_exposure %>% mutate(ID = as.character(ID)),
    by = "ID"
  ) %>%
  tidyr::drop_na()

# Build link key from sample info
link_key <- OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  mutate(ID = as.character(ID)) %>%
  select(ID, VisitID, Link) %>%
  arrange(ID, VisitID) %>%
  group_by(ID, VisitID) %>%
  slice(1) %>%
  ungroup()

# Add VisitID from Test and join methylation Link
merged_demo2_gear <- merged_demo_gear %>%
  mutate(
    ID = as.character(ID),
    VisitID = case_when(
      Test == 1 ~ "V01",
      Test == 2 ~ "V02",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(link_key, by = c("ID", "VisitID"))

################################################################################
## CLOCK CALCULATION (beta_noob -> DNAmAge -> meta_use2_gear)
################################################################################

# Check structure
str(beta_noob)

# Remove suffixes from CpG site names
row.names(beta_noob) <- sub("_.*", "", row.names(beta_noob))
head(rownames(beta_noob))

colnames(beta_noob)[1:5]
merged_demo2_gear$Link[1:5]

beta_ids <- colnames(beta_noob)

sum(is.na(merged_demo2_gear$Link))                         # how many missing Links
sum(merged_demo2_gear$Link %in% beta_ids, na.rm = TRUE)   # how many Links found in beta_noob
length(setdiff(na.omit(merged_demo2_gear$Link), beta_ids))# Links in metadata not in beta
length(setdiff(beta_ids, na.omit(merged_demo2_gear$Link)))# beta samples not in metadata

meta_use_gear <- merged_demo2_gear %>%
  filter(!is.na(Link)) %>%
  filter(Link %in% colnames(beta_noob))

beta_use_gear <- beta_noob[, meta_use_gear$Link, drop = FALSE]

identical(colnames(beta_use_gear), meta_use_gear$Link)

meta_use_gear %>%
  count(Link) %>%
  filter(n > 1)

####################################################
# DNAm methylation calculation
dnam_age_gear <- DNAmAge(beta_use_gear)

str(dnam_age_gear)
head(dnam_age_gear)

names(meta_use_gear)

meta_use2_gear <- meta_use_gear %>%
  left_join(
    dnam_age_gear %>% rename(Link = id),
    by = "Link"
  )

# Calculate DIFFERENCE-based age acceleration
meta_use2_gear <- meta_use2_gear %>%
  mutate(
    Horvath_AA_diff     = Horvath     - Age,
    Hannum_AA_diff      = Hannum      - Age,
    Levine_AA_diff      = Levine      - Age,
    skinHorvath_AA_diff = skinHorvath - Age,
    PedBE_AA_diff       = PedBE       - Age,
    Wu_AA_diff          = Wu          - Age,
    BLUP_AA_diff        = BLUP        - Age,
    EN_AA_diff          = EN          - Age
  )

summary(meta_use2_gear$Horvath_AA_diff)
hist(meta_use2_gear$Horvath_AA_diff, main = "Horvath age acceleration (difference)")

# Calculate residual-based age-adjusted age acceleration
meta_use2_gear <- meta_use2_gear %>%
  mutate(
    Horvath_AA_resid     = resid(lm(Horvath     ~ Age, data = .)),
    Hannum_AA_resid      = resid(lm(Hannum      ~ Age, data = .)),
    Levine_AA_resid      = resid(lm(Levine      ~ Age, data = .)),
    skinHorvath_AA_resid = resid(lm(skinHorvath ~ Age, data = .)),
    PedBE_AA_resid       = resid(lm(PedBE       ~ Age, data = .)),
    Wu_AA_resid          = resid(lm(Wu          ~ Age, data = .)),
    BLUP_AA_resid        = resid(lm(BLUP        ~ Age, data = .)),
    EN_AA_resid          = resid(lm(EN          ~ Age, data = .))
  )

summary(meta_use2_gear$Horvath_AA_resid)

# Add DunedinPACE
pace_list <- PACEProjector(beta_noob)
pace_vec  <- pace_list$DunedinPACE

pace_df <- tibble::tibble(
  Link = names(pace_vec),
  DunedinPACE = as.numeric(pace_vec)
)

meta_use2_gear <- meta_use2_gear %>%
  left_join(pace_df, by = "Link")

nrow(meta_use2_gear)

################################################################################
## GEAR PAH REGRESSIONS (Unadjusted and adjusted for years served)
################################################################################

# Calculate BMI and classify FF_status
dat_gear <- meta_use2_gear %>%
  mutate(
    BMI = `Weight (kg)` / (`Height (m)`^2),
    `Years served` = as.numeric(`Years served`)
  )

# treat FF_status as categorical
dat_gear <- dat_gear %>%
  mutate(FF_status = as.factor(FF_status))

###############################################
# name exposures
pah_vars_gear <- names(dat_gear)[grepl("\\(ng\\)|^PAH", names(dat_gear))]
extra_vars <- c("FF_status", "Years served", "BMI")
exposure_vars_gear <- c(pah_vars_gear, extra_vars)

##############################################
# name outcomes for loop analyses
outcomes_gear <- c(
  "Horvath","Hannum","Levine","BNN","skinHorvath","PedBE","Wu","TL","BLUP","EN",
  "Horvath_AA_diff","Hannum_AA_diff","Levine_AA_diff","skinHorvath_AA_diff",
  "PedBE_AA_diff","Wu_AA_diff","BLUP_AA_diff","EN_AA_diff",
  "Horvath_AA_resid","Hannum_AA_resid","Levine_AA_resid","skinHorvath_AA_resid",
  "PedBE_AA_resid","Wu_AA_resid","BLUP_AA_resid","EN_AA_resid","DunedinPACE"
)

# Keep only outcomes that exist in the data
outcomes_gear <- outcomes_gear[outcomes_gear %in% names(dat_gear)]

###################################################
# model set up

run_one_model <- function(df, outcome, exposure) {
  y <- df[[outcome]]
  x <- df[[exposure]]
  
  if (all(is.na(y)) || all(is.na(x))) return(NULL)
  if (dplyr::n_distinct(y, na.rm = TRUE) < 2) return(NULL)
  if (dplyr::n_distinct(x, na.rm = TRUE) < 2) return(NULL)
  
  fml <- as.formula(paste0("`", outcome, "` ~ `", exposure, "`"))
  fit <- lm(fml, data = df)
  
  broom::tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      n = stats::nobs(fit)
    ) %>%
    select(outcome, exposure, term, estimate, std.error, statistic, p.value, n)
}

run_models_grid <- function(df, outcomes, exposures, timepoint_label) {
  grid <- tidyr::crossing(outcome = outcomes, exposure = exposures)
  
  res <- purrr::pmap_dfr(
    list(grid$outcome, grid$exposure),
    ~ run_one_model(df, ..1, ..2)
  )
  
  res %>%
    mutate(timepoint = timepoint_label) %>%
    select(timepoint, everything())
}

#######################################################
# Split pre/post and run everything

dat_pre_gear  <- dat_gear %>% filter(Test == 1)
dat_post_gear <- dat_gear %>% filter(Test == 2)

res_pre_gear  <- run_models_grid(dat_pre_gear,  outcomes_gear, exposure_vars_gear, "Pre")
res_post_gear <- run_models_grid(dat_post_gear, outcomes_gear, exposure_vars_gear, "Post")

results_all_gear <- bind_rows(res_pre_gear, res_post_gear)

# View
results_all_gear %>% arrange(timepoint, outcome, p.value)

# Export
write.csv(results_all_gear, "GEAR_assoc_unadjusted_all_outcomes_pre_post.csv", row.names = FALSE)

################################################################################
## Adjusting for years of service
################################################################################

run_one_model_yrs <- function(df, outcome, exposure) {
  
  y <- df[[outcome]]
  x <- df[[exposure]]
  yrs <- df[["Years served"]]
  
  if (all(is.na(y)) || all(is.na(x)) || all(is.na(yrs))) return(NULL)
  if (dplyr::n_distinct(y, na.rm = TRUE) < 2) return(NULL)
  if (dplyr::n_distinct(x, na.rm = TRUE) < 2) return(NULL)
  
  fml <- as.formula(paste0("`", outcome, "` ~ `", exposure, "` + `Years served`"))
  fit <- lm(fml, data = df)
  
  broom::tidy(fit) %>%
    filter(term == exposure | term == paste0("`", exposure, "`")) %>%
    mutate(
      outcome = outcome,
      exposure = exposure,
      covariates = "Years served",
      n = stats::nobs(fit)
    ) %>%
    select(outcome, exposure, covariates,
           estimate, std.error, statistic, p.value, n)
}

run_models_grid_yrs <- function(df, outcomes, exposures, timepoint_label) {
  
  grid <- tidyr::crossing(outcome = outcomes, exposure = exposures)
  
  res <- purrr::pmap_dfr(
    list(grid$outcome, grid$exposure),
    ~ run_one_model_yrs(df, ..1, ..2)
  )
  
  res %>%
    mutate(timepoint = timepoint_label) %>%
    select(timepoint, everything())
}

res_pre_yrs_gear  <- run_models_grid_yrs(dat_pre_gear,  outcomes_gear, exposure_vars_gear, "Pre")
res_post_yrs_gear <- run_models_grid_yrs(dat_post_gear, outcomes_gear, exposure_vars_gear, "Post")

results_yrs_gear <- dplyr::bind_rows(res_pre_yrs_gear, res_post_yrs_gear)

# View
results_yrs_gear

# Export
write.csv(
  results_yrs_gear,
  "GEAR_assoc_adjusted_years_served_pre_post.csv",
  row.names = FALSE
)


################################################################################

### Results output for tables and figures ###





################################################################################
#Identify subjects with demographics and CpG data (pre or post)

beta_links <- colnames(beta_noob)

cpg_subjects <- OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  filter(Link %in% beta_links) %>%
  distinct(ID)

subjects_with_demo_and_cpg <- OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  filter(Link %in% colnames(beta_noob)) %>%
  filter(ID %in% demographics_0228$ID) %>%
  distinct(ID)

nrow(subjects_with_demo_and_cpg)    # N=22

# N subject pre and N subjects post
OKF_EPIC2_SampleInfoFile_SampleInfo_ %>%
  filter(Link %in% colnames(beta_noob)) %>%
  filter(ID %in% demographics_0228$ID) %>%
  count(VisitID)



# Table 1 demographics

# Baseline-only dataset
dat_bl <- meta_use2_pass %>%
  filter(Test == 1) %>%
  mutate(
    BMI = `Weight (kg)` / (`Height (m)`^2),
    Gender = factor(Gender),
    Race   = factor(Race)
  )


dat_bl <- dat_bl %>%
  mutate(
    `Years served` = as.numeric(`Years served`),
    `Years served` = ifelse(
      is.na(`Years served`) | `Years served` < 0,
      NA,
      `Years served`
    )
  )

# FF_status formatting
dat_bl <- dat_bl %>%
  mutate(
    FF_status = recode(
      FF_status,
      "Career" = "Career",
      "Volunteer" = "Volunteer"
    ),
    FF_status = factor(FF_status, levels = c("Career", "Volunteer"))
  )

# Gender formatting
dat_bl <- dat_bl %>%
  mutate(
    Gender = recode(
      Gender,
      "F" = "Female",
      "M" = "Male"
    ),
    Gender = factor(Gender, levels = c("Female", "Male"))
  )

# Race formatting
dat_bl <- dat_bl %>%
  mutate(
    Race = recode(
      Race,
      "NHW" = "Non-Hispanic White"
    ),
    Race = factor(Race)
  )



# Variables for Table 1 (demographics + clock ages only)
vars_table1 <- c(
  "Age",
  "Years served",
  "BMI",
  "Gender",
  "Race",
  "FF_status",
  "Horvath",
  "Hannum",
  "Levine",
  "skinHorvath",
  "BLUP",
  "EN",
  "DunedinPACE"
)

# Create Table 1 (overall)
t1_demo <- dat_bl %>%
  select(any_of(vars_table1)) %>%
  tbl_summary(
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

t1_demo

# Export to Word
ft_tbl1 <- as_flex_table(t1_demo)

doc <- read_docx() %>%
  body_add_par(
    paste0("Table 1. Baseline demographics and DNAm age estimates ", nrow(dat_bl), ")"),
    style = "heading 1"
  ) %>%
  body_add_flextable(ft_tbl1)

print(doc, target = "Table1_Baseline_Demographics_Clocks.docx")


################################################################################
# Table 2 PAH stats between Skin and Gear

names(skin_PAH_exposure)
names(gear_PAH_exposure)

# analytes by panel
shared_analytes <- c(
  "Acenaphthene (ng)", "Fluorene (ng)", "Naphthalene (ng)",
  "2-Methyl-Naphthalene (ng)", "Phenanthrene (ng)", "Pyrene (ng)"
)

gear_only_analytes <- c(
  "Acenaphthylene (ng)", "Anthracene (ng)", "Fluoranthene (ng)",
  "1-Methyl-Naphthalene (ng)"
)

# Keep ONE sum variable per matrix (avoid duplicates in the table)
skin_sum_var <- "PAHsumskin...9"
gear_sum_var <- "PAHsumgear...13"

skin_vars <- c(shared_analytes, skin_sum_var)
gear_vars <- c(shared_analytes, gear_only_analytes, gear_sum_var)

# Skin baseline + Gear post-fire
skin_long <- skin_PAH_exposure %>%
  filter(Test == 1) %>%
  mutate(ID = as.character(ID)) %>%
  select(ID, all_of(skin_vars)) %>%
  pivot_longer(
    cols = all_of(skin_vars),
    names_to = "PAH",
    values_to = "Concentration_ng"
  ) %>%
  mutate(SampleType = "Skin")

gear_long <- gear_PAH_exposure %>%
  filter(Test == 2) %>%
  mutate(ID = as.character(ID)) %>%
  select(ID, all_of(gear_vars)) %>%
  pivot_longer(
    cols = all_of(gear_vars),
    names_to = "PAH",
    values_to = "Concentration_ng"
  ) %>%
  mutate(SampleType = "Gear")

pah_long <- bind_rows(skin_long, gear_long) %>%
  filter(!is.na(Concentration_ng)) %>%
  mutate(
    # Rename sum variables to reader-friendly labels
    PAH = recode(
      PAH,
      !!skin_sum_var := "Sum PAHs (skin)",
      !!gear_sum_var := "Sum PAHs (gear)"
    )
  )

sum_analytes <- c("Sum PAHs (skin)", "Sum PAHs (gear)")

pah_long <- pah_long %>%
  mutate(
    Section = case_when(
      PAH %in% shared_analytes ~ "PAHs collected in both skin and gear",
      PAH %in% gear_only_analytes ~ "PAHs collected in gear only",
      PAH %in% sum_analytes ~ "Summed PAH metrics",
      TRUE ~ "Other"
    )
  )

table2_df <- pah_long %>%
  group_by(Section, PAH, SampleType) %>%
  summarise(
    med = median(Concentration_ng, na.rm = TRUE),
    p25 = quantile(Concentration_ng, 0.25, na.rm = TRUE),
    p75 = quantile(Concentration_ng, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(`Median (IQR)` = sprintf("%.2f (%.2f, %.2f)", med, p25, p75)) %>%
  select(Section, PAH, SampleType, `Median (IQR)`) %>%
  pivot_wider(names_from = SampleType, values_from = `Median (IQR)`)

table2_df <- table2_df %>%
  mutate(
    PAH = factor(PAH, levels = c(shared_analytes, gear_only_analytes, sum_analytes)),
    Section = factor(
      Section,
      levels = c(
        "PAHs collected in both skin and gear",
        "PAHs collected in gear only",
        "Summed PAH metrics"
      )
    )
  ) %>%
  arrange(Section, PAH) %>%
  mutate(
    PAH = as.character(PAH),
    Section = as.character(Section)
  )

ft_tbl2 <- flextable(table2_df) %>%
  merge_v(j = "Section") %>%
  valign(j = "Section", valign = "top") %>%
  bold(j = "Section", bold = TRUE) %>%
  set_header_labels(
    Section = "",
    PAH = "PAH analyte",
    Skin = "Skin\nMedian (IQR), ng",
    Gear = "Gear\nMedian (IQR), ng"
  ) %>%
  theme_vanilla() %>%
  autofit()

ft_tbl2

doc <- read_docx() %>%
  body_add_par("Table 2. PAH concentrations by sample type", style = "heading 1") %>%
  body_add_par(
    "Values are median (IQR), ng. Skin = baseline wipe; Gear = post-fire turnout gear wipe.",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_tbl2)

print(doc, target = "Table2_PAH_Skin_vs_Gear.docx")


################################################################################
################################################################################

# Creating Supp. Table 1 (unadjusted estimates)

############################################

#Unadjusted (univariate results Skin (Pre)



skin_unadj <- results_all_skin

SKIN_TIMEPOINT <- "Pre"

keep_outcomes <- c(
  "Horvath_AA_resid",
  "Hannum_AA_resid",
  "Levine_AA_resid",
  "skinHorvath_AA_resid",
  "BLUP_AA_resid",
  "EN_AA_resid",
  "DunedinPACE"
)

outcome_labels <- c(
  Horvath_AA_resid     = "Horvath",
  Hannum_AA_resid      = "Hannum",
  Levine_AA_resid      = "Levine",
  skinHorvath_AA_resid = "SkinHorvath",
  BLUP_AA_resid        = "BLUP",
  EN_AA_resid          = "EN",
  DunedinPACE          = "DunedinPACE"
)

skin_pah_exposures <- c(
  "Acenaphthene (ng)",
  "Fluorene (ng)",
  "Naphthalene (ng)",
  "2-Methyl-Naphthalene (ng)",
  "Phenanthrene (ng)",
  "Pyrene (ng)",
  "PAHsumskin...9"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_cell <- function(est, se, p) {
  ifelse(
    is.na(est), "",
    paste0(
      sprintf("%.2f (%.2f)", est, se),
      "\n",
      "p = ", fmt_p(p)
    )
  )
}

pretty_exposure <- function(x) {
  x %>%
    str_replace(" \\(ng\\)", "") %>%
    recode("PAHsumskin...9" = "Sum PAHs (skin)")
}

term_matches_exposure <- function(term, exposure) {
  term == exposure | term == paste0("`", exposure, "`")
}

skin_tbl <- skin_unadj %>%
  filter(timepoint == SKIN_TIMEPOINT) %>%
  filter(outcome %in% keep_outcomes) %>%
  filter(exposure %in% skin_pah_exposures) %>%
  rowwise() %>%
  filter(term_matches_exposure(term, exposure)) %>%
  ungroup() %>%
  mutate(
    Clock = unname(outcome_labels[outcome]),
    exposure_pretty = pretty_exposure(exposure),
    cell = fmt_cell(estimate, std.error, p.value)
  ) %>%
  select(exposure, exposure_pretty, Clock, cell)

table_unadj <- skin_tbl %>%
  select(-exposure) %>%
  pivot_wider(
    names_from = Clock,
    values_from = cell
  ) %>%
  mutate(
    exposure_pretty = factor(exposure_pretty, levels = pretty_exposure(skin_pah_exposures))
  ) %>%
  arrange(exposure_pretty) %>%
  mutate(exposure_pretty = as.character(exposure_pretty))

ft_sup_tab1 <- flextable(table_unadj) %>%
  set_header_labels(exposure_pretty = "PAH analyte") %>%
  theme_vanilla() %>%
  padding(padding = 1, part = "all") %>%
  border_outer(border = fp_border(color = "black", width = 0.75)) %>%
  border_inner_v(border = fp_border(color = "black", width = 0.5)) %>%
  border_inner_h(border = fp_border(color = "black", width = 0.25)) %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Supplementary Table S1. Unadjusted associations between skin PAH concentrations and epigenetic aging outcomes",
    style = "heading 1"
  ) %>%
  body_add_par(
    "Cells show beta (SE); p-value for the PAH term. Skin results use the Pre timepoint only.",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_sup_tab1)

print(doc, target = "Table_S1_Unadjusted_Skin_PAH_EAA_beta_SE.docx")



#################################################################################
################################################################################

# Table 3 adjusted linear regression results table
# Adjusted for Years of service

names(results_yrs_skin)

#################################################################################
#################################################################################

# Table 3 adjusted linear regression results table -- SKIN ONLY
# Adjusted for Years served

skin_adj <- results_yrs_skin
SKIN_TIMEPOINT <- "Pre"

keep_outcomes <- c(
  "Horvath_AA_resid",
  "Hannum_AA_resid",
  "Levine_AA_resid",
  "skinHorvath_AA_resid",
  "BLUP_AA_resid",
  "EN_AA_resid",
  "DunedinPACE"
)

outcome_labels <- c(
  Horvath_AA_resid     = "Horvath",
  Hannum_AA_resid      = "Hannum",
  Levine_AA_resid      = "Levine",
  skinHorvath_AA_resid = "SkinHorvath",
  BLUP_AA_resid        = "BLUP",
  EN_AA_resid          = "EN",
  DunedinPACE          = "DunedinPACE"
)

skin_pah_exposures <- c(
  "Acenaphthene (ng)",
  "Fluorene (ng)",
  "Naphthalene (ng)",
  "2-Methyl-Naphthalene (ng)",
  "Phenanthrene (ng)",
  "Pyrene (ng)",
  "PAHsumskin...9"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_cell <- function(est, se, p) {
  ifelse(
    is.na(est), "",
    paste0(
      sprintf("%.2f (%.2f)", est, se),
      "\n",
      "p = ", fmt_p(p)
    )
  )
}

pretty_exposure <- function(x) {
  x %>%
    str_replace(" \\(ng\\)", "") %>%
    recode("PAHsumskin...9" = "Sum PAHs (skin)")
}

skin_tbl <- skin_adj %>%
  filter(timepoint == SKIN_TIMEPOINT) %>%
  filter(outcome %in% keep_outcomes) %>%
  filter(exposure %in% skin_pah_exposures) %>%
  mutate(
    Clock = unname(outcome_labels[outcome]),
    exposure_pretty = pretty_exposure(exposure),
    cell = fmt_cell(estimate, std.error, p.value)
  ) %>%
  select(exposure, exposure_pretty, Clock, cell)

table_adj <- skin_tbl %>%
  select(-exposure) %>%
  pivot_wider(
    names_from = Clock,
    values_from = cell
  ) %>%
  mutate(
    exposure_pretty = factor(exposure_pretty, levels = pretty_exposure(skin_pah_exposures))
  ) %>%
  arrange(exposure_pretty) %>%
  mutate(exposure_pretty = as.character(exposure_pretty))

ft_tbl3 <- flextable(table_adj) %>%
  set_header_labels(exposure_pretty = "PAH analyte") %>%
  theme_vanilla() %>%
  padding(padding = 1, part = "all") %>%
  border_outer(border = fp_border(color = "black", width = 0.75)) %>%
  border_inner_v(border = fp_border(color = "black", width = 0.5)) %>%
  border_inner_h(border = fp_border(color = "black", width = 0.25)) %>%
  line_spacing(space = 1, part = "body") %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Table 3. Associations between skin PAH concentrations and epigenetic aging outcomes (adjusted for years served)",
    style = "heading 1"
  ) %>%
  body_add_par(
    "Cells show beta (SE) on the first line and p-value on the second line for the PAH term. Skin results use the Pre timepoint only. Models adjusted for years served.",
    style = "Normal"
  ) %>%
  body_add_flextable(ft_tbl3)

print(doc, target = "Table3_Adjusted_Skin_PAH_EAA_beta_SE.docx")




##################################################################
##################################################################
# Making scatter plots

# Making scatter plots

# Assign dataset
dat_skin_full <- dat_skin


make_partial_df <- function(df, x, y, adj = "Years served") {
  x_res <- resid(lm(as.formula(paste0("`", x, "` ~ `", adj, "`")),
                    data = df, na.action = na.exclude))
  y_res <- resid(lm(as.formula(paste0("`", y, "` ~ `", adj, "`")),
                    data = df, na.action = na.exclude))
  
  out <- df %>%
    transmute(
      x_raw = .data[[x]],
      y_raw = .data[[y]],
      yrs   = .data[[adj]]
    )
  
  out$x_adj <- as.numeric(x_res)
  out$y_adj <- as.numeric(y_res)
  
  out %>% filter(!is.na(x_adj) & !is.na(y_adj))
}


plot_partial <- function(df, x, y, title, xlab, ylab) {
  d <- make_partial_df(df, x = x, y = y, adj = "Years served")
  
  ggplot(d, aes(x = x_adj, y = y_adj)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}


skin_pre <- dat_skin_full %>% filter(Test == 1)

#####
# Two panels

p1 <- plot_partial(
  skin_pre,
  x = "PAHsumskin...9",
  y = "Hannum_AA_resid",
  title = "A. Sum Skin PAHs vs Hannum AA",
  xlab = "Sum Skin PAHs, adjusted for Years served",
  ylab = "Hannum AA residual"
)

p2 <- plot_partial(
  skin_pre,
  x = "PAHsumskin...9",
  y = "DunedinPACE",
  title = "B. Sum Skin PAHs vs DunedinPACE",
  xlab = "Sum Skin PAHs, adjusted for Years served",
  ylab = "DunedinPACE, adjusted for Years served"
)

fig_sumskin <- p1 | p2 +
  plot_annotation(
    title = "Figure 1. Adjusted associations of summed skin PAH exposure with epigenetic aging outcomes"
  )

# Export

ggsave("Figure1_SumSkin_Hannum_DunedinPACE.png",
       fig_sumskin,
       width = 10,
       height = 5,
       dpi = 400)



#################################################################################
#################################################################################

# Forest plot for the "Sum" variabnles in both skin and gear across all clocks

# adjusted results
skin_adj_res <- results_yrs_skin

#########################
# choose clocks

keep_outcomes <- c(
  "Horvath_AA_resid",
  "Hannum_AA_resid",
  "Levine_AA_resid",
  "skinHorvath_AA_resid",
  "BLUP_AA_resid",
  "EN_AA_resid",
  "DunedinPACE"
)

clock_labels <- c(
  Horvath_AA_resid     = "Horvath",
  Hannum_AA_resid      = "Hannum",
  Levine_AA_resid      = "Levine",
  skinHorvath_AA_resid = "SkinHorvath",
  BLUP_AA_resid        = "BLUP",
  EN_AA_resid          = "EN",
  DunedinPACE          = "DunedinPACE"
)

#########################
# Skin summed PAHs only

forest_df <- skin_adj_res %>%
  filter(timepoint == "Pre") %>%
  filter(grepl("^PAHsumskin", exposure)) %>%
  filter(outcome %in% keep_outcomes) %>%
  mutate(
    Clock = unname(clock_labels[outcome]),
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    significant = p.value < 0.05,
    Clock = factor(
      Clock,
      levels = c("Horvath","Hannum","Levine","SkinHorvath","BLUP","EN","DunedinPACE")
    )
  )

#########################
# Forest plot

p <- ggplot(forest_df, aes(x = estimate, y = Clock, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2,
    color = "gray50"
  ) +
  
  geom_point(size = 3) +
  
  scale_color_manual(
    values = c(
      "TRUE"  = "#08519C",   # darker blue
      "FALSE" = "#9ECAE1"    # light blue
    )
  ) +
  
  labs(
    title = "Adjusted associations of summed skin PAH exposure with epigenetic aging across clocks",
    x = "Beta (95% CI), adjusted for Years served",
    y = ""
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  "Figure2_Forest_SumSkinPAH_AcrossClocks.png",
  p,
  width = 7,
  height = 5,
  dpi = 400
)





###############################################################################################
###############################################################################################

# Leave one out approach to see which analytes are the most impactful

#For Skin

dat_pre_skin  <- dat_skin %>% filter(Test == 1) 
dat_post_skin <- dat_skin %>% filter(Test == 2)

df_pre <- dat_pre_skin



outcome <- "Hannum_AA_resid"
sum_var <- "PAHsumskin...9"   
covar   <- "Years served"

analytes <- c(
  "Acenaphthene (ng)",
  "Fluorene (ng)",
  "Naphthalene (ng)",
  "2-Methyl-Naphthalene (ng)",
  "Phenanthrene (ng)",
  "Pyrene (ng)"
)

full_fit <- lm(
  as.formula(paste0("`", outcome, "` ~ `", sum_var, "` + `", covar, "`")),
  data = df_pre
)

summary(full_fit)
beta_full <- coef(full_fit)[sum_var]
r2_full   <- summary(full_fit)$r.squared
p_full    <- summary(full_fit)$coefficients[sum_var, "Pr(>|t|)"]
n_full    <- nobs(full_fit)

###########################################
#run regression LOO analysis

loo_res <- map_dfr(analytes, function(a){
  
  df_tmp <- df_pre %>%
    mutate(sum_loo = .data[[sum_var]] - .data[[a]])
  
  fit <- lm(
    as.formula(paste0("`", outcome, "` ~ sum_loo + `", covar, "`")),
    data = df_tmp
  )
  
  beta_loo <- coef(fit)["sum_loo"]
  r2_loo   <- summary(fit)$r.squared
  p_loo    <- summary(fit)$coefficients["sum_loo", "Pr(>|t|)"]
  
  tibble(
    left_out = a,
    beta_full = beta_full,
    beta_loo = beta_loo,
    delta_beta = beta_loo - beta_full,
    pct_atten = 100 * (beta_full - beta_loo) / abs(beta_full),
    r2_full = r2_full,
    r2_loo = r2_loo,
    delta_r2 = r2_loo - r2_full,
    p_full = p_full,
    p_loo = p_loo,
    n = n_full
  )
}) %>% arrange(desc(pct_atten))

loo_res


###########################################################
# Create figure to show influence of each individual analyte on hannum clock acceleration

loo_plot_df <- loo_res %>%
  mutate(
    direction = ifelse(pct_atten >= 0,
                       "Primary contributor (removal reduces β) ",
                       "Offestting contributor (removal increases β)"),
    left_out = fct_reorder(left_out, pct_atten)
  )

hannum_LOO <- ggplot(loo_plot_df, aes(x = pct_atten, y = left_out, fill = direction)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0) +
  labs(
    title = "Leave-one-out analyte contribution (skin)\nOutcome: Hannum; Model: sumPAH + Years served",
    x = "% attenuation of summed-PAH beta when analyte is removed",
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

hannum_LOO

ggsave(
  filename = "LOO_significant_skin_Hannum_contribution.jpg",
  plot = hannum_LOO,
  width = 5,
  height = 7,
  units = "in",
  dpi = 600
)



########################################

# Same figure but now for DunedinPACE summed PAH skin exposure 



df_pre <- dat_pre_skin

outcome <- "DunedinPACE"     
sum_var <- "PAHsumskin...9"
covar   <- "Years served"

analytes <- c(
  "Acenaphthene (ng)",
  "Fluorene (ng)",
  "Naphthalene (ng)",
  "2-Methyl-Naphthalene (ng)",
  "Phenanthrene (ng)",
  "Pyrene (ng)"
)

full_fit <- lm(
  as.formula(paste0("`", outcome, "` ~ `", sum_var, "` + `", covar, "`")),
  data = df_pre
)

summary(full_fit)
beta_full <- coef(full_fit)[sum_var]
r2_full   <- summary(full_fit)$r.squared
p_full    <- summary(full_fit)$coefficients[sum_var, "Pr(>|t|)"]
n_full    <- nobs(full_fit)

###########################################
# run regression LOO analysis

loo_res <- map_dfr(analytes, function(a){
  
  df_tmp <- df_pre %>%
    mutate(sum_loo = .data[[sum_var]] - .data[[a]])
  
  fit <- lm(
    as.formula(paste0("`", outcome, "` ~ sum_loo + `", covar, "`")),
    data = df_tmp
  )
  
  beta_loo <- coef(fit)["sum_loo"]
  r2_loo   <- summary(fit)$r.squared
  p_loo    <- summary(fit)$coefficients["sum_loo", "Pr(>|t|)"]
  
  tibble(
    left_out  = a,
    beta_full = beta_full,
    beta_loo  = beta_loo,
    delta_beta = beta_loo - beta_full,
    pct_atten = 100 * (beta_full - beta_loo) / abs(beta_full),
    r2_full   = r2_full,
    r2_loo    = r2_loo,
    delta_r2  = r2_loo - r2_full,
    p_full    = p_full,
    p_loo     = p_loo,
    n         = n_full
  )
}) %>% arrange(desc(pct_atten))

loo_res


###########################################################
# Create figure: DunedinPACE LOO

loo_plot_df <- loo_res %>%
  mutate(
    direction = ifelse(
      pct_atten >= 0,
      "Primary contributor (removal reduces β)",
      "Offsetting contributor (removal increases β)"
    ),
    left_out = fct_reorder(left_out, pct_atten)
  )

pace_LOO <- ggplot(loo_plot_df, aes(x = pct_atten, y = left_out, fill = direction)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0) +
  labs(
    title = "Leave-one-out analyte contribution (skin)\nOutcome: DunedinPACE; Model: sumPAH + Years served",
    x = "% attenuation of summed-PAH beta when analyte is removed",
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

pace_LOO

ggsave(
  filename = "LOO_skin_DunedinPACE_contribution.jpg",
  plot = pace_LOO,
  width = 5,
  height = 7,
  units = "in",
  dpi = 600
)


##################################################################################################
##################################################################################################

# Sensitivity analyses for cell type proportions

names(df_pre)


# Cell proportions
cells <- OKF_EPIC2_houseman_cbc %>%
  rename(Link = `...1`) %>%
  mutate(Link = as.character(Link))

# Merge cells into skin dataset
dat_skin_cells <- dat_skin %>%
  mutate(Link = as.character(Link)) %>%
  left_join(cells, by = "Link")

# check
names(dat_skin_cells)

df_pre <- dat_skin_cells %>% filter(Test == 1)

sum_var <- "PAHsumskin...9"
cov_yrs <- "Years served"

out_hannum <- "Hannum_AA_resid"
out_pace   <- "DunedinPACE"

# Years served only
fit_hannum_yrs <- lm(
  as.formula(paste0("`", out_hannum, "` ~ `", sum_var, "` + `", cov_yrs, "`")),
  data = df_pre
)

fit_pace_yrs <- lm(
  as.formula(paste0("`", out_pace, "` ~ `", sum_var, "` + `", cov_yrs, "`")),
  data = df_pre
)

# Years served + cells
cell_terms <- "CD8T + CD4T + NK + Bcell + Mono + Gran"

fit_hannum_cells <- lm(
  as.formula(paste0("`", out_hannum, "` ~ `", sum_var, "` + `", cov_yrs, "` + ", cell_terms)),
  data = df_pre
)

fit_pace_cells <- lm(
  as.formula(paste0("`", out_pace, "` ~ `", sum_var, "` + `", cov_yrs, "` + ", cell_terms)),
  data = df_pre
)

summary(fit_hannum_yrs)
summary(fit_hannum_cells)

summary(fit_pace_yrs)
summary(fit_pace_cells)



get_pah_term <- function(fit, term = "PAHsumskin...9") {
  sm <- summary(fit)$coefficients
  out <- sm[term, c("Estimate", "Std. Error", "Pr(>|t|)")]
  data.frame(
    beta = out[["Estimate"]],
    se   = out[["Std. Error"]],
    p    = out[["Pr(>|t|)"]],
    n    = nobs(fit),
    row.names = NULL
  )
}

cell_sensitivity_results <- 
  bind_rows(
  cbind(outcome = out_hannum, model = "Years served",         get_pah_term(fit_hannum_yrs,   sum_var)),
  cbind(outcome = out_hannum, model = "Years served + cells", get_pah_term(fit_hannum_cells, sum_var)),
  cbind(outcome = out_pace,   model = "Years served",         get_pah_term(fit_pace_yrs,     sum_var)),
  cbind(outcome = out_pace,   model = "Years served + cells", get_pah_term(fit_pace_cells,   sum_var))
)


cell_sensitivity_results <- cell_sensitivity_results %>%
  mutate(
    beta = round(beta, 3),
    se   = round(se, 3),
    p    = ifelse(p < 0.001, "<0.001", round(p, 3))
  )

ft_cells <- flextable(cell_sensitivity_results) %>%
  set_header_labels(
    outcome = "Outcome",
    model   = "Model",
    beta    = "Beta",
    se      = "SE",
    p       = "p-value",
    n       = "N"
  ) %>%
  theme_vanilla() %>%
  autofit()


#####
# Export

doc <- read_docx() %>%
  body_add_par(
    "Supplementary Table S3. Sensitivity analyses adjusting for leukocyte composition",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft_cells)

print(doc, target = "Table_S2_Cell_Adjusted_Sensitivity.docx")


















