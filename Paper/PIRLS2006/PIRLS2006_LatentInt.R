library(haven)
library(dplyr)
library(tidyverse)
library(lavaan)
library(here)
library(semTools)

# Import Data
PIRLS2006 <- read_sav(here("Paper", "PIRLS2006", "asgusar2.sav"))
# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# Data Preprocessing
PIRLS_Data <- PIRLS2006 %>%
  select(ASBGRST1:ASBGRST6, ASRREA01:ASRREA05) %>%
  drop_na()

PIRLS_Data <- PIRLS_Data %>%
  mutate(across(
    .cols = c(ASBGRST2, ASBGRST3, ASBGRST5, ASBGRST6), # replace with your actual column names
    .fns = ~ recode(as.numeric(.x), `1` = 4, `2` = 3, `3` = 2, `4` = 1),
    .names = "{.col}_recode"
  )) %>%
  drop_na()

# Measurement Model
IM <- PIRLS_Data %>% select(ASBGRST1, ASBGRST4, ASBGRST6_recode)
EM <- PIRLS_Data %>% select(ASBGRST2_recode, ASBGRST3_recode, ASBGRST5_recode)

model <- "IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
          EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode"

summary(sem(model, PIRLS_Data))
fitmeasures(sem(model, PIRLS_Data))

# Dependent Vars
DVs <- c("ASRREA01", "ASRREA02", "ASRREA03", "ASRREA04", "ASRREA05")

# UPI
mod_upi <- "
          IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
          EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode

          DV ~ b1*IM + b2*EM + b3*IM:EM

          IM ~~ v1*IM
          EM ~~ v2*EM
          beta1 := b1*sqrt(v1)
          beta2 := b2*sqrt(v2)
          beta3 := b3*sqrt(v1)*sqrt(v2)
          "
# Fit the model using upi()
beta1_upi <- c()
beta2_upi <- c()
beta3_upi <- c()
se1_upi <- c()
se2_upi <- c()
se3_upi <- c()
p1_upi <- c()
p2_upi <- c()
p3_upi <- c()

for (dv in DVs) {
  mod_upi_temp <- gsub("DV", dv, mod_upi)  # Replace "DV" with the current dv in the loop
  fit_upi <- upi(mod_upi_temp, PIRLS_Data, mode = "match")
  beta1_upi <- unname(c(beta1_upi, coef(fit_upi, type = "user")["beta1"]))
  beta2_upi <- unname(c(beta2_upi, coef(fit_upi, type = "user")["beta2"]))
  beta3_upi <- unname(c(beta3_upi, coef(fit_upi, type = "user")["beta3"]))
  se1_upi <- c(se1_upi, sqrt(vcov(fit_upi, type = "user")["beta1", "beta1"]))
  se2_upi <- c(se2_upi, sqrt(vcov(fit_upi, type = "user")["beta2", "beta2"]))
  se3_upi <- c(se3_upi, sqrt(vcov(fit_upi, type = "user")["beta3", "beta3"]))
  par_upi <- parameterEstimates(fit_upi)
  p1_upi <- c(p1_upi, as.integer(par_upi %>% filter(label == "beta1") %>% select(pvalue)))
  p2_upi <- c(p2_upi, as.integer(par_upi %>% filter(label == "beta2") %>% select(pvalue)))
  p3_upi <- c(p3_upi, as.integer(par_upi %>% filter(label == "beta3") %>% select(pvalue)))
}

beta1_upi_avg <- mean(beta1_upi)
beta2_upi_avg <- mean(beta2_upi)
beta3_upi_avg <- mean(beta3_upi)

combine_se <- function(combined_path, se_val) {
  within_var <- mean(se_val^2)
  between_var <- var(combined_path)
  m <- length(combined_path)
  total_variance <- within_var + (1 + 1/m) * between_var
  return(sqrt(total_variance))
}

se1_upi_avg <- combine_se(beta1_upi, se1_upi)
se2_upi_avg <- combine_se(beta2_upi, se2_upi)
se3_upi_avg <- combine_se(beta3_upi, se3_upi)

combine_p_value <- function(avg_path, combined_se) {
  combined_z <- avg_path / combined_se
  combined_p_value <- 2 * (1 - pnorm(abs(combined_z)))
  return(combined_p_value)
}

p1_upi_avg <- combine_p_value(beta1_upi_avg, se1_upi_avg)
p2_upi_avg <- combine_p_value(beta2_upi_avg, se2_upi_avg)
p3_upi_avg <- combine_p_value(beta3_upi_avg, se3_upi_avg)

# RAPI
mod_rapi <- "
          IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
          EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode

          DV ~ b1*IM + b2*EM + b3*IM:EM

          IM ~~ v1*IM
          EM ~~ v2*EM
          beta1 := b1*sqrt(v1)
          beta2 := b2*sqrt(v2)
          beta3 := b3*sqrt(v1)*sqrt(v2)
          "
# Fit the model using upi()
beta1_rapi <- c()
beta2_rapi <- c()
beta3_rapi <- c()
se1_rapi <- c()
se2_rapi <- c()
se3_rapi <- c()
p1_rapi <- c()
p2_rapi <- c()
p3_rapi <- c()

for (dv in DVs) {
  mod_rapi_temp <- gsub("DV", dv, mod_rapi)  # Replace "DV" with the current dv in the loop
  fit_rapi <- rapi(mod_rapi_temp, PIRLS_Data)
  beta1_rapi <- c(beta1_rapi, coef(fit_rapi, type = "user")["beta1"])
  beta2_rapi <- c(beta2_rapi, coef(fit_rapi, type = "user")["beta2"])
  beta3_rapi <- c(beta3_rapi, coef(fit_rapi, type = "user")["beta3"])
  se1_rapi <- c(se1_rapi, sqrt(vcov(fit_rapi, type = "user")["beta1", "beta1"]))
  se2_rapi <- c(se2_rapi, sqrt(vcov(fit_rapi, type = "user")["beta2", "beta2"]))
  se3_rapi <- c(se3_rapi, sqrt(vcov(fit_rapi, type = "user")["beta3", "beta3"]))
  par_rapi <- parameterEstimates(fit_rapi)
  p1_rapi <- c(p1_rapi, as.integer(par_rapi %>% filter(label == "beta1") %>% select(pvalue)))
  p2_rapi <- c(p2_rapi, as.integer(par_rapi %>% filter(label == "beta2") %>% select(pvalue)))
  p3_rapi <- c(p3_rapi, as.integer(par_rapi %>% filter(label == "beta3") %>% select(pvalue)))
}

beta1_rapi_avg <- mean(beta1_rapi)
beta2_rapi_avg <- mean(beta2_rapi)
beta3_rapi_avg <- mean(beta3_rapi)

se1_rapi_avg <- combine_se(beta1_rapi, se1_rapi)
se2_rapi_avg <- combine_se(beta2_rapi, se2_rapi)
se3_rapi_avg <- combine_se(beta3_rapi, se3_rapi)

p1_rapi_avg <- combine_p_value(beta1_rapi_avg, se1_rapi_avg)
p2_rapi_avg <- combine_p_value(beta2_rapi_avg, se2_rapi_avg)
p3_rapi_avg <- combine_p_value(beta3_rapi_avg, se3_rapi_avg)

# 2SPA
mod_cfa <- '
IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode
'
# Obtain factor scores
fs_dat <- get_fs(PIRLS_Data, model = mod_cfa, method = "Bartlett", std.lv = TRUE)
# Obtain factor product
fs_dat$fs_int <- fs_dat$fs_IM * fs_dat$fs_EM
# Centering the product variable
fs_dat$fs_int <- fs_dat$fs_int - mean(fs_dat$fs_int)
fs_dat$fs_int_se <- sqrt(1 * fs_dat$fs_IM_se^2 + 1 * fs_dat$fs_EM_se^2 + fs_dat$fs_IM_se^2*fs_dat$fs_EM_se^2)

beta1_tspa <- c()
beta2_tspa <- c()
beta3_tspa <- c()
se1_tspa <- c()
se2_tspa <- c()
se3_tspa <- c()
p1_tspa <- c()
p2_tspa <- c()
p3_tspa <- c()

for (dv in DVs) {
  if ("ReadingScore" %in% names(fs_dat)) {
    fs_dat$ReadingScore <- NULL
  }
  fs_dat$ReadingScore <- PIRLS_Data[[dv]]
  fit_tspa <- tspa(model = "ReadingScore ~ b1*IM + b2*EM + b3*IM:EM
                          beta1 := b1 * sqrt(v1)
                          beta2 := b2 * sqrt(v2)
                          beta3 := b3 * sqrt(v1) * sqrt(v2)",
                   data = fs_dat,
                   se = list(IM = fs_dat$fs_IM_se[1],
                             EM = fs_dat$fs_EM_se[1]))
  beta1_tspa <- c(beta1_tspa, coef(fit_tspa, type = "user")["beta1"])
  beta2_tspa <- c(beta2_tspa, coef(fit_tspa, type = "user")["beta2"])
  beta3_tspa <- c(beta3_tspa, coef(fit_tspa, type = "user")["beta3"])
  se1_tspa <- c(se1_tspa, sqrt(vcov(fit_tspa, type = "user")["beta1", "beta1"]))
  se2_tspa <- c(se2_tspa, sqrt(vcov(fit_tspa, type = "user")["beta2", "beta2"]))
  se3_tspa <- c(se3_tspa, sqrt(vcov(fit_tspa, type = "user")["beta3", "beta3"]))
  par_tspa <- parameterEstimates(fit_tspa)
  p1_tspa <- c(p1_tspa, as.integer(par_tspa %>% filter(label == "beta1") %>% select(pvalue)))
  p2_tspa <- c(p2_tspa, as.integer(par_tspa %>% filter(label == "beta2") %>% select(pvalue)))
  p3_tspa <- c(p3_tspa, as.integer(par_tspa %>% filter(label == "beta3") %>% select(pvalue)))
}

beta1_tspa_avg <- mean(beta1_tspa)
beta2_tspa_avg <- mean(beta2_tspa)
beta3_tspa_avg <- mean(beta3_tspa)

se1_tspa_avg <- combine_se(beta1_tspa, se1_tspa)
se2_tspa_avg <- combine_se(beta2_tspa, se2_tspa)
se3_tspa_avg <- combine_se(beta3_tspa, se3_tspa)

p1_tspa_avg <- combine_p_value(beta1_tspa_avg, se1_tspa_avg)
p2_tspa_avg <- combine_p_value(beta2_tspa_avg, se2_tspa_avg)
p3_tspa_avg <- combine_p_value(beta3_tspa_avg, se3_tspa_avg)