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
  select(ASBGRST1:ASBGRST6, ASRREA01:ASRREA05)

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

# UPI
mod_upi <- "
          IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
          EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode

          ASRREA02 ~ b1*IM + b2*EM + b3*IM:EM

          IM ~~ v1*IM
          EM ~~ v2*EM
          beta1 := b1*sqrt(v1)
          beta2 := b2*sqrt(v2)
          beta3 := b3*sqrt(v1)*sqrt(v2)
          "
# Fit the model using upi()
fit_upi <- upi(mod_upi, PIRLS_Data, mode = "match")
summary(fit_upi, fit.measure = TRUE, standardized = TRUE)

# RAPI
mod_rapi <- '
          IM =~ ASBGRST1 + ASBGRST4 + ASBGRST6_recode
          EM =~ ASBGRST2_recode + ASBGRST3_recode + ASBGRST5_recode

          ASRREA02 ~ b1*IM + b2*EM + b3*IM:EM

          IM ~~ v1*IM
          EM ~~ v2*EM
          beta1 := b1*sqrt(v1)
          beta2 := b2*sqrt(v2)
          beta3 := b3*sqrt(v1)*sqrt(v2)
  '
rapi_fit <- rapi(data = PIRLS_Data,
                 model = mod_rapi)
summary(rapi_fit, fit.measure = TRUE, standardized = TRUE)

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
ReadingScore <- PIRLS_Data$ASRREA05
fs_dat <- cbind(fs_dat, ReadingScore)

fit_tspa <- tspa(model = "ReadingScore ~ b1*IM + b2*EM + b3*IM:EM
                          beta1 := b1 * sqrt(v1)
                          beta2 := b2 * sqrt(v2)
                          beta3 := b3 * sqrt(v1) * sqrt(v2)",
                 data = fs_dat,
                 se = list(IM = fs_dat$fs_IM_se[1],
                           EM = fs_dat$fs_EM_se[1]))
summary(fit_tspa, fit.measure = TRUE, standardized = TRUE)

# Write data for Mplus
write.table(PIRLS_Data, file = here("Paper", "PIRLS2006", "PIRLS_Data.dat"),
            sep = " ", row.names = FALSE, col.names = FALSE, na = ".")