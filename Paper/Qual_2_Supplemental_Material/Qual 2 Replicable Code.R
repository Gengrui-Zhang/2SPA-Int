# Import dataset
TA2019 <- zap_formats(zap_labels(read_sav("~/TA2019.sav")))

# Scale 1: ROSENBURG SELF-ESTEEM SCALE
# Recoded items: TA190106, TA190108, TA190111, TA190112, TA190113
# Strongly agree - 4; Agree - 3; Disagree - 2; Strongly Disagree - 1; Missing - 9.
SelfE <- TA2019 %>%
  select(TA190104:TA190113) %>%
  mutate(TA190106 = recode(TA190106, "1" = "4", "2" = "3", "3" = "2", "4" = "1", "9" = "9"),
         TA190108 = recode(TA190108, "1" = "4", "2" = "3", "3" = "2", "4" = "1", "9" = "9"),
         TA190111 = recode(TA190111, "1" = "4", "2" = "3", "3" = "2", "4" = "1", "9" = "9"),
         TA190112 = recode(TA190112, "1" = "4", "2" = "3", "3" = "2", "4" = "1", "9" = "9"),
         TA190113 = recode(TA190113, "1" = "4", "2" = "3", "3" = "2", "4" = "1", "9" = "9")) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_all(na_if, 9)
colnames(SelfE) <- paste0("SelfE", 1:10)

# Scale 2: Everyday Discrimination
# No recoded items
# Never - 1; Less than once a year - 2; A few times a year - 3; A few times a month - 4;
# At least once a week - 5; Almost every day - 6; Missing - 9.
PED <- TA2019 %>%
  select(TA192066:TA192072) %>%
  mutate_all(na_if, 9) %>%
  mutate_all(na_if, 8)
colnames(PED) <- paste0("PED", 1:7)

# Scale 3: PHQ-9 Depression Scale
# No recoded items
# Not at all - 1; Several days - 2; More than half the days - 3; Nearly every day - 4; Missing - 9.
PHQ <- TA2019 %>%
  select(TA190114:TA190122) %>%
  mutate_all(na_if, 9) %>%
  mutate_all(na_if, 8)
colnames(PHQ) <- paste0("PHQ", 1:9)

# Dimension of dat: 2,595 observations and 26 first-order indicators
dat <- cbind(PED, SelfE, PHQ)
# Mean-centering first-order indicators of PED and SelfE
dat.centered <- dat %>%
  mutate(across(.cols = everything(), .fns = ~ .x - mean(.x, na.rm = TRUE)))
dat.matchpair <- indProd(dat.centered,
                         var1 = c("PED6", "PED3", "PED7", "PED1", "PED5", "PED2", "PED4"),
                         var2 = c("SelfE10", "SelfE9", "SelfE6", "SelfE7", "SelfE5", "SelfE3", "SelfE8"),
                         match = TRUE,
                         meanC = FALSE,
                         residualC = FALSE,
                         doubleMC = TRUE)
# Model Specification
model.matchpair <- "# Measurement model
                      PHQ =~ PHQ1 + PHQ2 + PHQ3 + PHQ4 + PHQ5 + PHQ6 + PHQ7 + PHQ8 + PHQ9
                      PED =~ PED6 + PED3 + PED7 + PED1 + PED5 + PED2 + PED4
                      SelfE =~ SelfE10 + SelfE9 + SelfE6 + SelfE7 + SelfE5 + SelfE3 + SelfE8
                      PED.SelfE =~ PED6.SelfE10 + PED3.SelfE9 + PED7.SelfE6 + PED1.SelfE7 +
                                   PED5.SelfE5 + PED2.SelfE3 + PED4.SelfE8
                    # Latent variance
                      PED ~~ v1*PED
                      SelfE ~~ v2*SelfE
                      PED.SelfE ~~ v3*PED.SelfE
                    # Latent covariance
                      PED ~~ v12*SelfE
                      PED ~~ v13*PED.SelfE
                      SelfE ~~ v23*PED.SelfE
                    # Residual variance of DV
                      PHQ ~~ v4*PHQ
                    # Structural model
                      PHQ ~ g1*PED + g2*SelfE + g3*PED.SelfE
                    # Standardized
                      v_y := g1^2*v1 + g2^2*v2 + g3^2*v3 + 2*g1*g2*v12 +
                             2*g1*g3*v13 + 2*g2*g3*v23 + v4
                      gamma1 := g1*sqrt(v1)/sqrt(v_y)
                      gamma2 := g2*sqrt(v2)/sqrt(v_y)
                      gamma3 := g3*sqrt(v1)*sqrt(v2)/sqrt(v_y)"
# Model Fitting
fit.matchpair <- sem(data = dat.matchpair,
                     model = model.matchpair)
# Compute composite scores using first-order indicators
dat.centered <- dat.centered %>%
  mutate(
    PED.mean = rowMeans(select(., starts_with("PED")), na.rm = TRUE),
    SelfE.mean = rowMeans(select(., starts_with("SelfE")), na.rm = TRUE),
    PHQ.mean = rowMeans(select(., starts_with("PHQ")), na.rm = TRUE),
    PED.SelfE.mean = PED.mean*SelfE.mean - mean(PED.mean*SelfE.mean, na.rm = T)
  )
# Model Specification
model.rapi <- "# Measurement model
                 PHQ =~ 1*PHQ.mean
                 PED =~ 1*PED.mean
                 SelfE =~ 1*SelfE.mean
                 PED.SelfE =~ 1*PED.SelfE.mean
               # Error variance
                 PED.mean ~~ ev1*PED.mean
                 SelfE.mean ~~ ev2*SelfE.mean
                 PED.SelfE.mean ~~ ev3*PED.SelfE.mean
               # Latent variance
                 PED ~~ v1*PED
                 SelfE ~~ v2*SelfE
                 PED.SelfE ~~ v3*PED.SelfE
               # Error Constraints
                 ev1 == (1 - 0.8965932) * v1 / 0.8965932
                 ev2 == (1 - 0.8792078) * v2 / 0.8792078
                 ev3 == ev1 * v2 + ev2 * v1 + ev1 * ev2
               # Latent covariance
                 PED ~~ v12*SelfE
                 PED ~~ v13*PED.SelfE
                 SelfE ~~ v23*PED.SelfE
               # Residual variance of DV
                 PHQ ~~ v4*PHQ
               # Structural model
                 PHQ ~ g1*PED + g2*SelfE + g3*PED.SelfE
               # Standardized
                 v_y := g1^2*v1 + g2^2*v2 + g3^2*v3 + 2*g1*g2*v12 +
                        2*g1*g3*v13 + 2*g2*g3*v23 + v4
                 gamma1 := g1*sqrt(v1)/sqrt(v_y)
                 gamma2 := g2*sqrt(v2)/sqrt(v_y)
                 gamma3 := g3*sqrt(v1)*sqrt(v2)/sqrt(v_y)"
# Model Fitting
fit.rapi <- sem(data = dat.centered,
                model = model.rapi)
# Compute factor scores
model.fs <- "PHQ =~ PHQ1 + PHQ2 + PHQ3 + PHQ4 + PHQ5 + PHQ6 + PHQ7 + PHQ8 + PHQ9
             PED =~ PED1 + PED2 + PED3 + PED4 + PED5 + PED6 + PED7
             SelfE =~ SelfE1 + SelfE2 + SelfE3 + SelfE4 + SelfE5 +
                      SelfE6 + SelfE7 + SelfE8 + SelfE9 + SelfE10"
dat.fs <- get_fs(dat.centered,
                 model = model.fs,
                 method = "Bartlett",
                 std.lv = TRUE)
# obtain the single indicators
dat.fs <- dat.fs[ ,1:6]
colnames(dat.fs) <- gsub("_", ".", colnames(dat.fs))
# Obtain the factor scores as single indicators
dat.fs$fs.PED.SelfE <- dat.fs$fs.PED*dat.fs$fs.SelfE
dat.fs$fs.PED.SelfE <- dat.fs$fs.PED.SelfE - mean(dat.fs$fs.PED.SelfE)
# Compute the standard error of interaction
dat.fs$fs.PED.SelfE.se <- sqrt(1*dat.fs$fs.PED.se[1]^2 + 1*dat.fs$fs.SelfE.se[1]^2 +
                                 dat.fs$fs.PED.se[1]^2*dat.fs$fs.SelfE.se[1]^2)
# Model Specification
model.2spaint <- "# Measurement model
                    PHQ =~ 1*fs.PHQ
                    PED =~ 1*fs.PED
                    SelfE =~ 1*fs.SelfE
                    PED.SelfE =~ 1*fs.PED.SelfE
                  # Error variance
                    fs.PED ~~ 0.09875111*fs.PED
                    fs.SelfE ~~ 0.3397634*fs.SelfE
                    fs.PED.SelfE ~~ 0.22559*fs.PED.SelfE
                  # Latent variance
                    PED ~~ v1*PED
                    SelfE ~~ v2*SelfE
                    PED.SelfE ~~ v3*PED.SelfE
                  # Latent covariance
                    PED ~~ v12*SelfE
                    PED ~~ v13*PED.SelfE
                    SelfE ~~ v23*PED.SelfE
                  # Residual variance of DV
                    PHQ ~~ v4*PHQ
                  # Structural model
                    PHQ ~ b1*PED + b2*SelfE + b3*PED.SelfE
                  # Standardized
                    v_y := b1^2*v1 + b2^2*v2 + b3^2*v3 + 2*b1*b2*v12 +
                           2*b1*b3*v13 + 2*b2*b3*v23 + v4
                    beta1 := b1*sqrt(v1)/sqrt(v_y)
                    beta2 := b2*sqrt(v2)/sqrt(v_y)
                    beta3 := b3*sqrt(v1)*sqrt(v2)/sqrt(v_y)"
# Model Fitting
fit.2spaint <- sem(data = dat.fs,
                   model = model.2spaint)

# Measurement Model Fit
model.mm <- "# Measurement model
               PHQ =~ PHQ1 + PHQ2 + PHQ3 + PHQ4 + PHQ5 + PHQ6 + PHQ7 + PHQ8 + PHQ9
               PED =~ PED1 + PED2 + PED3 + PED4 + PED5 + PED6 + PED7
               SelfE =~ SelfE1 + SelfE2 + SelfE3 + SelfE4 + SelfE5 +
                        SelfE6 + SelfE7 + SelfE8 + SelfE9 + SelfE10"
fit.mm <- cfa(model.mm,
              data = dat.centered)
fitmeasures(fit.mm)
