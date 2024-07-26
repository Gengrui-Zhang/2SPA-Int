seed_value <- 66234
set.seed(seed_value)

library(semTools)
library(lavaan)
library(semPlot)
library(SimDesign)
library(MASS)
library(mnormt)
library(dplyr)
library(tidyverse)
# library(MBESS)
library(semlrtp)
# Directly call functions (call R2spa)
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/get_fscore.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/get_fsint.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/upi.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/rapi_new.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/tspa_old.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/parseInteractionTerms.R")
source("/Users/jimmy_z/R Projects/2S-PA-Int/R/parseIndicators.R")

# Data Generation

# Helper Function
generate_sem_data <- function(N, model, Alpha, Phi, Lambda, Gamma, Theta, SD_y) {
  # Generate scores for observed items: x1 - x3, m1 - m3
  eta_scores <- rmnorm(N, mean = Alpha, varcov = Phi) # Simulate latent scores
  # Generate eta scores with interaction
  eta_scores_int <- cbind(eta_scores, eta_scores[,1]*eta_scores[,2])
  delta <- rmnorm(N, mean = rep(0, length(diag(Theta))), varcov = Theta) # Errors/Residuals of Indicators
  item_scores <- tcrossprod(eta_scores, Lambda) + delta # Item/Indicator scores

  y_scores <- tcrossprod(eta_scores_int, Gamma) + rnorm(N, 0, SD_y) # Y scores / DV scores

  # Parsing formula
  indicator_terms <- unlist(strsplit(gsub(" ", "",
                                          unlist(strsplit(model, "\n"))[grep("=~",
                                                                             unlist(strsplit(model, "\n")))]),
                                     split = "=~"))
  indicator_vars <- indicator_terms[grepl("+", indicator_terms, fixed = TRUE) == "FALSE"]
  indicator_num <- as.vector(unlist(lapply(strsplit(indicator_terms[grepl("+", indicator_terms, fixed = TRUE) == "TRUE"],
                                                    split = "+", fixed = TRUE), length)))

  df <- as.data.frame(cbind(item_scores, y_scores))
  names <- list()
  for (n in seq(length(indicator_vars))) {
    names[[n]] <- tolower(paste0(indicator_vars[n], seq(indicator_num[n])))
  }
  colnames(df) <- c(unlist(names), "Y")
  return(df)
}

DESIGNFACTOR <- createDesign(
  N = c(100, 250, 500),
  cor_xm = c(0, 0.3, 0.6), # correlation between latent x and m / error variance of Y
  rel = c(0.7, 0.8, 0.9),
  gamma_xm = c(0, 0.3) # Two levels of the interaction effect
)

FIXED_PARAMETER <- list(model = '
                                  # Measurement Model
                                    X =~ x1 + x2 + x3
                                    M =~ m1 + m2 + m3
                                  # Structural Model
                                    Y ~ b1*X + b2*M + b3*X:M
                                  # Define Standardized Coefficients
                                    X ~~ v1*X
                                    M ~~ v2*M
                                    beta1 := b1*sqrt(v1)
                                    beta2 := b2*sqrt(v2)
                                    beta3 := b3*sqrt(v1)*sqrt(v2)
                                  ',
                        beta1 = 1,  # fixed
                        beta2 = 0.9,  # fixed
                        beta3 = 0.75,  # three conditions
                        mu_x = 0, # latent mean of x: fixed at 0
                        mu_m = 0, # latent mean of m: fixed at 0,
                        gamma_x = 0.3,
                        gamma_m = 0.3
)


GenData <- function (condition, fixed_objects = NULL) {
  N <- condition$N # Sample size
  mu_x <- fixed_objects$mu_x
  mu_m <- fixed_objects$mu_m
  beta1 <- fixed_objects$beta1 # beta 1: fixed at 1
  beta2 <- fixed_objects$beta2 # beta 2: fixed at 0.9
  beta3 <- fixed_objects$beta3 # beta 3: fixed at 0.75
  cor_xm <- condition$cor_xm # latent correlation: varied
  gamma_x <- fixed_objects$gamma_x
  gamma_m <- fixed_objects$gamma_m
  gamma_xm <- condition$gamma_xm
  rel = condition$rel

  # Compute disturbance variance
  sd_y <- sqrt(1 - (gamma_x^2 + gamma_m^2 + gamma_xm^2 + 2*gamma_x*gamma_m*cor_xm))

  # Compute error variance
  sum_error <- sum(c(beta1, beta2, beta3))^2*(1 - rel)/rel
  err_var <- sum_error*c(0.44, 0.33, 0.23)

  # Simulate SEM data
  Alpha <- c(mu_x, mu_m) # Latent means
  Phi <- matrix(c(1, cor_xm,
                  cor_xm, 1), nrow = 2) # latent var/cov
  Lambda <- cbind(c(beta1, beta2, beta3, rep(0, 3)),
                  c(rep(0, 3), beta1, beta2, beta3)) # factor loadings
  Theta <- diag(err_var,
                nrow = 6)
  Gamma <- matrix(c(gamma_x, gamma_m, gamma_xm), nrow = 1)
  SD_y <- sd_y

  generate_sem_data(N,
                    model = fixed_objects$model,
                    Alpha = Alpha,
                    Phi = Phi,
                    Lambda = Lambda,
                    Theta = Theta,
                    Gamma = Gamma,
                    SD_y = SD_y
  )
}

extract_res <- function (condition, dat, fixed_objects = NULL) {

  # Fit using moderated multiple regression
  dat_reg <- dat %>%
    mutate(x_c = rowSums(dat[c("x1", "x2", "x3")]) - mean(rowSums(dat[c("x1", "x2", "x3")])),
           m_c = rowSums(dat[c("m1", "m2", "m3")]) - mean(rowSums(dat[c("m1", "m2", "m3")])),
           xm = x_c*m_c)
  fit_reg <- sem(model = "Y ~ c0*x_c + c1*m_c + c2*xm",
                 data = dat_reg)
  if (lavInspect(fit_reg, what = "converged")) {
    std_col <- standardizedSolution(fit_reg)
    reg_est <- std_col[std_col$label == "c2", "est.std"]
    reg_se <- std_col[std_col$label == "c2", "se"]
    reg_lrtp <- lrtp(fit_reg)
    reg_lrtp_ci <- as.data.frame(reg_lrtp) %>%
      filter(op == "~" & label == "c2") %>%
      select(ci.lower, ci.upper)
  } else {
    reg_est <- NA
    reg_se <- NA
  }
  # Fit using rapi function
  fit_rapi <- rapi(model = fixed_objects$model,
                   data = dat)
  if (lavInspect(fit_rapi, what = "converged")) {
    rapi_est <- coef(fit_rapi, type = "user")["beta3"]
    rapi_se <- sqrt(vcov(fit_rapi, type = "user")["beta3", "beta3"])
    rapi_lrtp <- lrtp(fit_rapi)
    rapi_lrtp_ci <- as.data.frame(rapi_lrtp) %>%
      filter(op == "~" & label == "b3") %>%
      select(ci.lower, ci.upper)
  } else {
    rapi_est <- NA
    rapi_se <- NA
  }
  # Fit using upi function
  fit_upi <- upi(model = fixed_objects$model,
                 data = dat,
                 mode = "match")
  if (lavInspect(fit_upi, what = "converged")) {
    upi_est <- coef(fit_upi, type = "user")["beta3"]
    upi_se <- sqrt(vcov(fit_upi, type = "user")["beta3", "beta3"])
    upi_lrtp <- lrtp(fit_upi)
    upi_lrtp_ci <- as.data.frame(upi_lrtp) %>%
      filter(op == "~" & label == "b3") %>%
      select(ci.lower, ci.upper)
  } else {
    upi_est <- NA
    upi_se <- NA
  }
  # Fit using tspa function
  fs_dat <- get_fs(dat,
                   model ='
                             X =~ x1 + x2 + x3
                             M =~ m1 + m2 + m3
                             ',
                   method = "Bartlett",
                   std.lv = TRUE)
  Y <- dat$Y
  fs_dat <- cbind(fs_dat, Y)
  fit_tspa <- tspa(model = "Y ~ b1*X + b2*M + b3*X:M
                              beta1 := b1 * sqrt(v1)
                              beta2 := b2 * sqrt(v2)
                              beta3 := b3 * sqrt(v1) * sqrt(v2)",
                   data = fs_dat,
                   se = list(X = fs_dat$fs_X_se[1],
                             M = fs_dat$fs_M_se[1]))
  if (lavInspect(fit_tspa, what = "converged")) {
    tspa_est <- coef(fit_tspa, type = "user")["beta3"]
    tspa_se <- sqrt(vcov(fit_tspa, type = "user")["beta3", "beta3"])
    tspa_lrtp <- lrtp(fit_tspa)
    tspa_lrtp_ci <- as.data.frame(tspa_lrtp) %>%
      filter(op == "~" & label == "b3") %>%
      select(ci.lower, ci.upper)
  } else {
    tspa_est <- NA
    tspa_se <- NA
  }
  # Extract parameter estimates and standard errors
  paret <- c(reg_est, reg_se, reg_lrtp_ci$ci.lower, reg_lrtp_ci$ci.upper,
             rapi_est, rapi_se, rapi_lrtp_ci$ci.lower, rapi_lrtp_ci$ci.upper,
             upi_est, upi_se, upi_lrtp_ci$ci.lower, upi_lrtp_ci$ci.upper,
             tspa_est, tspa_se, tspa_lrtp_ci$ci.lower, tspa_lrtp_ci$ci.upper)
  names(paret) <- c("reg_yint_est", "reg_yint_se", "reg_lrtp_ci_lower", "reg_lrtp_ci_upper",
                    "rapi_yint_est", "rapi_yint_se", "rapi_lrtp_ci_lower", "rapi_lrtp_ci_upper",
                    "upi_yint_est", "upi_yint_se", "upi_lrtp_ci_lower", "upi_lrtp_ci_upper",
                    "tspa_yint_est", "tspa_yint_se", "tspa_lrtp_ci_lower", "tspa_lrtp_ci_upper")
  return(paret)
}

evaluate_res <- function (condition, results, fixed_objects = NULL) {

  # Population parameter
  pop_par <- condition$gamma_xm

  # Separate estimates and se
  results_est <- as.data.frame(results[colnames(results)[grepl("_est", colnames(results))]])
  results_se <- as.data.frame(results[colnames(results)[grepl("_se", colnames(results))]])
  lrtp_ci_lower <- as.data.frame(results[colnames(results)[grepl("_ci_lower", colnames(results))]])
  lrtp_ci_upper <- as.data.frame(results[colnames(results)[grepl("_ci_upper", colnames(results))]])

  # Descriptive of SEs
  descriptive <- function(est, se, type = NULL) {
    output <- numeric(ncol(est))
    if (type == "SD") {
      output <- apply(est, 2L, sd, na.rm = T)
    } else if (type == "MeanSE") {
      output <- apply(se, 2, mean, na.rm = T)
    } else if (type == "MedianSE") {
      output <- apply(se, 2, median, na.rm = TRUE)
    } else if (type == "MAD") {
      output <- apply(est, 2, function(x) mad(x, na.rm = TRUE))
    }
    return(output)
  }

  # Helper function: robust bias
  robust_bias <- function(est, se, pop_par, trim = 0, type = NULL) {
    output <- numeric(ncol(est))
    for (i in seq_len(ncol(est))) {
      if (type == "raw") {
        output[i] <- mean((est[,i] - pop_par), na.rm = TRUE)
      } else if (type == "standardized") {
        output[i] <- (mean(est[,i], na.rm = TRUE) - pop_par)/sd(est[,i], na.rm = TRUE)
      } else if (type == "trim") {
        output[i] <- mean(est[,i], trim = trim, na.rm = TRUE) - pop_par
      } else if (type == "median") {
        output[i] <- (median(est[,i], na.rm = TRUE) - pop_par) / mad(est[,i], na.rm = TRUE)
      } else {
        output[i] <- (mean(est[,i], trim = trim, na.rm = TRUE) - pop_par) / sd(est[,i], na.rm = TRUE)
      }
    }
    names(output) <- colnames(est)
    return(output)
  }

  # Helper function: relative SE bias
  rse_bias <- function(est, est_se, trim = 0, type = "raw") {
    if (type == "raw") {
      est_se <- as.matrix(est_se)
      est <- as.matrix(est)
      est_se_mean <- apply(est_se, 2, mean, na.rm = T)
      emp_sd <- apply(est, 2L, sd, na.rm = T)
      rse_bias <- est_se_mean / emp_sd - 1
    } else if (type == "median") {
      est_se <- as.matrix(est_se)
      est <- as.matrix(est)
      est_se_median <- apply(est_se, 2, median, na.rm = TRUE)
      emp_mad <- apply(est, 2, function(x) mad(x, na.rm = TRUE))
      rse_bias <- est_se_median / emp_mad - 1
    } else if (type == "trim") {
      est_se <- as.matrix(est_se)
      est <- as.matrix(est)
      est_se_mean <- apply(est_se, 2, mean, trim = trim, na.rm = TRUE)
      emp_sd <- apply(est, 2L, sd, na.rm = T)
      rse_bias <- est_se_mean / emp_sd - 1
    }
    return(rse_bias)
  }

  # Helper function: detecting outliers for SE
  outlier_se <- function(est_se) {
    results <- c()
    for(column in names(est_se)) {
      # Calculate Q1, Q3, and IQR
      Q1 <- quantile(est_se[[column]], 0.25, na.rm = TRUE)
      Q3 <- quantile(est_se[[column]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      # Determine outliers
      lower_bound <- (Q1 - 1.5 * IQR)
      upper_bound <- (Q3 + 1.5 * IQR)
      outliers <- est_se[[column]][est_se[[column]] < lower_bound | est_se[[column]] > upper_bound]
      # Calculate the percentage of outliers
      percentage <- length(outliers) / sum(!is.na(est_se[[column]])) * 100
      results[column] <- percentage
    }
    return(results)
  }

  # Helper function for calculating coverage rate, Type I error rate, and power
  ci_stats <- function(est, est_se, par, stats_type, lrt_lo_95 = NULL, lrt_hi_95 = NULL) {
    est_se <- as.matrix(est_se)
    est <- as.matrix(est)

    # Calculate the confidence intervals
    lo.95 <- est - qnorm(.975) * est_se
    hi.95 <- est + qnorm(.975) * est_se
    ci_est <- vector("list", length = ncol(est))
    names(ci_est) <- colnames(est)

    # Construct confidence intervals for each method
    for (i in seq_len(ncol(est))) {
      ci_est[[i]] <- cbind(lo.95[,i], hi.95[,i])
    }

    # Extract LRT CIs
    if (!is.null(lrt_lo_95) && !is.null(lrt_hi_95)) {
      lrt_lo_95 <- as.matrix(lrt_lo_95)
      lrt_hi_95 <- as.matrix(lrt_hi_95)
      ci_lrt <- vector("list", length = ncol(est))
      names(ci_lrt) <- colnames(est)

      for (i in seq_len(ncol(est))) {
        ci_lrt[[i]] <- cbind(lrt_lo_95[,i], lrt_hi_95[,i])
      }
    }

    # Determine which statistic to calculate
    if (stats_type == "Coverage") {
      return(sapply(ci_est, function(ci) mean(ci[,1] <= par & ci[,2] >= par)))
    } else if (stats_type == "TypeI") {
      return(sapply(ci_est, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
    } else if (stats_type == "Lrt_TypeI") {
      return(sapply(ci_lrt, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
    } else if (stats_type == "Power") {
      return(sapply(ci_est, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
    } else if (stats_type == "Lrt_Power") {
      return(sapply(ci_lrt, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
    } else {
      return("Invalid stats_type specified. Please choose from 'Coverage', 'TypeI', or 'Power'.")
    }
  }

  # Helper function for convergence rate
  convergence_rate <- function(est) {
    apply(est, 2, function(x) 1-(sum(is.na(x)) / length(x)))
  }

  c(raw_bias = robust_bias(results_est,
                           results_se,
                           pop_par,
                           type = "raw"),
    std_bias = robust_bias(results_est,
                           results_se,
                           pop_par,
                           type = "standardized"),
    trim_bias = robust_bias(results_est,
                            results_se,
                            pop_par,
                            trim = 0.2,
                            type = "trim"), # 20% trimmed mean
    stdMed_bias = robust_bias(results_est,
                              results_se,
                              pop_par,
                              type = "median"),
    coverage = ci_stats(results_est, results_se, pop_par, "Coverage"),
    type1 = ci_stats(results_est, results_se, pop_par, "TypeI"),
    type1_lrt = ci_stats(results_est, results_se, pop_par, "Lrt_TypeI",
                         lrt_lo_95 = lrtp_ci_lower,
                         lrt_hi_95 = lrtp_ci_upper),
    power = ci_stats(results_est, results_se, pop_par, "Power"),
    power_lrt = ci_stats(results_est, results_se, pop_par, "Lrt_Power",
                         lrt_lo_95 = lrtp_ci_lower,
                         lrt_hi_95 = lrtp_ci_upper),
    rmse = RMSE(na.omit(results_est),
                parameter = pop_par),
    raw_rse_bias = rse_bias(results_est,
                            results_se,
                            type = "raw"),
    stdMed_rse_bias = rse_bias(results_est,
                               results_se,
                               type = "median"),
    SD = descriptive(results_est,
                     results_se,
                     type = "SD"),
    MeanSE = descriptive(results_est,
                         results_se,
                         type = "MeanSE"),
    MedianSE = descriptive(results_est,
                           results_se,
                           type = "MedianSE"),
    MAD = descriptive(results_est,
                      results_se,
                      type = "MAD"),
    trim_rse_bias = rse_bias(results_est,
                             results_se,
                             trim = 0.2,
                             type = "trim"),
    outlier_se = outlier_se(results_se),
    convergence_rate = convergence_rate(results_est)
  )
}

# Run 2000 replications

Match_07012024 <- runSimulation(design = DESIGNFACTOR,
                                replications = 2000,
                                generate = GenData,
                                analyse = extract_res,
                                summarise = evaluate_res,
                                fixed_objects = FIXED_PARAMETER,
                                save = TRUE,
                                save_results = TRUE,
                                filename = "Match_07012024",
                                control = list(allow_na = TRUE),
                                parallel = TRUE,
                                ncores = min(4L, parallel::detectCores() - 1))

Match_07012024$seed_value <- seed_value
