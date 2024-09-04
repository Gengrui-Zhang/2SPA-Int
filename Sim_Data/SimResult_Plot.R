library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(ggplot2)
library(rlang)

# Read Data
data_path <- here("Sim_Data", "Match_07012024.rds")
df <- readRDS(data_path)

# Helper Function
process_data <- function(data) {
  sim_results <- data %>%
    gather("var", "val", raw_bias.rapi_yint_est:convergence_rate.tspa_yint_est) %>%
    dplyr::select(-c(SIM_TIME:WARNINGS)) %>%
    separate(col = var, into = c("stats", "parmet"), sep = "\\.") %>%
    separate(col = parmet, into = c("method", "par", "result"),  sep = "_") %>%
    dplyr::select(-result) %>%
    spread(stats, val) %>%
    relocate(REPLICATIONS, .after = last_col()) %>%
    mutate(N_lab = as_factor(paste0("italic(N) == ", N)),
           gammaint_lab = as_factor(paste0("\\beta_{xm} == ", gamma_xm)),
           cor_xm_lab = as_factor(paste0("corr[XM] == ", cor_xm)),
           rel_lab = as_factor(paste0("rho == ", rel)))

  return(sim_results)
}

# Plot function
plot_results <- function(data, x_var, y_var, x_label, y_label, y_limits = NULL, color_var = "method", shape_var = "method", point_size = 3, add_lines = FALSE) {
  plot <- data %>%
    ggplot(aes(x = factor(!!sym(x_var)), y = !!sym(y_var), color = !!sym(color_var), shape = !!sym(shape_var))) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = point_size, alpha = 0.8) +
    facet_grid(cor_xm_lab ~ rel_lab, labeller = label_parsed) +
    scale_shape_manual(values = c("reg" = 17, "rapi" = 16, "upi" = 18, "tspa" = 15),
                       labels = c("reg" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int")) +
    scale_color_manual(values = c("reg" = "pink", "rapi" = "#E69F00", "upi" = "#F0E442", "tspa" = "#009E73"),
                       labels = c("reg" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int")) +
    labs(x = x_label, y = y_label, color = "Method", shape = "Method") +
    theme_minimal() +
    theme(legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "right")
  
  if (add_lines) {
    plot <- plot + 
      geom_line(aes(group = !!sym(shape_var), linetype = !!sym(shape_var)), 
                position = position_dodge(width = 0.75), 
                size = 0.7, alpha = 0.6) +
      scale_linetype_manual(name = "Line",  # Set the legend name to "Line"
                            values = c("reg" = "dotted", "rapi" = "solid", "upi" = "dashed", "tspa" = "twodash"),
                            labels = c("reg" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int"))
  }
  
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)
  }
  
  return(plot)
}

# Split the data
type1data <- df %>% filter(gamma_xm == 0)
powerdata <- df %>% filter(gamma_xm == 0.3)

type1_pd <- process_data(type1data)
power_pd <- process_data(powerdata)

# Write the processed data
write_csv(type1_pd, "Sim_Data/type1_07012024.csv")
write_csv(power_pd, "Sim_Data/power_07012024.csv")

# Plot the results
type1_plot <- read.csv("Sim_Data/type1_07012024.csv")
# type1_plot <- type1_plot %>% filter(method != "reg")
power_plot <- read.csv("Sim_Data/power_07012024.csv")
# power_plot <- power_plot %>% filter(method != "reg")

# Standard Bias
sd_bias_plot <- plot_results(power_plot,
                             "N",
                             "std_bias",
                             "Sample Size (N)",
                             "Standardized Bias",
                             c(-0.5, 0.5),
                             add_lines = TRUE)
sd_bias_plot <- sd_bias_plot +
  geom_hline(yintercept = c(-0.4, 0.4), linetype = "dashed", color = "red")

# Median-Mad Relative SE Bias
rse_bias_plot <- plot_results(power_plot,
                              "N",
                              "stdMed_rse_bias",
                              "Sample Size (N)",
                              "Relative SE Bias",
                              c(-0.3, 0.3))
rse_bias_plot <- rse_bias_plot +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "red")

# Coverage rate
cov_bias_plot <- plot_results(power_plot,
                              "N",
                              "coverage",
                              "Sample Size (N)",
                              "Coverage Rate of 95% CI")
cov_bias_plot <- cov_bias_plot +
  geom_hline(yintercept = c(0.91), linetype = "dashed", color = "red")

# RMSE
rmse_bias_plot <- plot_results(power_plot,
                               "N",
                               "rmse",
                               "Sample Size (N)",
                               "Root Mean Squre Error",
                               add_lines = TRUE)

# Type I error rate
type1_bias_plot <- plot_results(type1_plot,
                              "N",
                              "type1_lrt",
                              "Sample Size (N)",
                              "Type I Error Rate",
                              c(0, 0.2))

# Statistical Power
power_bias_plot <- plot_results(power_plot,
                                "N",
                                "power_lrt",
                                "Sample Size (N)",
                                "Statistical Power")
power_bias_plot <- power_bias_plot +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "red")
