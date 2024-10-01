# Read Data
ML <- "/Users/jimmy_z/R Projects/2SPA-Int/Sim_Data/continuous_bound_09252024-results_marklai-HP-Z4-G4-Workstation"
MLR <- "/Users/jimmy_z/R Projects/2SPA-Int/Sim_Data/continuous_boundMLR_09252024-results_marklai-HP-Z4-G4-Workstation"

ML_list <- list.files(path = ML, pattern = "results-row-.*\\.rds", full.names = TRUE)
MLR_list <- list.files(path = MLR, pattern = "results-row-.*\\.rds", full.names = TRUE)

ML_results_list <- lapply(ML_list, function(file) {
  data <- readRDS(file)
})
MLR_results_list <- lapply(MLR_list, function(file) {
  data <- readRDS(file)
})

ML_data_list <- lapply(ML_list, function(file) {
  data <- readRDS(file)
  data[["results"]]  
})
MLR_data_list <- lapply(MLR_list, function(file) {
  data <- readRDS(file)
  data[["results"]]  
})

ML_condition_list <- do.call(rbind, lapply(ML_list, function(file) {
  data <- readRDS(file)
  data[["condition"]]  
}))
MLR_condition_list <- do.call(rbind, lapply(MLR_list, function(file) {
  data <- readRDS(file)
  data[["condition"]]  
}))

# Inspect Case:
# Condition 1: N = 100, Cor = 0, Rel = 0.7, gamma = 0
# Condition 34: N = 100, Cor = 0.3, Rel = 0.7, gamma = 0
# Condition 52: N = 100, Cor = 0.6, Rel = 0.7, gamma = 0
# Condition 21: N = 100, Cor = 0, Rel = 0.7, gamma = 0.3
# Condition 25: N = 100, Cor = 0.3, Rel = 0.7, gamma = 0.3
# Condition 28: N = 100, Cor = 0.6, Rel = 0.7, gamma = 0.3

# Examine SE
extract_columns <- function(data, column_names, suffix) {
  columns <- setNames(
    lapply(column_names, function(name) {
      col_name <- paste0(name, suffix)
      if (col_name %in% names(data)) {
        return(data[[col_name]])
      } else {
        warning(paste("Column", col_name, "not found in the data. Returning NA."))
        return(rep(NA, nrow(data)))
      }
    }),
    paste0(column_names, suffix)
  )
  data.frame(columns)
}
column_names <- c("mmr", "upi", "rapi", "tspa")
ML_est_list <- lapply(ML_data_list, extract_columns, column_names = column_names, suffix = ".est")
ML_se_list <- lapply(ML_data_list, extract_columns, column_names = column_names, suffix = ".se_std")
MLR_est_list <- lapply(MLR_data_list, extract_columns, column_names = column_names, suffix = ".est")
MLR_se_list <- lapply(MLR_data_list, extract_columns, column_names = column_names, suffix = ".se_std")

ML_se_median <- do.call(rbind, lapply(ML_se_list, function(se) apply(se, 2, median, na.rm = TRUE)))
MLR_se_median <- do.call(rbind, lapply(MLR_se_list, function(se) apply(se, 2, median, na.rm = TRUE)))

round(MLR_se_median - ML_se_median, 3)

# Density Plot
par(mfrow = c(2, 2))
for (i in 1:ncol(se)) {
  col_sd <- sd(est[[i]], na.rm = TRUE)
  plot(density(se[[i]], na.rm = TRUE), 
       main = paste("Density Plot of SE -", colnames(se)[i]), 
       xlab = "SE", 
       col = "blue", 
       lwd = 2)
  abline(v = col_sd, col = "red", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))

# Mean and Median
apply(se, 2, mean)
apply(se, 2, median)
apply(est, 2, sd)
apply(est, 2, mad)
# Boxplot
par(mfrow = c(2, 2))
for (i in 1:ncol(se)) {
  boxplot(se[[i]], 
          main = paste("Boxplot of SE -", colnames(se)[i]), 
          ylab = "SE", 
          col = "lightgreen")
}
par(mfrow = c(1, 1))

# Outliers
Q1 <- quantile(se, 0.25, na.rm = TRUE)
Q3 <- quantile(se, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1
lower_bound <- (Q1 - 1.5 * IQR)
upper_bound <- (Q3 + 1.5 * IQR)
outliers <- se[se < lower_bound | se > upper_bound]
percentage <- length(outliers) / sum(!is.na(se)) * 100




