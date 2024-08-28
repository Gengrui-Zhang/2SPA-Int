library(MplusAutomation)
library(dplyr)
library(parallel)
library(pbapply)

# Define the directory paths
output_dir <- "/Users/jimmy_z/R Projects/2SPA-Int/Sim_Script/tspa_int_continous_mplus/"
data_dir <- "/Users/jimmy_z/desktop/Temp_Data/"
NREPS <- 500

# Conditions
observations <- c(100, 250, 500)
cor_xm <- c(0, 0.3, 0.6)
int_path <- c("typeO", "power")
reliabilities <- c(0.7, 0.8, 0.9)

# Function to calculate residual variance for y
calc_resid_y <- function(gamma_xm, cor_xm) {
  return(1 - (0.3^2 + 0.3^2 + gamma_xm^2 + 2 * 0.3 * 0.3 * cor_xm))
}

# Function to clean up temporary files
cleanup_files <- function(files) {
  lapply(files, function(f) {
    if (file.exists(f)) file.remove(f)
  })
}

# Function to run a single replication
run_mplus_simulation <- function(obs, cor, path, rel, rep, output_dir, data_dir) {
  sum_error <- (1 + 0.9 + 0.75)^2 * (1 - rel) / rel
  err_var <- c(0.44, 0.33, 0.23) * sum_error
  gamma_xm <- ifelse(path == "typeO", 0, 0.3)
  resid_y <- calc_resid_y(gamma_xm, cor)
  
  temp_data_file <- paste0(data_dir, "temp_data_", rep, ".dat")
  
  # Generation script with saving temporary data
  gen_script <- paste0(
    "TITLE: Data Generation for SEM with Latent Interaction\n\n",
    "MONTECARLO:\n",
    "  NAMES ARE x1 x2 x3 m1 m2 m3 y;\n",
    "  NOBSERVATIONS = ", obs, ";\n",
    "  NREPS = 1;\n",
    "  SEED = ", sample(1:1000000, 1), ";\n",
    "  SAVE = ", temp_data_file, ";\n\n",
    "MODEL POPULATION:\n",
    "  eta1 BY x1*1 x2*0.9 x3*0.75;\n",
    "  eta2 BY m1*1 m2*0.9 m3*0.75;\n",
    "  [eta1@0 eta2@0];\n",
    "  eta1*1 eta2*1;\n",
    "  eta3 | eta1 XWITH eta2;\n",
    "  eta1 WITH eta2*", cor, ";\n",
    "  y ON eta1*0.3 eta2*0.3 eta3*", gamma_xm, ";\n",
    "  x1*", err_var[1], " x2*", err_var[2], " x3*", err_var[3], ";\n",
    "  m1*", err_var[1], " m2*", err_var[2], " m3*", err_var[3], ";\n",
    "  y*", resid_y, ";\n\n",
    "ANALYSIS:\n",
    "  TYPE = RANDOM;\n",
    "  ESTIMATOR = MLR;\n",
    "  ALGORITHM = INTEGRATION;\n",
    "  INTEGRATION = MONTECARLO;\n\n",
    "OUTPUT:\n  TECH9;\n\n"
  )
  
  gen_filename <- paste0(output_dir, "temp_gen_", rep, ".inp")
  writeLines(gen_script, gen_filename)
  model_run <- runModels(gen_filename)
  
  if (model_run) {
    # Analysis script using the temporary data file
    analysis_script <- paste0(
      "TITLE: Data Analysis for SEM with Latent Interaction\n\n",
      "DATA:\n  FILE = ", temp_data_file, ";\n  FORMAT = FREE;\n\n",
      "VARIABLE:\n  NAMES ARE x1 x2 x3 m1 m2 m3 y;\n  USEVARIABLES ARE x1 x2 x3 m1 m2 m3 y;\n\n",
      "ANALYSIS:\n",
      "  TYPE = RANDOM;\n",
      "  ESTIMATOR = MLR;\n",
      "  ALGORITHM = INTEGRATION;\n",
      "  INTEGRATION = MONTECARLO;\n\n",
      "MODEL:\n",
      "  eta1 BY x1 x2 x3;\n",
      "  eta2 BY m1 m2 m3;\n",
      "  eta3 | eta1 XWITH eta2;\n",
      "  y ON eta1 eta2 eta3;\n\n",
      "OUTPUT:\n  STANDARDIZED;\n  TECH1;\n  TECH4;\n\n"
    )
    
    analysis_filename <- paste0(output_dir, "temp_analysis_", rep, ".inp")
    writeLines(analysis_script, analysis_filename)
    model_run <- runModels(analysis_filename)
    
    if (model_run) {
      analysis_output <- paste0(output_dir, "temp_analysis_", rep, ".out")
      par_est <- readModels(analysis_output, what = "parameters")$parameters
      
      std_coef <- par_est$stdyx.standardized
      usd_coef <- par_est$unstandardized
      
      # Check if results exist
      if (!is.null(std_coef) && !is.null(usd_coef)) {
        est_std <- std_coef[std_coef$paramHeader == "Y.ON" & std_coef$param == "ETA3", ]$est
        se_std <- std_coef[std_coef$paramHeader == "Y.ON" & std_coef$param == "ETA3", ]$se
        se_usd <- usd_coef[usd_coef$paramHeader == "Y.ON" & usd_coef$param == "ETA3", ]$se
        
        # Clean up temporary files
        cleanup_files(c(gen_filename, analysis_filename, temp_data_file))
        
        return(data.frame(obs, cor, path, rel, rep, est_std, se_std, se_usd))
      }
    }
  }
  
  # If the model run fails, return NA for the parameters but still include the replication
  cleanup_files(c(gen_filename, analysis_filename, temp_data_file))
  return(data.frame(obs, cor, path, rel, rep, est_std = NA, se_std = NA, se_usd = NA))
}

# Split the replications into batches
batch_size <- 50  # Adjust based on your system
rep_batches <- split(1:NREPS, ceiling(seq_along(1:NREPS) / batch_size))

# Process each batch sequentially with progress tracking
final_results <- list()

for (obs in observations) {
  for (cor in cor_xm) {
    for (path in int_path) {
      for (rel in reliabilities) {
        cat("Starting batch for conditions: obs =", obs, ", cor =", cor, ", path =", path, ", rel =", rel, "\n")
        
        for (batch_num in seq_along(rep_batches)) {
          cat("Starting batch", batch_num, "of", length(rep_batches), "\n")
          
          batch_results <- pblapply(rep_batches[[batch_num]], function(rep) {
            run_mplus_simulation(obs = obs, cor = cor, path = path, rel = rel, rep = rep, output_dir = output_dir, data_dir = data_dir)
          }, cl = detectCores() - 1)  # Leaving one core free for system processes
          
          batch_results <- batch_results[!sapply(batch_results, is.null)]  # Remove NULL results
          
          if (length(batch_results) > 0) {
            if (!exists("results_df")) {
              results_df <- do.call(rbind, batch_results)
            } else {
              results_df <- rbind(results_df, do.call(rbind, batch_results))
            }
          }
          
          cat("Finished batch", batch_num, "of", length(rep_batches), "\n")
        }
        
        cat("Finished all batches for conditions: obs =", obs, ", cor =", cor, ", path =", path, ", rel =", rel, "\n")
      }
    }
  }
}

# Save the accumulated results to a CSV file if any results were obtained
if (exists("results_df") && nrow(results_df) > 0) {
  write.csv(results_df, file = paste0(output_dir, "LMS_Standardized_Results.csv"), row.names = FALSE)
  cat("Results saved to LMS_Standardized_Results.csv\n")
} else {
  cat("No valid results to save.\n")
}

# Clean up
if (exists("results_df")) {
  rm(results_df)
}
