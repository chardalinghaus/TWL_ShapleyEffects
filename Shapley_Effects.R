rm(list = ls())
while (!is.null(dev.list())) dev.off()

# Set the working directory
setwd("C:/")
file <- "TWL_MCsimul_Fut1"

# Load necessary packages
library(sensitivity)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(gtools)
library(boot)
library(RANN)
library(corpcor)
library(whitening)
library(TSP)
library(Cairo)

#Calculate Sobol Indices and Shapley effects

# Read the future data
twl_samples <- read.table(paste0(file, ".txt"), header = TRUE)

# Convert the "Date" column to a proper date format
twl_samples$Date <- as.Date(twl_samples$Date, format = "%Y-%m-%d")

# Extract the year from the "Date" column
twl_samples$Year <- format(twl_samples$Date, "%Y")

# Separate rows per year into a list of data frames
twl_per_year <- split(twl_samples, twl_samples$Year)

# Create a directory to save the results
output_directory <- file  # Change this to the desired directory path
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

#Create a model - TWL
TWL <- function(X) {
  return(rowSums(X[, c("SLR", "AT", "SS", "Runup")]))
}

# Iterate through each year and compute Sobol Indices and Shapley effects
for (year in names(twl_per_year)) {
  year_twl <- twl_per_year[[year]]
  # Compute Sobol and Shapley effects using KNN method
  X <- year_twl[, c("SLR", "AT", "SS", "Runup", "Setup_Eq", "Swash_Eq", "βf", "Scenario", "Model")]
  Y <- year_twl[, c("TWL")]
  
  results <- shapleysobol_knn(model=TWL, X, method = "knn", n.knn = 20, n.limit = 2000, U = NULL,
                              n.perm = NULL, noise = F, rescale = T, parl=10)
  #If U = 0, calculates total Sobol indices. If U = 1, calculates first-order Sobol indices.

  # Calculate Shapley effects
  shapley_effects <- tell(results,Y)
  print(shapley_effects)
  
  #Plot and save Shapley
  plot_shap <- file.path(output_directory, paste0("shapley_plot_", year, ".png"))
  CairoPNG(filename = plot_shap, width = 800, height = 600)
  plot(results)
  dev.off()
  
  # Calculate Sobol indices
  sobol_indices <- extract(results)
  print(sobol_indices)
  
  # Save Shapley effects and Sobol indices as CSV files
  shapley_values <- shapley_effects$Shap
  shapley_effects_file <- file.path(output_directory, paste0("shapley_effects_", year, ".csv"))
  write.csv(shapley_values, shapley_effects_file, row.names = FALSE)
  sobol_indices_file <- file.path(output_directory, paste0("sobol_indices_", year, ".csv"))
  write.csv(sobol_indices, sobol_indices_file, row.names = FALSE)
  
}