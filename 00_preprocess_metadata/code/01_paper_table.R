# ==============================================================================
# Script Title: 01_paper_table.R
# Description: This script generates a table for the paper.
# ==============================================================================
# Author/s: @DiegoFE94 (GitHub)
# Affiliation: Machine Learning for Life Sciences Laboratory (MALL)
# Email: diego.fedreira@udc.es
# Date Created: 2024 (Last Update on January 2025)
# ==============================================================================
# 0. Set Up
# ==============================================================================
# 0.1 Clean Environment
rm(list = ls())
set.seed(1965)

# 0.2 Load Packages
library(readxl)
library(dplyr)
library(tableone)
library(nortest)

# 0.3 Declare or Load Custom Functions
# Write a function that generates the appropriate summary based on normality
summarize_data <- function(df, variables) {
  require(nortest)
  summary_list <- list()
  for (var in variables) {
    # Check for normality using Anderson-Darling test (from nortest package)
    normality_test <- ad.test(df[[var]])
    if (normality_test$p.value > 0.05) {
      # If it's normally distributed, return mean ± SD
      summary_list[[var]] <- sprintf("%.2f ± %.2f",
                                     mean(df[[var]], na.rm = TRUE),
                                     sd(df[[var]], na.rm = TRUE))
    } else {
      # If not normally distributed, return median (min-max)
      summary_list[[var]] <- sprintf("%.2f (%.2f - %.2f)",
                                     median(df[[var]], na.rm = TRUE),
                                     min(df[[var]], na.rm = TRUE),
                                     max(df[[var]], na.rm = TRUE))
    }
  }
  return(summary_list)
}

# 0.4 Declare Variables
# Inputs
input_path <- "extdata/datos_voluntarios.xlsx"
# Outputs

# Arguments

# ==============================================================================
# 1. Load data
# ==============================================================================
# 1.1 Load the data
my_data <- read_excel(input_path)
colnames(my_data) <- make.names(colnames(my_data))
colnames(my_data)

# 1.2 Prepare data
cols <- c("Edad..años.", "talla..cm.", "Peso..Kg.", "IMC", "Perimetro..cintura..cm.",
          "GLUCOSA.INICIAL", "GLUCOSA.2H","TAS","TAD","HbA1c", "INSULINA.BASAL.µUl.mL",
          "HOMA", "TRIGLICERIDOS..mg.dL.",  "COLESTEROL...HDL")
my_data <- my_data %>%
  mutate(CLASIF.IMC = ifelse(CLASIF.IMC == "Normopeso", "Normopeso", "Obeso"))
table(my_data$CLASIF.IMC)

# 1.3 Normal weight filter
normopeso_df <- my_data %>%
  filter(CLASIF.IMC == "Normopeso") %>%
  select(all_of(cols)) %>%
  mutate(across(where(is.character), as.numeric))
# Call the function on your dataframe
np_results <- summarize_data(normopeso_df, cols)

# 1.4 Obesity filter
obeso_df <- my_data %>%
  filter(CLASIF.IMC == "Obeso") %>%
  select(all_of(cols)) %>%
  mutate(across(where(is.character), as.numeric))

# Call the function on your dataframe
ob_results <- summarize_data(obeso_df, cols)

# 1.5 Print the results
np_results
ob_results

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()