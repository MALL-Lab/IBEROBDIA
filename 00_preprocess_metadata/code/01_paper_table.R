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
# Uncomment and add your packages here
# library()

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables

# Inputs

# Outputs

# Arguments

# ==============================================================================
# 1. Section Title 1
# ==============================================================================
# Add your code for Section 1 here.

# ==============================================================================
# 2. Section Title 2
# ==============================================================================
# Add your code for Section 2 here.

# ==============================================================================
# 3. Section Title 3
# ==============================================================================
# Add your code for Section 3 here.

# ==============================================================================
# 4. Section Title 4
# ==============================================================================
# Add your code for Section 4 here.

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()
library(readxl)
library(dplyr)
library(tableone)
library(nortest)

my_data <- read_excel("extdata/datos_voluntarios.xlsx")
colnames(my_data) <- make.names(colnames(my_data))
colnames(my_data)
cols <- c("Edad..años.", "talla..cm.", "Peso..Kg.", "IMC", "Perimetro..cintura..cm.",
          "GLUCOSA.INICIAL", "GLUCOSA.2H","TAS","TAD","HbA1c", "INSULINA.BASAL.µUl.mL",
          "HOMA", "TRIGLICERIDOS..mg.dL.",  "COLESTEROL...HDL")

my_data <- 
  my_data %>%
  mutate(CLASIF.IMC = ifelse(CLASIF.IMC == "Normopeso", "Normopeso", "Obeso"))
table(my_data$CLASIF.IMC)

# Filtramos los normopesos
normopeso_df <- my_data %>%
  filter(CLASIF.IMC == "Normopeso") %>%
  select(all_of(cols)) %>%
  mutate(across(where(is.character), as.numeric))

# Filtramos los obesos
obeso_df <- my_data %>%
  filter(CLASIF.IMC == "Obeso") %>%
  select(all_of(cols)) %>%
  mutate(across(where(is.character), as.numeric))


# Write a function that generates the appropriate summary based on normality
summarize_data <- function(df, variables) {
  
  summary_list <- list()
  
  for (var in variables) {
    # Check for normality using Anderson-Darling test (from nortest package)
    normality_test <- ad.test(df[[var]])
    
    if (normality_test$p.value > 0.05) {
      # If it's normally distributed, return mean ± SD
      summary_list[[var]] <- sprintf("%.2f ± %.2f", mean(df[[var]], na.rm = TRUE), sd(df[[var]], na.rm = TRUE))
    } else {
      # If not normally distributed, return median (min-max)
      summary_list[[var]] <- sprintf("%.2f (%.2f - %.2f)", median(df[[var]], na.rm = TRUE),
                                     min(df[[var]], na.rm = TRUE), max(df[[var]], na.rm = TRUE))
    }
  }
  
  return(summary_list)
}

# Call the function on your dataframe
np_results <- summarize_data(normopeso_df, cols)
ob_results <- summarize_data(obeso_df, cols)

# Print the results
np_results
ob_results
