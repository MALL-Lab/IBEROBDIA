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
