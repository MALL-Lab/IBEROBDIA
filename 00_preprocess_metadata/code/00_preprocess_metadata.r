# Clinical Data
## Clinical Data Extraction
require(tibble)
require(readxl)
require(dplyr)
path = "extdata/FASTQ/"
names = list.files(path = path, pattern = ".fastq.gz")
library(stringr)
library(data.table)

# Clinical Data
vars = list()
for (i in seq_along(names)) {
  vars[i] = str_remove(names[i], ".fastq.gz")
  vars[i] = str_split(vars[i], "_", n = Inf, simplify = FALSE)
}
metalist = list()
for (i in seq_along(vars)) {
  if (length(vars[[i]]) == 3){
    metalist[[i]] = c(vars[[i]][1],vars[[i]][2],"Healthy",vars[[i]][3])
  }else{
    metalist[[i]] = vars[[i]]
  }
}
clinical = as.data.frame(do.call(rbind, metalist))
colnames(clinical)= c("Sample_ID","CI", "DT2_All","Gender")
clinical <-  as_tibble(x = clinical, rownames = NA)

# Mutate to cerate new Variables
clinical <- clinical %>% 
  select(Sample_ID, CI, Gender, DT2_All) %>%
  mutate( DT2_P_H = ifelse(DT2_All == "GAA", "PreDT2",
                              ifelse(DT2_All == "TGA", "PreDT2",
                                     DT2_All))) %>%
  mutate( Aff_H = ifelse(DT2_All == "GAA", "Affected",
                           ifelse(DT2_All == "TGA", "Affected",
                                  ifelse(DT2_All == "DT2", "Affected",
                                  DT2_All))))

# Load New Clinical Data
clinical_v2 <- readxl::read_xlsx(path = "extdata/clinical_data.xlsx")
clinical_v2 <- as_tibble(clinical_v2)
clinical_v2 <- clinical_v2 %>%
  select(Sample_ID, MS, GS, BMI, IR)

full_df <- dplyr::full_join(clinical, clinical_v2, by = "Sample_ID")

metadata <- full_df %>%
  # Delete 084CMLES sample (No FASTAQ)
  filter(Sample_ID != "084CMLES") %>%
  # Change Metabolic Syndrome variable
  mutate(MS = ifelse(MS == "No", "No_MS", MS)) %>% # No Metabolic Syndrome
  mutate(MS = ifelse(MS == "Sí", "MS", MS)) %>% # Metabolic Syndrome
  # Change Status of Glucosevariable
  mutate(GS = ifelse(GS == "No", "NG", GS)) %>% # Normal Glucose
  mutate(GS = ifelse(GS == "Sí", "AG", GS)) %>% # Altered glucose
  # Change Body Mass Index variable
  mutate(BMI = ifelse(BMI == "Sobrepeso 2", "OB", BMI)) %>% # Obesity
  mutate(BMI = ifelse(BMI == "Obesidad", "OB", BMI)) %>% # Obesity
  mutate(BMI = ifelse(BMI == "Normopeso", "NW", BMI)) %>% # Normal Weight
  # Change Insuline Resistanece variable
  mutate(IR = ifelse(IR == "Resistencia significativa", "SR", IR)) %>% # Significant Resistance
  mutate(IR = ifelse(IR == "Resistencia temprana", "ER", IR)) %>% #Early Resistance 
  mutate(IR = ifelse(IR == "Resistencia Temprana", "ER", IR)) %>% #Early Resistance 
  mutate(IR = ifelse(IR == "Rango Saludable", "No_IR", IR)) %>%   # No resistance
  mutate(IR = ifelse(IR == "Rango saludable", "No_IR", IR)) %>%   # No resistance
  column_to_rownames(var = "Sample_ID")

saveRDS(object = metadata, file = "00_preprocess_metadata/data/clinical_data_v2.rds")
