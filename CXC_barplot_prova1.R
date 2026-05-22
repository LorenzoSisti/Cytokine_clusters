pacman::p_load(
  bio3d,
  tidyverse,   
  readxl,
  skimr,
  janitor,
  FactoMineR,
  factoextra,
  cluster,
  clValid,
  pvclust # questo magari lo uso alla fine
)

# Setting the directories

#work_dir <- "/Users/lorenzosisti/Cytokine_clusters/"

# Functions (Questa parte la devo modificare bene, vorrei creare una pipeline decente che permetta automaticamente di sostituire i valori problematici.)


# Importing the citokine files (manually cleaned from Excel)

ck_patients <- read_xlsx("ck_patients.xlsx")
ck_ctrl <- read_xlsx("ck_ctrl.xlsx")

# Tolgo la citochina problematica per valori OOR

ck_patients <- ck_patients %>%
  mutate(`6Ckine/CCL21` = ifelse(`6Ckine/CCL21` == "OOR >", 
                                 48432.02, 
                                 as.numeric(`6Ckine/CCL21`)))


# Unione dei due dataframe (pz + ctrl)
merged_df <- bind_rows(ck_patients, ck_ctrl)
#clean_names(merged_df)

# 1. Creazione variabile binaria (OK)
merged_df <- merged_df %>%
  mutate(Diagnosis_Binary = ifelse(Diagnosis == "LLC", "LLC", "NHL"))


# Genero una tabella solo colonne numeriche
#merged_df <- pca_df_merge %>% select(-c(Sex_final, Age_final, Patient_Donor,-institute, -`Patient ID`))
merged_df <- merged_df %>% select(-c(Gender, Age, Date_birth,Diagnosis, `ENA-78/CXCL5`))


# Calcolo una media per i valori duplicati facendo un group_by 
avg_df <- merged_df %>%
  group_by(Patient) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), # Mantiene i metadati non numerici
    across(where(is.character), first)
  ) %>% 
  ungroup()  # chiude il contesto di gruppo
