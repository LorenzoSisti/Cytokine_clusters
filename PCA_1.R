pacman::p_load(
  bio3d,
  tidyverse,   
  readxl,
  skimr,
  janitor,
  FactoMineR,
  factoextra
)

# Setting the directories

work_dir <- "/Users/lorenzosisti/Citokine_ISS/"

# Importing the citokine files (cleaned from Excel)

ck_patients <- read_xlsx(paste(work_dir, "ck_patients.xlsx", sep = ""))
ck_ctrl <- read_xlsx(paste(work_dir, "ck_ctrl.xlsx", sep = ""))

# Tolgo la citochina problematica per valori OOR

ck_patients$`ENA-78/CXCL5` <- NA

ck_patients <- ck_patients %>%
  mutate(`6Ckine/CCL21` = ifelse(`6Ckine/CCL21` == "OOR >", 
                                 48432.02, 
                                 as.numeric(`6Ckine/CCL21`)))


# Unione dei due dataframe
df_num <- bind_rows(ck_patients, ck_ctrl)
#clean_names(df_num)

# Genero una tabella solo colonne numeriche
#df_num <- pca_df_merge %>% select(-c(Sex_final, Age_final, Patient_Donor,-institute, -`Patient ID`))
df_num <- df_num %>% select(-c(Gender, Age, Date_birth,Diagnosis, `ENA-78/CXCL5`))


# Calcolo una media per i valori duplicati facendo un group_by (questa cosa si può sostituire con una funzione if ci sono pazienti duplicati then fai una media)
df_avg <- df_num %>%
  group_by(Patient) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), # Mantiene i metadati non numerici
   across(where(is.character), first)
  ) %>% 
  ungroup()  # chiude il contesto di gruppo

#df_avg<- df_num
data_log <- df_avg %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))


# Trasformazione log2 
#data_log <- log2(df_avg + 1)

# converto in data.frame
data_log_df <- as.data.frame(data_log)

# assegno rownames
rownames(data_log_df) <- make.unique(as.character(data_log_df$Patient))

# prcomp
pca_res <- prcomp(subset(data_log_df, select = -Patient), scale. = TRUE)

# Scree plot (genera una percentuale di varianza spiegata)
fviz_eig(pca_res)

# coordinate dei punti per poter assegnare le etichette in modo più mirato
pca_coords <- as.data.frame(pca_res$x)
#data_log$Patient <- rownames(pca_coords)
#pca_coords$Institute <- pca_df_merge$institute

# soglia per outlier (1 deviazioni standard)
sd_mult <- 1
outliers <- pca_coords %>%
  filter(abs(PC1) > sd_mult * sd(PC1) | abs(PC2) > sd_mult * sd(PC2))

outliers_lab <- outliers %>%
  rownames_to_column("Patient")

# Genero il grafico della PCA 
p <- fviz_pca_ind(
  pca_res,
  geom = "point",
  #col.ind = pca_coords$Institute,  # colori
  pointsize = 2,
  repel = TRUE,
  mean.point = FALSE
)

p + 
  ggrepel::geom_text_repel(
    data = outliers_lab,
    aes(x = PC1, y = PC2, label = Patient),#, color = Institute),
    size = 3
  )


ggsave("BIORAD_Citokynes/PCA_before_remove_BatchEffect_log2.png", width = 8, height = 6, dpi = 300)