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

merged_df$`ENA-78/CXCL5` <- NA

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

data_log <- avg_df %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))

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

# varianza spiegata
var_explained <- (pca_res$sdev)^2

# percentuale
percent_var <- var_explained / sum(var_explained) * 100

percent_var

cum_percent_var <- cumsum(percent_var)

cum_percent_var

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


###############

# Cose sul clustering con clValid che aggiungo in data 19/02

clmethods <- c("hierarchical","kmeans")
clmetric <- c("euclidean","correlation")
agglomeration_method <- c("ward","single","complete","average")

intern <- clValid(pca_coords[, 1:5], nClust = 2:6, clMethods = clmethods, validation = "internal", method = "ward")

summary(intern)

prova <- list()

prova[["kmeans"]] <- clValid(
  pca_coords[, 1:5], nClust = 2:6,
  clMethods = "kmeans",
  validation = "internal",
  metric = "euclidean"
)

for (i in agglomeration_method) {
  prova[[paste0("hier_", i)]] <- clValid(
    pca_coords[, 1:5], nClust = 2:6,
    clMethods = "hierarchical",
    validation = "internal",
    metric = "euclidean",
    method = i
  )
}

# stampare tutti
for (nm in names(prova)) {
  cat("\n====================\n", nm, "\n")
  print(summary(prova[[nm]]))
}

###############

# usa solo le prime 5 PC (spiegano 75% della varianza)
dist_rows <- dist(pca_coords[, 1:5], method = "euclidean")

hc_rows <- hclust(dist_rows, method = "average")

clusters_rows <- cutree(hc_rows, k = 2)

### ora senza pca

X_rows <- data_log_df %>%
  select(-Patient) %>%  # togli la colonna testo
  as.matrix() %>%
  scale()

dist_rows <- dist(X_rows, method = "euclidean")

hc_rows <- hclust(dist_rows, method = "average")

clusters_rows <- cutree(hc_rows, k = 2)

### QUI DEVO INSERIRE UNA FUNZIONE CHE MI CALCOLI IL MIGLIOR MODO DI AGGLOMERARE I DATI (PARTENDO DALLA SILHOUETTE)
### CON LA PCA E PRENDENDO LE PRIME 5 PC LA SILHOUETTE AUMENTA DA 0.27 A 0.38, QUINDI BENE



######################### 
### ORA FACCIO UN CLUSTERING SULLE COLONNE 
########################


X_cols <- data_log_df %>%
  select(-Patient) %>%  # togli la colonna testo
  as.matrix() %>%
  t()

X_cols_scaled <- t(scale(t(X_cols)))


prova_cols <- list()
agglomeration_method <- c("ward","single","complete","average")

# kmeans
prova_cols[["kmeans"]] <- clValid(
  X_cols_scaled,
  nClust = 2:12,
  clMethods = "kmeans",
  validation = "internal",
  metric = "euclidean"
)

# hierarchical con metodi diversi
for (i in agglomeration_method) {
  prova_cols[[paste0("hier_", i)]] <- clValid(
    X_cols_scaled,
    nClust = 2:12,
    clMethods = "hierarchical",
    validation = "internal",
    metric = "euclidean",
    method = i
  )
}

for (nm in names(prova_cols)) {
  cat("\n====================\n", nm, "\n")
  print(summary(prova_cols[[nm]]))
}

### proviamo la correlazione per evitare variazioni di scala assoluta

cor_cluster <- clValid(
  X_cols,
  nClust = 2:12,
  clMethods = "hierarchical",
  validation = "internal",
  metric = "correlation",
  method = "average"
)

summary(cor_cluster)

dist_cols <- dist(X_cols_scaled)
hc_cols <- hclust(dist_cols, method = "average")
clusters_cols <- cutree(hc_cols, k = 2)

### C'è davvero pochissima tendenza al clustering sulle colonne.

dist_cols <- dist(t(subset(data_log_df, select = -Patient)))
hc_cols <- hclust(dist_cols, method = "ward.D2") 

# 2. Definizione dei tagli (k rimane lo stesso)
#clusters_rows <- cutree(hc_rows, k = 3)
clusters_cols <- cutree(hc_cols, k = 11) # ANCORA DEVO CAPIRE A COSA SERVA QUESTO, AD OGNI MODO POI RIPETO LA PIPELINE PER CAPIRE IL MIGLIOR CLUSTERING SULLE COLONNE (CITOCHINE)

# 3. Preparazione annotazioni 
annotation_row <- data.frame(
  Patient_Donor = factor(df_merge$Patient_Donor),
  Gender = factor(df_merge$Sex_final),
  AgeGroup = as.numeric(df_merge$Age_final),
  Cluster_row = factor(clusters_rows, levels = 1:3)
)
rownames(annotation_row) <- df_merge$`Patient ID`

annotation_col <- data.frame(
  Cluster_col = factor(clusters_cols, levels = 1:11)
)
rownames(annotation_col) <- colnames(numeric_data)

# 4. FIX Palette Colori (Mapping esplicito per evitare buchi bianchi)
cols_cluster_row <- brewer.pal(3, "Set1")
names(cols_cluster_row) <- as.character(1:3)

# Per 11 cluster palette più estesa
cols_cluster_col <- colorRampPalette(brewer.pal(11, "Set3"))(11)
names(cols_cluster_col) <- as.character(1:11)

ann_colors <- list(
  Patient_Donor = setNames(brewer.pal(length(levels(annotation_row$Patient_Donor)), "Set2"), levels(annotation_row$Patient_Donor)),
  Gender = setNames(brewer.pal(length(levels(annotation_row$Gender)), "Pastel1"), levels(annotation_row$Gender)),
  Cluster_row = cols_cluster_row,
  Cluster_col = cols_cluster_col
)

# 5. Heatmap finale
pheatmap(
  numeric_data,
  color = color_palette,
  breaks = breaks,
  cluster_rows = hc_rows,
  cluster_cols = hc_cols,
  cutree_rows = 3,
  cutree_cols = 11,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  angle_col = 45,
  border_color = NA,
  fontsize_col = 12,
  fontsize_row = 10, # Abbassato leggermente per leggibilità se hai molti pazienti
  cellwidth = 20,
  cellheight = 15,
  main = "Heatmap (Metodo Ward.D2) - Log2 Data",
  filename = "BIORAD_Citokynes/Heatmap_WardD2_ALL.png")