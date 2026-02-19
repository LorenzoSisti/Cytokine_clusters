pacman::p_load(
  bio3d,
  tidyverse,   
  readxl,
  skimr,
  janitor,
  FactoMineR,
  factoextra,
  cluster,
  clValid
)

# Setting the directories

work_dir <- "/Users/lorenzosisti/Cytokine_clusters/"

# Functions (Questa parte la devo modificare bene, vorrei creare una pipeline decente che permetta automaticamente di sostituire i valori problematici.)

check_OOR_columns <- function(df) {
  
  for (j in seq_along(df)) {
    col <- df[[j]]
    
    if (any(grepl("^OOR\\s*[<>]", 
                  toupper(trimws(as.character(col)))))) {
      cat("⚠️  Colonna con OOR:", names(df)[j], "\n")
    }
  }
  
}

check_OOR_columns <- function(df) {
  
  threshold <- (nrow(df) / 2) + 1
  
  for (j in seq_along(df)) {
    col <- df[[j]]
    
    is_oor <- grepl("^OOR\\s*[<>]",
                    toupper(trimws(as.character(col))))
    
    n_oor <- sum(is_oor, na.rm = TRUE)
    
    if (n_oor > threshold) {
      cat("⚠️  Colonna con OOR > (nrow/2)+1:", 
          names(df)[j],
          "| n_OOR =", n_oor, "\n")
    }
  }
  
}


# Importing the citokine files (manually cleaned from Excel)

ck_patients <- read_xlsx(paste(work_dir, "ck_patients.xlsx", sep = ""))
ck_ctrl <- read_xlsx(paste(work_dir, "ck_ctrl.xlsx", sep = ""))

check_OOR_columns(ck_patients)
check_OOR_columns(ck_ctrl)

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


# Calcolo una media per i valori duplicati facendo un group_by (questa cosa si può sostituire con una funzione if ci sono pazienti duplicati then fai una media)
avg_df <- merged_df %>%
  group_by(Patient) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), # Mantiene i metadati non numerici
    across(where(is.character), first)
  ) %>% 
  ungroup()  # chiude il contesto di gruppo

#avg_df<- merged_df
data_log <- avg_df %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))


# Trasformazione log2 
#data_log <- log2(avg_df + 1)

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

# Questa parte qui sotto va rifatta bene, potrebbe essere ottimizzata 

clmethods <- c("hierarchical","kmeans")
clmetrics <- c("euclidean","correlation")
agglos <- c("ward","single","complete","average")

res <- list()

for (m in clmetrics) {
  # kmeans non usa "method", quindi ha senso separare
  res[[paste0("kmeans_", m)]] <- clValid(
    pca_coords, nClust = 2:6,
    clMethods = "kmeans",
    validation = "internal",
    metric = m
  )
  
  for (a in agglos) {
    res[[paste0("hier_", a, "_", m)]] <- clValid(
      pca_coords, nClust = 2:6,
      clMethods = "hierarchical",
      validation = "internal",
      metric = m,
      method = a
    )
  }
}

summary(res)

###########
# 1. Estraiamo le matrici delle metriche da ogni oggetto clValid
# Usiamo slot(x, "measures") perché clValid è un oggetto S4
list_of_df <- lapply(names(res), function(name) {
  # Estraiamo la matrice delle misure (metriche x numero cluster)
  mat <- res[[name]]@measures
  
  # La trasformiamo in dataframe "lungo" per renderla leggibile
  df <- as.data.frame(mat)
  df$Metric_Type <- rownames(mat)
  df$Configuration <- name
  
  return(df)
})

# 2. Uniamo tutto in un unico grande dataframe
final_results <- do.call(rbind, list_of_df)

# 3. Riordiniamo le colonne per leggibilità
final_results <- final_results[, c("Configuration", "Metric_Type", "2", "3", "4", "5", "6")]

# Visualizziamo il risultato
print(final_results)

############




# esempio: guarda un oggetto
res[["hier_ward_euclidean"]]
summary(res[["hier_ward_euclidean"]])


### Clustering based on euclidean distance between observtions

# Partitioning method: k-means with silhouette evaluation
k_means_partitioning <- fviz_nbclust(pca_coords, kmeans, method = "silhouette", k.max = 20) 

# Hierarchical method (AGNES, average linkage)
hierarchical_partitioning <- agnes(pca_coords, metric = "euclidean", stand = FALSE, method = "average")
fviz_dend(hierarchical_partitioning, show_labels = TRUE)

cophenetic_distance <- cophenetic(hierarchical_partitioning) # Evaluate goodness of clustering: cophenetic correlation
cor(dist(pca_coords, method = "euclidean"), cophenetic_distance)

# In our case, average linkage shows highest cophenetic correlation and was therefore chosen as the default linking method (this is true for all the clustering analysis)

### Clustering based on Pearson correlation distance

# Partitioning method: k-means with Pearson-based distance
k_means_partitioning_pearson <- fviz_nbclust(pca_coords, kmeans, method = "silhouette", diss = get_dist(pca_coords, method = "pearson"), k.max = 20) #La matrice delle distanze tra scaled_data è creata usando la metrica euclidea

# Hierarchical method (average linkage)
# We don't use AGNES anymore because it does not support correlation-based distances
dissimilarity_matrix <- get_dist(pca_coords, method = "pearson")
hierarchical_partitioning_pearson <- hclust(dissimilarity_matrix, method = "average")
fviz_dend(hierarchical_partitioning_pearson, show_labels = FALSE)

# Cophenetic correlation for Pearson-based clustering
cor(get_dist(scaled_data, method = "pearson"), cophenetic(hierarchical_partitioning_pearson))


##############





# usa solo le prime 5 PC (come prima)
dist_rows <- dist(pca_coords[, 1:5], method = "euclidean")

hc_rows <- hclust(dist_rows, method = "average")

clusters_rows <- cutree(hc_rows, k = 2)








# 1. Clustering con metodo Ward.D2
dist_rows <- dist(subset(data_log_df, select = -Patient))
hc_rows <- hclust(dist_rows, method = "ward.D2") 

dist_cols <- dist(t(subset(data_log_df, select = -Patient)))
hc_cols <- hclust(dist_cols, method = "ward.D2") 

# 2. Definizione dei tagli (k rimane lo stesso)
clusters_rows <- cutree(hc_rows, k = 3)
clusters_cols <- cutree(hc_cols, k = 11)

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