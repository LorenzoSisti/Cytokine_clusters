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
  pvclust,
  ggpubr
)

# 1. IMPORT UNA TIE-IN DEI FILE
ck_patients <- read_xlsx("ck_patients.xlsx")
ck_ctrl <- read_xlsx("ck_ctrl.xlsx")

# FIX: Assegniamo esplicitamente la diagnosi ai controlli prima dell'unione
ck_ctrl <- ck_ctrl %>% mutate(Diagnosis = "Donatore Sano")

# Tolgo la citochina problematica per valori OOR
ck_patients <- ck_patients %>%
  mutate(`6Ckine/CCL21` = ifelse(`6Ckine/CCL21` == "OOR >", 
                                 48432.02, 
                                 as.numeric(`6Ckine/CCL21`)))

# Unione dei due dataframe (pz + ctrl)
merged_df <- bind_rows(ck_patients, ck_ctrl)

# Gestione colonna vuota
merged_df$`ENA-78/CXCL5` <- NA

# FIX GRUPPI: Creiamo i 3 macro-gruppi puliti lasciando intatta la diagnosi originale per i colori dei sottotipi
merged_df <- merged_df %>%
  mutate(
    Diagnosis_Group = case_when(
      Diagnosis == "LLC" ~ "LLC",
      Diagnosis == "Donatore Sano" ~ "Donatore Sano",
      TRUE ~ "Linfoma" # Tutti gli altri saranno i vari sottotipi di linfoma (NHL, ecc.)
    )
  )

# 2. SELEZIONE COLONNE (Rimuovi i dati sensibili)
merged_df <- merged_df %>% 
  select(-c(Gender, Age, Date_birth, `ENA-78/CXCL5`))

# 3. AGGREGAZIONE
avg_df <- merged_df %>%
  group_by(Patient) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), 
    Diagnosis_Group = first(Diagnosis_Group),
    Diagnosis_Original = first(Diagnosis) 
  ) %>% 
  ungroup()

# 4. TRASFORMAZIONE LOG2 (Escludendo le colonne di testo)
data_log <- avg_df %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))

# 5. CONVERSIONE IN DATAFRAME E ROWNAMES
data_log_df <- as.data.frame(data_log)
rownames(data_log_df) <- make.unique(as.character(data_log_df$Patient))


# ==========================================================================
# PARTE 1: ANALISI PCA (Aggiornata con i 3 gruppi reali)
# ==========================================================================

# Esegui la PCA (escludendo le colonne di testo memorizzate nel DF)
pca_res <- prcomp(subset(data_log_df, select = -c(Patient, Diagnosis_Group, Diagnosis_Original)), scale. = TRUE)

# Scree plot
fviz_eig(pca_res)

pca_coords <- as.data.frame(pca_res$x)
var_explained <- (pca_res$sdev)^2
percent_var <- var_explained / sum(var_explained) * 100

# Individuazione Outliers
sd_mult <- 1
outliers <- pca_coords %>%
  filter(abs(PC1) > sd_mult * sd(PC1) | abs(PC2) > sd_mult * sd(PC2))

outliers_lab <- outliers %>%
  rownames_to_column("Patient")

# Grafico della PCA aggiornato a 3 colori (LLC, Linfomi, Sani)
p <- fviz_pca_ind(
  pca_res,
  geom = "point",
  col.ind = data_log_df$Diagnosis_Group, 
  palette = c("#E41A1C", "#377EB8", "#4DAF4A"), # Rosso=LLC, Blu=Linfoma, Verde=Sani
  addEllipses = TRUE,                                
  ellipse.type = "convex",
  legend.title = "Gruppo Diagnostico",
  pointsize = 2.5,
  repel = TRUE
)

p + ggrepel::geom_text_repel(
  data = outliers_lab,
  aes(x = PC1, y = PC2, label = Patient),
  size = 3
)

ggsave("PCA_before_remove_BatchEffect_log2.png", width = 8, height = 6, dpi = 300)


# ==========================================================================
# PARTE 2: NUOVO MODULO - BOXPLOT AUTOMATICI PER OGNI CITOCHINA
# ==========================================================================

# 1. Trasformiamo i dati in formato "Long" (ottimale per ggplot2)
data_long <- data_log %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Cytokine",
    values_to = "Expression"
  )

# 2. Creiamo la cartella di destinazione per i Boxplots se non esiste
#dir.create("BIORAD_Citokynes/Boxplots", showWarnings = FALSE, recursive = TRUE)

# 3. Estraiamo la lista unica delle citochine presenti
lista_citochine <- unique(data_long$Cytokine)

# 4. Loop per generare e salvare i grafici singolarmente
for (citochina in lista_citochine) {
  
  # Filtriamo il dataset per la singola citochina del ciclo corrente
  df_sub <- data_long %>% filter(Cytokine == citochina)
  
  # Generazione del Boxplot + Jitter delle singole diagnosi
  p_box <- ggplot(df_sub, aes(x = Diagnosis_Group, y = Expression)) +
    # Boxplot macro-gruppo (trasparente per far vedere i punti sotto)
    geom_boxplot(aes(fill = Diagnosis_Group), outlier.shape = NA, alpha = 0.3, show.legend = FALSE) +
    # Punti dei singoli pazienti colorati per diagnosi specifica (Linfomi diversi = Colori diversi)
    geom_jitter(aes(color = Diagnosis_Original), width = 0.2, size = 2.5, alpha = 0.8) +
    # Palette di colori automatica e accattivante per i punti
    scale_color_brewer(palette = "Set1") + 
    scale_fill_manual(values = c("LLC" = "#E41A1C", "Linfoma" = "#377EB8", "Donatore Sano" = "#4DAF4A")) +
    labs(
      title = paste("Espressione di", citochina),
      x = "Gruppo di Studio",
      y = "Log2 (Espressione + 1)",
      color = "Diagnosi Specifica"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # Rimuoviamo caratteri speciali/slash dal nome della citochina per evitare errori nel file system
  nome_file_pulito <- gsub("[[:punct:][:space:]]", "_", citochina)
  
  # Salvataggio del grafico in formato PNG ad alta risoluzione
  ggsave(
    filename = paste0("Boxplots/", nome_file_pulito, "_boxplot.png"),
    plot = p_box,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ==========================================================================
# NUOVO: GRAFICO CUMULATIVO (TUTTI I BOXPLOT IN UN'UNICA IMMAGINE)
# ==========================================================================

# Calcoliamo quante citochine abbiamo per calcolare proporzioni sensate
num_citochine <- length(lista_citochine)

p_all <- ggplot(data_long, aes(x = Diagnosis_Group, y = Expression)) +
  # Boxplot di sfondo per il macro-gruppo
  geom_boxplot(aes(fill = Diagnosis_Group), outlier.shape = NA, alpha = 0.3, show.legend = FALSE) +
  # Punti colorati per sottotipo di linfoma
  geom_jitter(aes(color = Diagnosis_Original), width = 0.2, size = 1, alpha = 0.7) +
  # Dividiamo il grafico in sotto-grafici per ogni Citochina
  # scales = "free_y" permette ad ogni citochina di avere la sua scala di espressione corretta
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 5) + 
  scale_color_brewer(palette = "Set1") + 
  scale_fill_manual(values = c("LLC" = "#E41A1C", "Linfoma" = "#377EB8", "Donatore Sano" = "#4DAF4A")) +
  labs(
    title = "Profilo Espressivo Globale delle Citochine",
    x = "Gruppo di Studio",
    y = "Log2 (Espressione + 1)",
    color = "Diagnosi Specifica"
  ) +
  theme_bw() + # Sfondo bianco, ottimo per griglie dense
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 8), # Testo dei titoli delle citochine
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7), # Ruotiamo i testi sull'asse X per non sovrapporli
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Salvataggio dell'immagine panoramica
# Nota: Dinamicamente impostiamo l'altezza in base a quante citochine hai nel file
righe_necessarie <- ceiling(num_citochine / 5)

ggsave(
  filename = "Boxplots/00_TUTTE_le_citochine_panel.png",
  plot = p_all,
  width = 15,                                # Molto largo per ospitare 5 colonne
  height = max(4, righe_necessarie * 2.5),   # Si allunga automaticamente in base al numero di righe
  dpi = 300
)

# ==========================================================================
# NUOVO: GRAFICO CUMULATIVO CON MACRO-CATEGORIA NHL (SENZA SOTTOTIPI)
# ==========================================================================

# 1. Modifichiamo leggermente i livelli nel dataset per avere "NHL" invece di "Linfoma"
data_long <- data_long %>%
  mutate(
    Diagnosis_Group = ifelse(Diagnosis_Group == "Linfoma", "NHL", Diagnosis_Group)
  )

# Calcoliamo quante citochine abbiamo per calcolare proporzioni sensate
lista_citochine <- unique(data_long$Cytokine)
num_citochine <- length(lista_citochine)

p_all <- ggplot(data_long, aes(x = Diagnosis_Group, y = Expression)) +
  # Boxplot di sfondo per il macro-gruppo
  geom_boxplot(aes(fill = Diagnosis_Group), outlier.shape = NA, alpha = 0.3, show.legend = FALSE) +
  # Punti colorati ESCLUSIVAMENTE in base alla macro-categoria (LLC, NHL, Donatore Sano)
  geom_jitter(aes(color = Diagnosis_Group), width = 0.2, size = 1, alpha = 0.7) +
  # Dividiamo il grafico in sotto-grafici per ogni Citochina con scala libera sulla Y
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 5) + 
  # Palette colori coerente per Boxplot e Punti
  scale_color_manual(values = c("LLC" = "#4DAF4A", "NHL" = "#377EB8", "Donatore Sano" = "#E41A1C")) +
  scale_fill_manual(values = c("LLC" = "#4DAF4A", "NHL" = "#377EB8", "Donatore Sano" = "#E41A1C")) +
  labs(
    title = "Profilo Espressivo Globale delle Citochine",
    x = "Gruppo Diagnostico",
    y = "Log2 (Espressione + 1)",
    color = "Gruppo Diagnostico"
  ) +
  theme_bw() + # Sfondo bianco, pulito per la griglia cumulativa
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 8), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Salvataggio dell'immagine panoramica
righe_necessarie <- ceiling(num_citochine / 5)

ggsave(
  filename = "Boxplots/00_TUTTE_le_citochine_panel_nuovo.png",
  plot = p_all,
  width = 15,                                
  height = max(4, righe_necessarie * 2.5),   
  dpi = 300
)

# ==========================================================================
# GRAFICO CUMULATIVO CON ASTERISCHI DI SIGNIFICATIVITÀ
# ==========================================================================

# [Assicurati che data_long abbia i 3 gruppi come nel passaggio precedente]
data_long <- data_long %>%
  mutate(Diagnosis_Group = ifelse(Diagnosis_Group == "Linfoma", "NHL", Diagnosis_Group))

# Definiamo esplicitamente SOLO i confronti richiesti
# Questo dice a ggplot di tracciare le barre di significatività solo tra queste coppie
miei_confronti <- list(
  c("LLC", "Donatore Sano"),
  c("NHL", "Donatore Sano")
)

# Calcoliamo quante citochine abbiamo
lista_citochine <- unique(data_long$Cytokine)
num_citochine <- length(lista_citochine)

p_all_stat <- ggplot(data_long, aes(x = Diagnosis_Group, y = Expression)) +
  # Boxplot di sfondo
  geom_boxplot(aes(fill = Diagnosis_Group), outlier.shape = NA, alpha = 0.3, show.legend = FALSE) +
  # Punti (pazienti)
  geom_jitter(aes(color = Diagnosis_Group), width = 0.2, size = 0.8, alpha = 0.6, show.legend = FALSE) +
  
  # ------------------------------------------------------------------------
# CORE STATISTICO: Aggiunge barre e asterischi per i confronti a coppie
# metodo = "wilcox.test" (Non-parametrico, consigliato per citochine)
# Se preferisci il t-test classico parametrico, cambia in: method = "t.test"
# ------------------------------------------------------------------------
stat_compare_means(
  comparisons = miei_confronti, 
  method = "wilcox.test", 
  label = "p.signif",          # Mostra gli asterischi (*, **, ***, ns) invece del valore numerico
  hide.ns = FALSE,             # Mostra "ns" se la differenza non è significativa
  size = 3                     # Dimensione del testo degli asterischi
) +
  
  # Aggiunge opzionalmente il p-value globale dell'ANOVA/Kruskal-Wallis in alto a ogni quadratino
  stat_compare_means(method = "kruskal.test", label.y = max(data_long$Expression) * 1.1, size = 2.5) +
  
  # Dividiamo in pannelli per ogni citochina
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 5) + 
  
  # Palette colori coerente
  scale_color_manual(values = c("LLC" = "#4DAF4A", "NHL" = "#377EB8", "Donatore Sano" = "#E41A1C")) +
  scale_fill_manual(values = c("LLC" = "#4DAF4A", "NHL" = "#377EB8", "Donatore Sano" = "#E41A1C")) +
  
  labs(
    title = "Analisi Comparativa dell'Espressione delle Citochine",
    subtitle = "Confronti statistici eseguiti tramite Wilcoxon Test (Mann-Whitney)",
    x = "Gruppo Diagnostico",
    y = "Log2 (Espressione + 1)"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic"),
    strip.text = element_text(face = "bold", size = 8), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

# Salvataggio dell'immagine panoramica statistica
righe_necessarie <- ceiling(num_citochine / 5)

ggsave(
  filename = "Boxplots/00_TUTTE_le_citochine_STATISTICHE.png",
  plot = p_all_stat,
  width = 16,                                
  height = max(5, righe_necessarie * 3.2),   # Leggermente più alto per fare spazio alle barre statistiche
  dpi = 300
)
