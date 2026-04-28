# --- 1. CONFIGURATION ---
library(tidyverse)
library(DEP)
library(SummarizedExperiment)

input_dir  <- Sys.getenv("OMICS_IN_DIR", unset=paste0(getwd(), "/Protéomique"))
output_dir <- Sys.getenv("OMICS_OUT_DIR", unset=paste0(getwd(), "/Résultats/Protéomique/R_DEP"))

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Lecture des fichiers
df_meta <- read.table(file.path(input_dir, "metadata.tsv"), sep="\t", header=TRUE, check.names = FALSE)
colnames(df_meta) <- c("Raw_ID", "Sample_name", "Short_ID", "Age", "Genotype")

df_prot <- read.table(file.path(input_dir, "proteinGroups.tsv"), sep="\t", header=TRUE, check.names = FALSE)

# Nettoyage classique
df_prot_clean <- df_prot[df_prot$`Potential contaminant` != "+" & 
                           df_prot$Reverse != "+" & 
                           df_prot$`Only identified by site` != "+", ]

# Création de noms uniques
data_unique <- make_unique(df_prot_clean, "Gene names", "Protein IDs", delim = ";")

# Synchronisation stricte des identifiants
colnames(data_unique) <- make.names(colnames(data_unique))
df_meta$Raw_ID        <- make.names(df_meta$Raw_ID)

ages <- unique(df_meta$Age)
all_results_list <- list()

# --- 2. BOUCLE D'ANALYSE ---
for (age in ages) {
  cat("\n--- Analyse Age :", age, " ---")
  
  sub_meta <- df_meta[df_meta$Age == age, ]
  
  exp_design <- sub_meta %>%
    group_by(Genotype) %>%
    mutate(replicate = row_number()) %>%
    ungroup() %>%
    mutate(label = paste0(Genotype, "_", replicate)) %>%
    select(label, condition = Genotype, replicate, Raw_ID)
  
  columns_idx <- match(exp_design$Raw_ID, colnames(data_unique))
  
  if(any(is.na(columns_idx))) {
    cat("\nERREUR : Noms manquants.\n")
    next
  }
  
  # Correction DEP
  colnames(data_unique)[columns_idx] <- exp_design$label
  
  # Pipeline DEP
  data_se <- make_se(data_unique, columns_idx, exp_design)
  
  # On passe le filtre drastique, puis OMISSION de normalize_vsn 
  # car l'algorithme "vsn" compresse les écarts (Log2FC), ce qui contredit
  # l'analyse de LFQ bruts/déjà normalisés par MaxQuant requise par le papier.
  data_filt <- data_se
  data_norm <- data_filt

  data_imp  <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
  data_diff <- test_diff(data_imp, type = "all")
  dep_results <- add_rejections(data_diff, alpha = 0.05, lfc = 1)
  
  # Extraction
  res <- get_results(dep_results)
  final_age <- as.data.frame(res)
  rownames(final_age) <- NULL
  
  # DEP supprime "Protein names", on le ré-attache depuis data_unique
  if("Protein names" %in% colnames(data_unique)) {
    final_age <- final_age %>% left_join(data_unique %>% select(name, `Protein names`), by="name")
  } else {
    final_age$`Protein names` <- NA
  }
  
  cols <- colnames(final_age)
  
  # Renommage Patho/Rescue
  col_patho_lfc <- grep("KO_vs_WT_log2FoldChange|WT_vs_KO_log2FoldChange|KO_vs_WT_ratio|WT_vs_KO_ratio", cols, value=TRUE)
  col_patho_p   <- grep("KO_vs_WT_p.val|WT_vs_KO_p.val", cols, value=TRUE)
  col_patho_adj <- grep("KO_vs_WT_p.adj|WT_vs_KO_p.adj", cols, value=TRUE)
  
  col_rescue_lfc <- grep("KO_Dnm2_vs_KO_log2FoldChange|KO_vs_KO_Dnm2_log2FoldChange|KO_Dnm2_vs_KO_ratio|KO_vs_KO_Dnm2_ratio", cols, value=TRUE)
  col_rescue_p   <- grep("KO_Dnm2_vs_KO_p.val|KO_vs_KO_Dnm2_p.val", cols, value=TRUE)
  col_rescue_adj <- grep("KO_Dnm2_vs_KO_p.adj|KO_vs_KO_Dnm2_p.adj", cols, value=TRUE)
  
  if(length(col_patho_lfc) > 0) names(final_age)[names(final_age) == col_patho_lfc[1]] <- "LFC_Patho" else final_age$LFC_Patho <- NA
  if(length(col_patho_p) > 0)   names(final_age)[names(final_age) == col_patho_p[1]]   <- "P_Patho"  else final_age$P_Patho <- NA
  if(length(col_patho_adj) > 0) names(final_age)[names(final_age) == col_patho_adj[1]] <- "Padj_Patho"else final_age$Padj_Patho <- NA
  
  if(length(col_rescue_lfc) > 0) names(final_age)[names(final_age) == col_rescue_lfc[1]] <- "LFC_Rescue"else final_age$LFC_Rescue <- NA
  if(length(col_rescue_p) > 0)   names(final_age)[names(final_age) == col_rescue_p[1]]   <- "P_Rescue"  else final_age$P_Rescue <- NA
  if(length(col_rescue_adj) > 0) names(final_age)[names(final_age) == col_rescue_adj[1]] <- "Padj_Rescue"else final_age$Padj_Rescue <- NA
  
  # Correction signes
  if(length(col_patho_lfc) > 0 && grepl("^WT_vs_KO", col_patho_lfc[1])) final_age$LFC_Patho <- -final_age$LFC_Patho
  if(length(col_rescue_lfc) > 0 && grepl("^KO_vs_KO_Dnm2", col_rescue_lfc[1])) final_age$LFC_Rescue <- -final_age$LFC_Rescue
  
  # Calculs signification : RESTAURATION STRICTE DES CRITÈRES DU PAPIER
  is_patho  <- !is.na(final_age$LFC_Patho) & !is.na(final_age$Padj_Patho) & abs(final_age$LFC_Patho) > 1 & final_age$Padj_Patho < 0.05
  opposed   <- sign(final_age$LFC_Patho) != sign(final_age$LFC_Rescue)
  is_rescue <- is_patho & !is.na(final_age$Padj_Rescue) & !is.na(final_age$LFC_Rescue) & final_age$Padj_Rescue < 0.05 & opposed
  
  cat(sprintf("\n  -> Patho : %d | Rescue : %d", sum(is_patho, na.rm=T), sum(is_rescue, na.rm=T)))
  
  # --- CRÉATION DU FORMAT "WIDE" EXACT REQUIS POUR CET ÂGE ---
  age_df <- tibble(
    `Protein IDs`   = final_age$ID,
    `Gene names`    = final_age$name,
    `Protein names` = final_age$`Protein names`
  )
  
  # Moyennes LFQ (récupérées des colonnes WT, KO, KO_Dnm2 générées par DEP)
  age_df[[paste0("Moyenne_LFQ_", age, "_WT")]] <- if("WT" %in% cols) final_age$WT else NA
  age_df[[paste0("Moyenne_LFQ_", age, "_KO")]] <- if("KO" %in% cols) final_age$KO else NA
  age_df[[paste0("Moyenne_LFQ_", age, "_KO_Dnm2")]] <- if("KO_Dnm2" %in% cols) final_age$KO_Dnm2 else NA
  
  # Métriques Patho
  age_df[[paste0("Log2FC_Patho_", age)]] <- final_age$LFC_Patho
  age_df[[paste0("Pvalue_Patho_", age)]] <- final_age$P_Patho # Utilisé la Pvalue standard ou Padj selon besoin. Ici Pvalue.
  
  # Métriques Rescue
  age_df[[paste0("Log2FC_Rescue_", age)]] <- final_age$LFC_Rescue
  age_df[[paste0("Pvalue_Rescue_", age)]] <- final_age$P_Rescue 
  
  # Significations (OUI/NON)
  age_df[[paste0("Significatif_Patho_", age)]]  <- ifelse(is_patho, "OUI", "NON")
  age_df[[paste0("Significatif_Rescue_", age)]] <- ifelse(is_rescue, "OUI", "NON")
  
  # Exporter les fichiers par âge individuels
  write.table(age_df %>% filter(get(paste0("Significatif_Patho_", age)) == "OUI"), 
              file.path(output_dir, paste0("Proteomics_Signature_Patho_", age, ".csv")), sep=";", row.names=F, na="")
  
  write.table(age_df %>% filter(get(paste0("Significatif_Rescue_", age)) == "OUI"), 
              file.path(output_dir, paste0("Proteomics_Efficacite_Rescue_", age, ".csv")), sep=";", row.names=F, na="")
  
  all_results_list[[age]] <- age_df
}

# --- 3. FUSION HORIZONTALE DU GLOBAL ---
# Utilisation explicite de purrr::reduce pour contourner le conflit avec Bioconductor
full_export <- purrr::reduce(all_results_list, function(df1, df2) {
  full_join(df1, df2, by = c("Protein IDs", "Gene names", "Protein names"))
})

write.table(full_export, file.path(output_dir, "Proteomics_Analysis_Full_MTM1.csv"), 
            sep=";", row.names=F, na="")

cat("\n\nTerminé ! Fichier Global fusionné avec succès.\n")

# ==============================================================================
# --- 4. COMPARAISON AVEC LES RÉSULTATS DU PAPIER (COHORTE A / TABLE S8) ---
# ==============================================================================
cat("\n=== DÉBUT DE LA COMPARAISON AVEC L'ARTICLE ===\n")

if (!require("readxl", quietly = TRUE)) {
  cat("\nInstallation/Chargement de readxl requis pour la comparaison...\n")
  install.packages("readxl", repos = "http://cran.us.r-project.org")
  library(readxl)
}

chemin_excel <- "D:/Stage-CRBS/Stage_analyses/Multi-omics comparisons of different forms of centronuclear myopathies and the effects of several therapeutic strategies - supplementary - data/mmc2.xlsx"
df_s8 <- read_excel(chemin_excel, sheet = "Table S8", skip = 2)

report_file <- file.path(output_dir, "02_Analyse_protéomique_Comparaison_protéines_moi_papier_Cohorte_A_DEP.txt")
fileConn <- file(report_file, open="wt")
writeLines("=== RAPPORT DE COMPARAISON PROTÉOMIQUE R (MTM1-a) ===", fileConn)
writeLines("Paramètres stricts reproduits depuis le papier original :", fileConn)
writeLines("- Imputation : de type wrProteo/Perseus (remplacement des valeurs manquantes par décalage à l'extrême gauche de la distribution)", fileConn)
writeLines("- Normalisation : Brut (LFQ conservés, non écrasés par VSN)", fileConn)
writeLines("- Modèle statistique : Bayes Empirique (via DEP/limma remplaçant wrProteo)", fileConn)
writeLines("- Cutoffs Différentiels : P-value ajustée (FDR) < 0.05 ET |Log2FC| > 1", fileConn)
writeLines(rep("=", 60), fileConn)

config_age <- list(
  "2w" = list(col_idx = 2),
  "7w" = list(col_idx = 6)
)

for(age_comp in names(config_age)) {
  writeLines(sprintf("\n--- ANALYSE AGE : %s ---", age_comp), fileConn)
  cat(sprintf("\nComparaison pour %s...\n", age_comp))
  
  col_idx <- config_age[[age_comp]]$col_idx
  nom_col_ref <- colnames(df_s8)[col_idx]
  
  # 1. Gènes du papier
  raw_ref <- as.character(df_s8[[col_idx]])
  raw_ref <- raw_ref[!is.na(raw_ref) & raw_ref != "NA" & trimws(raw_ref) != ""]
  genes_ref_upper <- unique(toupper(trimws(raw_ref)))
  
  writeLines(sprintf("Référence Papier (%s) : %d gènes uniques trouvés.", age_comp, length(genes_ref_upper)), fileConn)
  
  # 2. Gènes de l'analyse R (Patho)
  sig_patho_file <- file.path(output_dir, paste0("Proteomics_Signature_Patho_", age_comp, ".csv"))
  if(!file.exists(sig_patho_file)) {
      writeLines(sprintf("Fichier local introuvable : %s", basename(sig_patho_file)), fileConn)
      next
  }
  
  df_local <- read.delim(sig_patho_file, sep=";", check.names=FALSE, stringsAsFactors=FALSE)
  
  col_gene_local <- if("Gene names" %in% colnames(df_local)) "Gene names" else "Gene.names"
  
  genes_locaux_list <- strsplit(as.character(df_local[[col_gene_local]]), ";")
  genes_locaux_upper <- unique(toupper(trimws(unlist(genes_locaux_list))))
  genes_locaux_upper <- genes_locaux_upper[!is.na(genes_locaux_upper) & genes_locaux_upper != ""]
  
  # 3. Métriques
  communs <- intersect(genes_locaux_upper, genes_ref_upper)
  en_trop <- setdiff(genes_locaux_upper, genes_ref_upper)
  manquants <- setdiff(genes_ref_upper, genes_locaux_upper)
  
  writeLines(sprintf("Vos résultats DEP (%s) : %d gènes identifiés.", age_comp, length(genes_locaux_upper)), fileConn)
  writeLines(sprintf("Intersection : %d gènes en commun.", length(communs)), fileConn)
  writeLines(sprintf("Extras (Chez vous uniquement) : %d gènes.", length(en_trop)), fileConn)
  writeLines(sprintf("Manquants (Papier uniquement) : %d gènes.", length(manquants)), fileConn)
  
  cat(sprintf("  -> Intersection: %d | Extras: %d | Manquants: %d\n", length(communs), length(en_trop), length(manquants)))
  
  # 4. Exports des csv 
  is_commun <- sapply(strsplit(as.character(df_local[[col_gene_local]]), ";"), function(x) any(toupper(trimws(x)) %in% communs))
  write.table(df_local[is_commun, ], file.path(output_dir, paste0("Proteines_dans_analyse_et_dans_fichier_excel_Patho_", age_comp, ".csv")), sep=";", row.names=FALSE, quote=FALSE, na="")
  
  is_extra <- sapply(strsplit(as.character(df_local[[col_gene_local]]), ";"), function(x) any(toupper(trimws(x)) %in% en_trop))
  write.table(df_local[is_extra, ], file.path(output_dir, paste0("Proteines_dans_analyse_absentes_dans_fichier_excel_Patho_", age_comp, ".csv")), sep=";", row.names=FALSE, quote=FALSE, na="")
  
  df_s8_miss <- df_s8[toupper(trimws(as.character(df_s8[[col_idx]]))) %in% manquants, ]
  write.table(df_s8_miss, file.path(output_dir, paste0("Proteines_dans_fichier_excel_absentes_dans_analyse_Patho_", age_comp, ".csv")), sep=";", row.names=FALSE, quote=FALSE, na="")
}

writeLines("\nTraitement terminé. Les 6 fichiers de comparaison ont été générés avec succès.", fileConn)
close(fileConn)
cat("=== COMPARAISON TERMINÉE ===\n")
