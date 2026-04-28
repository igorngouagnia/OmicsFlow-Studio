# Setup Local Library Path for permissions
local_lib <- file.path(getwd(), "R_libs")
if(!dir.exists(local_lib)) dir.create(local_lib)
.libPaths(c(local_lib, .libPaths()))

# Setup Bioconductor repositories manually
bioc_repos <- c(
  CRAN = "https://cloud.r-project.org",
  BioCsoft = "https://bioconductor.org/packages/3.22/bioc",
  BioCann = "https://bioconductor.org/packages/3.22/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.22/data/experiment"
)
options(repos = bioc_repos)

if (!require("limma", quietly = TRUE, lib.loc = local_lib)) {
  install.packages("limma", lib = local_lib, quiet=TRUE)
  library(limma, lib.loc = local_lib)
}

if (!require("vsn", quietly = TRUE, lib.loc = local_lib)) {
  install.packages("vsn", lib = local_lib, quiet=TRUE)
  library(vsn, lib.loc = local_lib)
}

input_dir  <- Sys.getenv("OMICS_IN_DIR", "D:/Stage-CRBS/Stage_analyses/Protéomique")
base_out   <- Sys.getenv("OMICS_OUT_DIR", "D:/Stage-CRBS/Stage_analyses/Résultats/Protéomique/Supriya")
output_dir <- file.path(base_out, "R")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Data loading with base R for stability
df_prot <- read.table(file.path(input_dir, "proteinGroups.tsv"), sep="\t", header=TRUE, check.names = FALSE, quote="", fill=TRUE)

# Cleaning
df_prot_clean <- df_prot[df_prot$`Potential contaminant` != "+" & 
                           df_prot$Reverse != "+" & 
                           df_prot$`Only identified by site` != "+" & 
                           df_prot$`Unique peptides` >= 2, ]

data_unique <- df_prot_clean
uniprot_id <- sapply(strsplit(as.character(data_unique$`Fasta headers`), ";"), `[`, 1)
uniprot_id[is.na(uniprot_id)] <- paste0("unknown_", seq_len(sum(is.na(uniprot_id))))
data_unique$uniprot_id <- make.unique(uniprot_id)
rownames(data_unique) <- data_unique$uniprot_id

run_pipeline <- function(samples_wt, samples_ko, samples_rescue, age_label) {
  cat(paste0("🚀 Running R pipeline for ", age_label, "...\n"))
  
  all_samples <- c(samples_wt, samples_ko, samples_rescue)
  cols <- paste0("LFQ intensity ", all_samples)
  
  # Check if columns exist
  missing_cols <- cols[!(cols %in% colnames(data_unique))]
  if(length(missing_cols) > 0) {
    cat(paste0("⚠️ Warning: Missing columns for ", age_label, ": ", paste(missing_cols, collapse=", "), "\n"))
    return(NULL)
  }
  
  mat <- as.matrix(data_unique[, cols])
  mat[mat == 0] <- NA
  
  # Grouping
  conds <- c(rep("WT", length(samples_wt)), 
             rep("KO", length(samples_ko)), 
             rep("Rescue", length(samples_rescue)))
  
  # Filtering
  keep <- apply(mat, 1, function(x) {
    p_wt <- sum(!is.na(x[conds == "WT"])) / length(samples_wt)
    p_ko <- sum(!is.na(x[conds == "KO"])) / length(samples_ko)
    p_re <- if(length(samples_rescue) > 0) sum(!is.na(x[conds == "Rescue"])) / length(samples_rescue) else 0
    return(p_wt >= 0.1 | p_ko >= 0.1 | p_re >= 0.1)
  })
  
  mat_filtered <- mat[keep, ]
  if(nrow(mat_filtered) == 0) return(NULL)
  
  # Normalization
  mat_norm <- normalizeVSN(mat_filtered)
  
  # Limma
  condition <- factor(conds, levels=c("WT", "KO", "Rescue"))
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("WT", "KO", "Rescue")
  fit_lm <- lmFit(mat_norm, design)
  
  # Contrasts
  if(length(samples_rescue) > 0) {
    cont_matrix <- makeContrasts(Patho = KO - WT, Rescue = Rescue - KO, levels=design)
  } else {
    cont_matrix <- makeContrasts(Patho = KO - WT, levels=design)
  }
  
  fit2 <- contrasts.fit(fit_lm, cont_matrix)
  fit2 <- eBayes(fit2)
  
  # Results extraction
  res_p <- topTable(fit2, coef="Patho", sort.by="none", n=Inf)
  
  out <- data.frame(
    `Protein Annotation` = sapply(strsplit(as.character(data_unique[rownames(res_p), "Majority protein IDs"]), ";"), `[`, 1),
    `Unnamed: 1` = rownames(res_p),
    `Unnamed: 2` = sapply(strsplit(as.character(data_unique[rownames(res_p), "Gene names"]), ";"), `[`, 1),
    `Unnamed: 3` = data_unique[rownames(res_p), "Protein names"],
    check.names = FALSE
  )
  out <- cbind(out, mat_norm)
  
  out$logFC_Patho <- res_p$logFC
  out$adj.P.Val_Patho <- res_p$adj.P.Val
  
  if(length(samples_rescue) > 0) {
    res_r <- topTable(fit2, coef="Rescue", sort.by="none", n=Inf)
    out$logFC_Rescue <- res_r$logFC
    out$adj.P.Val_Rescue <- res_r$adj.P.Val
  }
  
  # Export
  write.csv(out, file.path(output_dir, paste0("Proteomics_Analysis_Full_R_", tolower(age_label), ".csv")), row.names=FALSE)
  cat(paste0("✅ Done for ", age_label, "\n"))
}

# --- SAMPLE DEFINITIONS ---
# E18.5
s_wt_e <- paste0("180808_AM_Ech", sprintf("%02d", 1:4), "_rep", rep(1:3, each=4))
s_ko_e <- paste0("180808_AM_Ech", sprintf("%02d", 5:8), "_rep", rep(1:3, each=4))
s_re_e <- paste0("180808_AM_Ech", sprintf("%02d", 9:12), "_rep", rep(1:3, each=4))
run_pipeline(s_wt_e, s_ko_e, s_re_e, "E18.5")

# 2W
s_wt_2 <- c("180808_AM_Ech13_rep1", "180808_AM_Ech13_rep2", "180808_AM_Ech13_rep3", 
            "180808_AM_Ech14_rep1", "180808_AM_Ech14_rep2", "180808_AM_Ech14_rep3")
s_ko_2 <- paste0("180808_AM_Ech", 17:20, "_rep", rep(1:3, each=4))
s_re_2 <- paste0("180808_AM_Ech", 21:24, "_rep", rep(1:3, each=4))
run_pipeline(s_wt_2, s_ko_2, s_re_2, "2W")

# 7W
s_wt_7 <- paste0("180808_AM_Ech", 26:29, "_rep", rep(1:3, each=4))
s_ko_7 <- paste0("180808_AM_Ech", 30:33, "_rep", rep(1:3, each=4))
s_re_7 <- c(paste0("180808_AM_Ech", 34:36, "_rep", rep(1:3, each=3)), 
            "180808_AM_Ech25_rep1", "180808_AM_Ech25_rep2", "180808_AM_Ech25_rep3")
run_pipeline(s_wt_7, s_ko_7, s_re_7, "7W")

cat("\nPipeline R terminé avec succès.\n")
