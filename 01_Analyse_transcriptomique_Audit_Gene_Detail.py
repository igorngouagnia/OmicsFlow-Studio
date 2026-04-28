import pandas as pd
import numpy as np
import os

# --- CONFIGURATION ---
METADATA_PATH = r"D:\Stage-CRBS\Stage_analyses\Transcriptomique\metadata.txt"
FILE_PATH = r"D:\Stage-CRBS\Stage_analyses\Transcriptomique\GSE160079_Raw_gene_counts_matrix_Cohort_MTM1-a.txt"
OUTPUT_DIR = r"D:\Stage-CRBS\Stage_analyses\Résultats\Transcriptomique\Audits"

# Cibles (Modifie ces deux valeurs selon tes besoins)
GENE_CIBLE = "ENSMUSG00000000001"
COHORTE = "A"

def run_detailed_audit():
    # Création du dossier d'audit s'il n'existe pas
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # 1. Chargement des données
    meta = pd.read_csv(METADATA_PATH, sep="\t")
    meta['Sample_name'] = meta['Sample_name'].str.strip()
    cohort_meta = meta[meta['Cohort'] == COHORTE]

    df = pd.read_csv(FILE_PATH, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    
    # 2. Calcul des Library Sizes
    library_sizes = df.iloc[:, 1:].sum()

    # 3. Extraction de la ligne du gène
    gene_data = df[df.iloc[:, 0] == GENE_CIBLE]
    if gene_data.empty:
        print(f"Erreur : Le gène {GENE_CIBLE} est introuvable.")
        return

    # 4. Construction du rapport texte
    report_lines = []
    report_lines.append(f"=== RAPPORT D'AUDIT DETAILLE : {GENE_CIBLE} ===")
    report_lines.append(f"Cohorte : {COHORTE}")
    report_lines.append(f"Source : {os.path.basename(FILE_PATH)}\n")
    report_lines.append("-" * 60)

    groupes = ['WT', 'WT treated', 'Disease', 'Rescue']
    final_means = {}

    for grp in groupes:
        samples = cohort_meta[cohort_meta['Group'] == grp]['Sample_name'].tolist()
        report_lines.append(f"\nGROUPE : {grp}")
        report_lines.append(f"{'Echantillon':<20} | {'Reads':<10} | {'Lib Size':<15} | {'CPM':<10}")
        report_lines.append("-" * 60)
        
        cpm_list = []
        for s in samples:
            if s in df.columns:
                raw_count = gene_data[s].values[0]
                total_reads = library_sizes[s]
                cpm = (raw_count / total_reads) * 1_000_000
                cpm_list.append(cpm)
                report_lines.append(f"{s:<20} | {raw_count:<10} | {total_reads:<15,} | {cpm:.4f}")
        
        if cpm_list:
            moyenne = np.mean(cpm_list)
            final_means[grp] = moyenne
            report_lines.append(f"-> MOYENNE {grp} : {moyenne:.4f} CPM")
        else:
            report_lines.append(f"-> AUCUNE DONNEE")

    # 5. Calcul des Ratios (Simulation Log2FC)
    report_lines.append("\n" + "=" * 60)
    report_lines.append("CALCULS DE VALIDATION (Log2FC avec pseudocount 0.1)")
    
    if 'WT' in final_means and 'Disease' in final_means:
        ratio = (final_means['Disease'] + 0.1) / (final_means['WT'] + 0.1)
        l2fc = np.log2(ratio)
        report_lines.append(f"Log2FC Patho (Disease vs WT) : {l2fc:.4f}")

    # 6. Sauvegarde du fichier
    file_name = f"Audit_{COHORTE}_{GENE_CIBLE}.txt"
    full_path = os.path.join(OUTPUT_DIR, file_name)
    
    with open(full_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    
    print(f"Audit terminé. Rapport généré : {full_path}")

if __name__ == "__main__":
    run_detailed_audit()