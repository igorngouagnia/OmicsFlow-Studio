import pandas as pd
import os

def comparer_donnees_cohorte():
    chemin_excel = r"D:\Stage-CRBS\Stage_analyses\Tamoxifen improves muscle structure and function of Bin1- and Dnm2-related centronuclear myopathies - Supplementary - data\brain-2021-02002-File011.xlsx"
    base_res = r"D:\Stage-CRBS\Stage_analyses\Résultats\Transcriptomique"
    
    # 1. Lecture Référence
    df_ref = pd.read_excel(chemin_excel, sheet_name="Table S2", skiprows=1)
    col_id_ref = "Ensembl Gene ID"
    ids_ref = set(df_ref[col_id_ref].astype(str).str.strip().unique())

    configs = [
        {
            "path": os.path.join(base_res, "Deseq2", "Deseq2_Genes_Signature_Patho_Cohort_G.csv"),
            "label": "deseq2",
            "out_dir": os.path.join(base_res, "Deseq2")
        },
        {
            "path": os.path.join(base_res, "Python", "Genes_Signature_Patho_Cohort_G.csv"),
            "label": "python",
            "out_dir": os.path.join(base_res, "Python")
        }
    ]

    for config in configs:
        if not os.path.exists(config["path"]): continue
        
        out_dir = config["out_dir"]
        label = config["label"]
        report_path = os.path.join(out_dir, f"01_Analyse_transcriptomique_Comparaison_G_{label}.txt")
        
        df_loc = pd.read_csv(config["path"], sep=";")
        col_id = next(c for c in ["Ensembl_Gene_ID", "Ensembl Gene ID"] if c in df_loc.columns)
        ids_loc = set(df_loc[col_id].astype(str).str.strip().unique())
        
        # Calculs
        communs = ids_loc & ids_ref
        extras = ids_loc - ids_ref
        missing = ids_ref - ids_loc
        match_rate = (len(communs) / len(ids_ref)) * 100 if ids_ref else 0

        # Export CSV
        df_loc[df_loc[col_id].astype(str).str.strip().isin(communs)].to_csv(os.path.join(out_dir, f"genes_dans_{label}_et_dans_fichier_excel_Patho_cohorte_G.csv"), index=False, sep=';')
        df_loc[df_loc[col_id].astype(str).str.strip().isin(extras)].to_csv(os.path.join(out_dir, f"genes_dans_{label}_absents_dans_fichier_excel_Patho_cohorte_G.csv"), index=False, sep=';')
        df_ref[df_ref[col_id_ref].astype(str).str.strip().isin(missing)].to_csv(os.path.join(out_dir, f"genes_dans_fichier_excel_absents_dans_{label}_Patho_cohorte_G.csv"), index=False, sep=';')

        # Export Rapport TXT
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(f"=== RAPPORT DE COMPARAISON TRANSCRIPTOMIQUE - MOTEUR {label.upper()} ===\n\n")
            f.write(f"Référence : Papier (Table S2) - {len(ids_ref)} gènes\n")
            f.write(f"Analyse locale : {label} - {len(ids_loc)} gènes détectés\n\n")
            f.write(f"RÉSULTATS :\n")
            f.write(f"- Gènes Communs (Intersection) : {len(communs)}\n")
            f.write(f"- Gènes Extras (Analyse seule) : {len(extras)}\n")
            f.write(f"- Gènes Manquants (Papier seul) : {len(missing)}\n\n")
            f.write(f"SCORE DE REPRODUCTIBILITÉ : {match_rate:.2f}%\n")
            f.write(f"\nFichiers CSV générés dans : {out_dir}\n")

if __name__ == "__main__":
    comparer_donnees_cohorte()