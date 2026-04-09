import pandas as pd
import os
import sys

def comparer_donnees_cohorte():
    # --- CONFIGURATION DES CHEMINS PORTABLES ---
    root_dir = os.getcwd()
    chemin_excel = os.environ.get("OMICS_REF_FILE", os.path.join(root_dir, "Multi-omics comparisons of different forms of centronuclear myopathies and the effects of several therapeutic strategies - supplementary - data", "mmc2.xlsx"))
    dossier_entree = os.environ.get("OMICS_OUT_DIR", os.path.join(root_dir, "Résultats", "Transcriptomique"))
    dossier_sortie = os.environ.get("OMICS_OUT_DIR", os.path.join(root_dir, "Résultats", "Transcriptomique"))
    
    # Nom du fichier de rapport à générer
    nom_rapport_txt = "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.txt"
    chemin_rapport = os.path.join(dossier_sortie, nom_rapport_txt)

    # Mapping des fichiers d'entrée vers leurs identifiants pour le nommage de sortie
    fichiers_config = {
        "Deseq2_Genes_Signature_Patho_Cohort_G.csv": "deseq2_standard",
        "Deseq2_Genes_Signature_Patho_Cohort_G_ref_WT_treated.csv": "deseq2_ref_WT_treated",
        "Genes_Signature_Patho_Cohort_G.csv": "python_manual"
    }

    # Fonction interne pour imprimer ET écrire dans le fichier
    def logger(message, file_ptr):
        print(message)
        file_ptr.write(message + "\n")

    try:
        with open(chemin_rapport, "w", encoding="utf-8") as f_txt:
            logger("=== RAPPORT DE COMPARAISON TRANSCRIPTOMIQUE ===", f_txt)
            
            # 1. LECTURE DE LA RÉFÉRENCE (TABLE S2)
            logger("Lecture de la Table S2 dans l'Excel...", f_txt)
            df_ref = pd.read_excel(chemin_excel, sheet_name="Table S2", skiprows=1)
            
            col_id_ref = "Ensembl Gene ID"
            ids_ref = set(df_ref[col_id_ref].astype(str).str.strip().unique())
            logger(f"-> {len(ids_ref)} gènes uniques trouvés dans le fichier Excel (S2).", f_txt)

            # 2. TRAITEMENT DE CHAQUE FICHIER
            for nom_fichier, label in fichiers_config.items():
                path_csv = os.path.join(dossier_entree, nom_fichier)
                
                if not os.path.isfile(path_csv):
                    logger(f"Fichier introuvable : {nom_fichier}", f_txt)
                    continue

                logger(f"\nAnalyse de : {nom_fichier} ({label})", f_txt)
                
                # Lecture du CSV (détection auto du séparateur)
                df_local = pd.read_csv(path_csv, sep=None, engine='python')
                
                # Identification de la colonne ID (cherche les variantes communes)
                col_id_local = next((c for c in ["Ensembl_Gene_ID", "Ensembl Gene ID", "gene_id"] if c in df_local.columns), None)
                
                if not col_id_local:
                    logger(f"Erreur : Colonne ID non trouvée dans {nom_fichier}", f_txt)
                    continue

                ids_locaux = set(df_local[col_id_local].astype(str).str.strip().unique())
                
                # Calcul des différences
                en_trop = ids_locaux - ids_ref
                manquants = ids_ref - ids_locaux

                # 3. EXPORT DES 3 FICHIERS PAR CONFIGURATION
                nom_sortie_communs = f"Genes_dans_analyse_et_dans_fichier_excel_Patho_cohorte_G_{label}.csv"
                df_communs = df_local[df_local[col_id_local].astype(str).str.strip().isin(ids_locaux & ids_ref)]
                df_communs.to_csv(os.path.join(dossier_sortie, nom_sortie_communs), index=False, sep=';')

                nom_sortie_extras = f"genes_dans_{label}_absents_dans_fichier_excel_Patho_cohorte_G.csv"
                df_en_trop = df_local[df_local[col_id_local].astype(str).str.strip().isin(en_trop)]
                df_en_trop.to_csv(os.path.join(dossier_sortie, nom_sortie_extras), index=False, sep=';')
                
                nom_sortie_manquants = f"genes_dans_fichier_excel_absents_dans_{label}_Patho_cohorte_G.csv"
                df_manquants = df_ref[df_ref[col_id_ref].astype(str).str.strip().isin(manquants)]
                df_manquants.to_csv(os.path.join(dossier_sortie, nom_sortie_manquants), index=False, sep=';')

                logger(f"   -> Créé : {nom_sortie_communs} ({len(communs)} gènes)", f_txt)
                logger(f"   -> Créé : {nom_sortie_extras} ({len(en_trop)} gènes)", f_txt)
                logger(f"   -> Créé : {nom_sortie_manquants} ({len(manquants)} gènes)", f_txt)

            logger("\nTraitement terminé. Les 6 fichiers ont été générés avec succès.", f_txt)
            logger(f"Le rapport a été enregistré dans : {nom_rapport_txt}", f_txt)

    except Exception as e:
        print(f"Une erreur est survenue : {e}")

if __name__ == "__main__":
    comparer_donnees_cohorte()