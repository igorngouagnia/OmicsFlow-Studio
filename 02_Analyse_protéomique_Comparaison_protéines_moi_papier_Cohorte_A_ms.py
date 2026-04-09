import pandas as pd
import os

def comparer_proteomique_mtm1a_ms():
    # --- CONFIGURATION DES CHEMINS ---
    chemin_excel = r"D:\Stage-CRBS\Stage_analyses\Multi-omics comparisons of different forms of centronuclear myopathies and the effects of several therapeutic strategies - supplementary - data\mmc2.xlsx"
    # Dossier où se trouvent les fichiers _ms.csv générés précédemment
    dossier_entree = r"D:\Stage-CRBS\Stage_analyses\Résultats\Protéomique"
    dossier_sortie = r"D:\Stage-CRBS\Stage_analyses\Résultats\Protéomique"
    
    # Nom du rapport mis à jour pour MS/MS
    nom_rapport_txt = "02_Analyse_protéomique_Comparaison_MSMS_moi_papier_Cohorte_A.txt"
    chemin_rapport = os.path.join(dossier_sortie, nom_rapport_txt)

    # Configuration des comparaisons utilisant les fichiers avec le suffixe _ms
    config_age = {
        "2w": {"fichier": "Proteomics_Signature_Patho_2w_ms.csv", "col_excel_idx": 1},
        "7w": {"fichier": "Proteomics_Signature_Patho_7w_ms.csv", "col_excel_idx": 5}
    }

    def logger(message, file_ptr):
        print(message)
        file_ptr.write(message + "\n")

    try:
        with open(chemin_rapport, "w", encoding="utf-8") as f_txt:
            logger("=== RAPPORT DE COMPARAISON PROTÉOMIQUE (MODE MS/MS COUNT) ===", f_txt)
            logger(f"Date de l'analyse : {pd.Timestamp.now()}\n", f_txt)
            
            # 1. LECTURE DE L'EXCEL (TABLE S8)
            df_s8 = pd.read_excel(chemin_excel, sheet_name="Table S8", skiprows=2)

            for age, info in config_age.items():
                logger(f"\n--- ANALYSE AGE : {age} (Fichier source : {info['fichier']}) ---", f_txt)
                
                # Extraction et normalisation des gènes de référence (Papier)
                nom_col_ref = df_s8.columns[info['col_excel_idx']]
                raw_ref = df_s8[nom_col_ref].dropna().astype(str).str.strip()
                genes_ref_upper = set(raw_ref.str.upper().unique())
                
                logger(f"Référence Papier ({age}) : {len(genes_ref_upper)} gènes uniques trouvés.", f_txt)

                # 2. LECTURE DE VOTRE FICHIER MS/MS
                path_csv = os.path.join(dossier_entree, info['fichier'])
                if not os.path.isfile(path_csv):
                    logger(f"ERREUR : Fichier local introuvable : {info['fichier']}", f_txt)
                    continue

                df_local = pd.read_csv(path_csv, sep=';')
                col_gene_local = next((c for c in ["Gene names", "Gene Identifier"] if c in df_local.columns), None)
                
                if not col_gene_local:
                    logger(f"Erreur : Pas de colonne de gènes dans {info['fichier']}", f_txt)
                    continue

                # Extraction de vos gènes (en gérant le séparateur ';' habituel de MaxQuant)
                genes_locaux_upper = set()
                for g in df_local[col_gene_local].dropna().astype(str):
                    parts = [p.strip().upper() for p in g.split(';')]
                    genes_locaux_upper.update(parts)

                # 3. COMPARAISON (Intersection, Extras, Manquants)
                communs = genes_locaux_upper.intersection(genes_ref_upper)
                en_trop = genes_locaux_upper - genes_ref_upper
                manquants = genes_ref_upper - genes_locaux_upper

                logger(f"Vos résultats MS/MS ({age}) : {len(genes_locaux_upper)} gènes identifiés.", f_txt)
                logger(f"Intersection : {len(communs)} gènes en commun.", f_txt)
                logger(f"Extras (Chez vous uniquement) : {len(en_trop)} gènes.", f_txt)
                logger(f"Manquants (Papier uniquement) : {len(manquants)} gènes.", f_txt)

                # 4. EXPORT DES 3 FICHIERS PAR ÂGE (Suffixe _ms conservé)
                
                # A. Intersection
                mask_communs = df_local[col_gene_local].apply(
                    lambda x: any(n.strip().upper() in communs for n in str(x).split(';'))
                )
                df_communs = df_local[mask_communs]
                df_communs.to_csv(os.path.join(dossier_sortie, f"Proteines_analyse_MS_et_excel_Patho_{age}_ms.csv"), sep=';', index=False)
                
                # B. Extras
                mask_extra = df_local[col_gene_local].apply(
                    lambda x: any(n.strip().upper() in en_trop for n in str(x).split(';'))
                )
                df_extra = df_local[mask_extra]
                df_extra.to_csv(os.path.join(dossier_sortie, f"Proteines_analyse_MS_absentes_excel_Patho_{age}_ms.csv"), sep=';', index=False)
                
                # C. Manquants (ceux du papier absents dans votre analyse MS)
                df_miss = df_s8[df_s8[nom_col_ref].astype(str).str.upper().isin(manquants)]
                df_miss.to_csv(os.path.join(dossier_sortie, f"Proteines_excel_absentes_analyse_MS_Patho_{age}_ms.csv"), sep=';', index=False)

            logger("\nComparaison MS/MS terminée. Les fichiers de résultats ont été générés.", f_txt)

    except Exception as e:
        print(f"Une erreur est survenue lors de l'analyse : {e}")

if __name__ == "__main__":
    comparer_proteomique_mtm1a_ms()