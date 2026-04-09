import pandas as pd
import numpy as np
from scipy import stats
import os

def analyse_metabolomique_complete():
    # --- 1. CONFIGURATION DES CHEMINS ---
    root_dir = os.getcwd()
    dossier_source = os.environ.get("OMICS_IN_DIR", os.path.join(root_dir, "Métabolomique"))
    dossier_sortie = os.environ.get("OMICS_OUT_DIR", os.path.join(root_dir, "Résultats", "Métabolomique"))
    fichier_excel = os.environ.get("OMICS_META_PATH", os.path.join(dossier_source, "IGBM-01-21VW MUSCLE DATA TABLES.XLSX"))
    
    if not os.path.exists(dossier_sortie):
        os.makedirs(dossier_sortie)

    try:
        print("Chargement des données...")
        # Chargement des feuilles avec les noms exacts
        df_log = pd.read_excel(fichier_excel, sheet_name="Log Transformed Data")
        df_meta = pd.read_excel(fichier_excel, sheet_name="Sample Meta Data")
        df_annot = pd.read_excel(fichier_excel, sheet_name="Chemical Annotation")

        # --- 2. DÉFINITION DES GROUPES ---
        samples_wt = df_meta[df_meta['TREATMENT'] == 'WT']['PARENT_SAMPLE_NAME'].tolist()
        samples_dis = df_meta[df_meta['TREATMENT'] == 'Disease']['PARENT_SAMPLE_NAME'].tolist()
        samples_res = df_meta[df_meta['TREATMENT'] == 'Rescue']['PARENT_SAMPLE_NAME'].tolist()

        metabolite_ids = [col for col in df_log.columns if col != 'PARENT_SAMPLE_NAME']
        results = []

        print(f"Analyse de {len(metabolite_ids)} métabolites en cours...")

        # --- 3. CALCULS STATISTIQUES ---
        for m_id in metabolite_ids:
            # Extraction des valeurs
            val_wt = df_log[df_log['PARENT_SAMPLE_NAME'].isin(samples_wt)][m_id].dropna()
            val_dis = df_log[df_log['PARENT_SAMPLE_NAME'].isin(samples_dis)][m_id].dropna()
            val_res = df_log[df_log['PARENT_SAMPLE_NAME'].isin(samples_res)][m_id].dropna()

            # Moyennes (en échelle Log)
            m_wt, m_dis, m_res = val_wt.mean(), val_dis.mean(), val_res.mean()

            # Comparaison Patho (Disease vs WT)
            log2fc_patho = m_dis - m_wt
            _, p_patho = stats.ttest_ind(val_dis, val_wt, equal_var=False) if len(val_dis)>1 else (0, 1)

            # Comparaison Rescue (Rescue vs Disease)
            log2fc_rescue = m_res - m_dis
            _, p_rescue = stats.ttest_ind(val_res, val_dis, equal_var=False) if len(val_res)>1 else (0, 1)

            results.append({
                'CHEM_ID': int(m_id),
                'Mean_WT': m_wt,
                'Mean_Disease': m_dis,
                'Mean_Rescue': m_res,
                'Log2FC_Patho': log2fc_patho,
                'p_value_Patho': p_patho,
                'Log2FC_Rescue': log2fc_rescue,
                'p_value_Rescue': p_rescue
            })

        # --- 4. JOINTURE ET FILTRAGE STRICT ---
        df_stats = pd.DataFrame(results)
        df_final = pd.merge(df_annot[['CHEM_ID', 'CHEMICAL_NAME', 'SUPER_PATHWAY', 'SUB_PATHWAY', 'HMDB', 'KEGG']], 
                            df_stats, on='CHEM_ID')

        # === CORRECTION FDR ET CORRECTION DU TAUX DE FAUX POSITIFS ===
        from statsmodels.stats.multitest import multipletests
        
        # FDR Patho
        mask_p = ~df_final['p_value_Patho'].isna()
        df_final['Padj_Patho'] = 1.0
        if mask_p.any():
            _, padj_p, _, _ = multipletests(df_final.loc[mask_p, 'p_value_Patho'], method='fdr_bh')
            df_final.loc[mask_p, 'Padj_Patho'] = padj_p
            
        # FDR Rescue
        mask_r = ~df_final['p_value_Rescue'].isna()
        df_final['Padj_Rescue'] = 1.0
        if mask_r.any():
            _, padj_r, _, _ = multipletests(df_final.loc[mask_r, 'p_value_Rescue'], method='fdr_bh')
            df_final.loc[mask_r, 'Padj_Rescue'] = padj_r

        # Indicateurs de base
        df_final['Significatif_Patho'] = df_final['Padj_Patho'] < 0.05
        df_final['Significatif_Rescue_Brut'] = df_final['Padj_Rescue'] < 0.05

        # Filtre STRICT pour le Rescue : Significatif dans Patho ET significatif dans Rescue
        df_final['Significatif_Rescue'] = (df_final['Significatif_Patho']) & (df_final['Significatif_Rescue_Brut'])

        # --- 5. EXPORT DES 3 FICHIERS CSV ---
        df_final.to_csv(os.path.join(dossier_sortie, "Metabolomics_Full_Analysis.csv"), sep=';', index=False)
        
        df_patho = df_final[df_final['Significatif_Patho']].copy()
        df_patho.to_csv(os.path.join(dossier_sortie, "Metabolomics_Signature_Patho.csv"), sep=';', index=False)
        
        df_rescue = df_final[df_final['Significatif_Rescue']].copy()
        df_rescue.to_csv(os.path.join(dossier_sortie, "Metabolomics_Signature_Rescue.csv"), sep=';', index=False)

        # --- 6. GÉNÉRATION DU RAPPORT TEXTE ---
        chemin_rapport = os.path.join(dossier_sortie, "Bilan_Analyse_Metabolomique.txt")
        with open(chemin_rapport, "w", encoding="utf-8") as f:
            f.write("=== BILAN DE L'ANALYSE MÉTABOLOMIQUE MTM1 ===\n\n")
            f.write(f"Nombre total de métabolites analysés : {len(df_final)}\n")
            f.write(f"Nombre d'échantillons WT      : {len(samples_wt)}\n")
            f.write(f"Nombre d'échantillons Disease : {len(samples_dis)}\n")
            f.write(f"Nombre d'échantillons Rescue  : {len(samples_res)}\n")
            f.write("-" * 40 + "\n")
            f.write(f"1. SIGNATURE PATHO (Disease vs WT):\n")
            f.write(f"   - Métabolites significatifs (p < 0.05) : {len(df_patho)}\n")
            f.write(f"   - Dont augmentés dans Disease : {len(df_patho[df_patho['Log2FC_Patho'] > 0])}\n")
            f.write(f"   - Dont diminués dans Disease : {len(df_patho[df_patho['Log2FC_Patho'] < 0])}\n\n")
            f.write(f"2. SIGNATURE RESCUE (Filtre Strict):\n")
            f.write(f"   - Métabolites restaurés (Signif Patho ET Signif Rescue) : {len(df_rescue)}\n")
            
            # Calcul du taux de restauration
            restauration_ok = df_rescue[np.sign(df_rescue['Log2FC_Patho']) != np.sign(df_rescue['Log2FC_Rescue'])]
            f.write(f"   - Dont inversion de signe (Restauration réelle) : {len(restauration_ok)}\n")
            f.write("\nAnalyse terminée avec succès.")

        print(f"Analyse terminée. Rapport disponible : {chemin_rapport}")

    except Exception as e:
        print(f"Erreur : {e}")

if __name__ == "__main__":
    analyse_metabolomique_complete()