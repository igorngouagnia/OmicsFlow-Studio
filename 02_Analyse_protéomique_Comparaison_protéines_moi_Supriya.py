import pandas as pd
import os

# Utilisation de variables d'environnement pour la portabilité (Streamlit App)
base_res = os.environ.get("OMICS_OUT_DIR", r"D:\Stage-CRBS\Stage_analyses\Résultats\Protéomique\Supriya")
R_DIR = os.path.join(base_res, "R")
PY_DIR = os.path.join(base_res, "Python")
REF_DIR = os.environ.get("OMICS_REF_DIR", r"D:\Stage-CRBS\Stage_analyses\Protéomique\Données-Supriya")

def compare(age):
    ref_file = os.path.join(REF_DIR, f"results_{age}.xlsx")
    ref_df = pd.read_excel(ref_file, sheet_name=0, header=3)
    
    col_fc = f"KO_{age.upper()}_vs_WT_{age.upper()} logFC"
    if col_fc not in ref_df.columns:
        col_fc = [c for c in ref_df.columns if 'logFC' in str(c)][0]
        
    ref_sig = ref_df[(ref_df['adj.P.Val'] < 0.05) & (ref_df[col_fc].abs() > 0.25)]
    ref_genes = set(ref_sig['Genename'].dropna().str.upper().str.strip())

    
    for pipeline, directory in [('R', R_DIR), ('Python', PY_DIR)]:
        loc_file = os.path.join(directory, f"Proteomics_Signature_Patho_{pipeline}_{age}.csv")
        if not os.path.exists(loc_file):
            print(f"File not found: {loc_file}")
            continue
            
        loc_sig = pd.read_csv(loc_file)
        loc_genes = set(loc_sig['Unnamed: 2'].dropna().str.upper().str.strip())
        
        intersection = loc_genes.intersection(ref_genes)
        only_loc = loc_genes - ref_genes
        only_ref = ref_genes - loc_genes
        
        inter_pct_ref = len(intersection) / len(ref_genes) * 100 if ref_genes else 0
        inter_pct_loc = len(intersection) / len(loc_genes) * 100 if loc_genes else 0
        
        report_file = os.path.join(directory, f"Comparaison_Bilan_Patho_{age}.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(f"=== Bilan Comparaison {pipeline} vs Supriya (Reference) pour {age} ===\n")
            f.write(f"- Implémente le filtre 10% valeurs manquantes per condition\n")
            f.write(f"- (R) VSN + Limma  ou  (Python) Log2 + Median Center + t-test/BH\n")
            f.write(f"Gènes significatifs dans la référence (Supriya) : {len(ref_genes)}\n")
            f.write(f"Gènes significatifs dans notre analyse {pipeline} : {len(loc_genes)}\n")
            f.write(f"Intersection (Gènes communs) : {len(intersection)}\n")
            f.write(f"\n---> Pourcentage d'intersection (sur réf.) : {inter_pct_ref:.2f}%\n")
            f.write(f"---> Pourcentage d'intersection (sur loc.) : {inter_pct_loc:.2f}%\n")
            f.write(f"\nGènes uniquement dans notre analyse ({len(only_loc)} gènes)\n")
            f.write(f"Gènes uniquement dans la référence ({len(only_ref)} gènes)\n")
        
        loc_sig_filtered = loc_sig[loc_sig['Unnamed: 2'].str.upper().str.strip().isin(only_loc)]
        loc_sig_filtered.to_csv(os.path.join(directory, f"Proteines_dans_analyse_absentes_dans_Supriya_Patho_{age}.csv"), index=False)
        
        loc_sig_inter = loc_sig[loc_sig['Unnamed: 2'].str.upper().str.strip().isin(intersection)]
        loc_sig_inter.to_csv(os.path.join(directory, f"Proteines_dans_analyse_et_dans_Supriya_Patho_{age}.csv"), index=False)
        
        ref_sig_filtered = ref_sig[ref_sig['Genename'].str.upper().str.strip().isin(only_ref)]
        ref_sig_filtered.to_csv(os.path.join(directory, f"Proteines_dans_Supriya_absentes_dans_analyse_Patho_{age}.csv"), index=False)
        print(f"[{pipeline} - {age}] Intersect vs Ref: {inter_pct_ref:.2f}% | Validated ? {'Yes' if inter_pct_ref > 90 else 'No'}")

compare('2w')
compare('7w')
