import pandas as pd
import numpy as np
import os
from statsmodels.stats.multitest import multipletests
import pingouin as pg

INPUT_DIR = os.environ.get("OMICS_IN_DIR", "Protéomique")
base_out = os.environ.get("OMICS_OUT_DIR", "Résultats/Protéomique/Supriya")
OUTPUT_DIR = os.path.join(base_out, "Python")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Loading and cleaning data...")
df_prot = pd.read_csv(os.path.join(INPUT_DIR, "proteinGroups.tsv"), sep='\t', low_memory=False)
df_prot = df_prot[df_prot['Potential contaminant'] != '+']
df_prot = df_prot[df_prot['Reverse'] != '+']
df_prot = df_prot[df_prot['Only identified by site'] != '+']
df_prot = df_prot[df_prot['Unique peptides'] >= 2]

def run_pipeline(samples_wt, samples_ko, samples_rescue, age_label):
    print(f"Running Python pipeline for {age_label} (Patho & Rescue)...")
    
    all_samples = samples_wt + samples_ko + samples_rescue
    cols = ["LFQ intensity " + s for s in all_samples]
    mat = df_prot[cols].replace(0, np.nan)
    
    # Filtrage : au moins 10% de valeurs présentes dans l'un des groupes
    def get_prop(row, smp_list):
        s_cols = ["LFQ intensity " + s for s in smp_list]
        vals = row[s_cols]
        return np.sum(~np.isnan(vals)) / len(vals) if len(vals) > 0 else 0

    keep = []
    for i, row in mat.iterrows():
        if get_prop(row, samples_wt) >= 0.1 or get_prop(row, samples_ko) >= 0.1 or get_prop(row, samples_rescue) >= 0.1:
            keep.append(i)
            
    mat = mat.loc[keep]
    idx = mat.index
    
    # Log2 and Median Normalization
    mat_log2 = np.log2(mat.astype(float))
    mat_norm = mat_log2 - mat_log2.median() + mat_log2.median().mean()
    
    def get_stats(df_n, s1, s2):
        ps, fcs = [], []
        c1 = ["LFQ intensity " + s for s in s1]
        c2 = ["LFQ intensity " + s for s in s2]
        for _, row in df_n.iterrows():
            v1 = row[c1].dropna()
            v2 = row[c2].dropna()
            if len(v1) > 1 and len(v2) > 1:
                try:
                    res = pg.ttest(v2, v1)
                    ps.append(res['p-val'].values[0])
                    fcs.append(np.mean(v2) - np.mean(v1))
                except:
                    ps.append(1.0); fcs.append(0.0)
            else:
                ps.append(1.0); fcs.append(0.0)
        _, adj_p, _, _ = multipletests(ps, method='fdr_bh')
        return fcs, ps, adj_p

    # Comparaison Patho: KO vs WT
    lfc_p, p_p, padj_p = get_stats(mat_norm, samples_wt, samples_ko)
    # Comparaison Rescue: Rescue vs KO
    lfc_r, p_r, padj_r = get_stats(mat_norm, samples_ko, samples_rescue)
    
    out = pd.DataFrame()
    out['Protein Annotation'] = df_prot.loc[idx, 'Majority protein IDs'].str.split(';').str[0]
    out['Unnamed: 1'] = df_prot.loc[idx, 'Fasta headers'].str.split(';').str[0]
    out['Unnamed: 2'] = df_prot.loc[idx, 'Gene names'].str.split(';').str[0]
    out['Unnamed: 3'] = df_prot.loc[idx, 'Protein names']
    
    for c in cols:
        out[c.replace("LFQ intensity ", "")] = mat_norm[c]
        
    out['logFC'] = lfc_p # Pour compatibilité ascendante (Patho par défaut)
    out['adj.P.Val'] = padj_p
    
    out['logFC_Patho'] = lfc_p
    out['P.Value_Patho'] = p_p
    out['adj.P.Val_Patho'] = padj_p
    
    out['logFC_Rescue'] = lfc_r
    out['P.Value_Rescue'] = p_r
    out['adj.P.Val_Rescue'] = padj_r
    
    out.to_csv(os.path.join(OUTPUT_DIR, f"Proteomics_Analysis_Full_Python_{age_label.lower()}.csv"), index=False)
    
    # Sauvegarde des signatures
    sig_p = out[(out['adj.P.Val_Patho'] < 0.05) & (out['logFC_Patho'].abs() > 0.25)]
    sig_p.to_csv(os.path.join(OUTPUT_DIR, f"Proteomics_Signature_Patho_Python_{age_label.lower()}.csv"), index=False)
    
    if len(samples_rescue) > 0:
        sig_r = out[(out['adj.P.Val_Rescue'] < 0.05) & (out['logFC_Rescue'].abs() > 0.25)]
        sig_r.to_csv(os.path.join(OUTPUT_DIR, f"Proteomics_Signature_Rescue_Python_{age_label.lower()}.csv"), index=False)

# E18.5 Samples
s_wt_e = [f"180808_AM_Ech{i:02d}_rep{j}" for i in range(1, 5) for j in range(1, 4)]
s_ko_e = [f"180808_AM_Ech{i:02d}_rep{j}" for i in range(5, 9) for j in range(1, 4)]
s_re_e = [f"180808_AM_Ech{i:02d}_rep{j}" for i in range(9, 13) for j in range(1, 4)]
run_pipeline(s_wt_e, s_ko_e, s_re_e, "E18.5")

# 2W Samples
s_wt_2 = ["180808_AM_Ech13_rep1", "180808_AM_Ech13_rep2", "180808_AM_Ech13_rep3", 
          "180808_AM_Ech14_rep1", "180808_AM_Ech14_rep2", "180808_AM_Ech14_rep3"]
s_ko_2 = [f"180808_AM_Ech{i}_rep{j}" for i in range(17, 21) for j in range(1, 4)]
s_re_2 = [f"180808_AM_Ech{i}_rep{j}" for i in range(21, 25) for j in range(1, 4)]
run_pipeline(s_wt_2, s_ko_2, s_re_2, "2W")

# 7W Samples
s_wt_7 = [f"180808_AM_Ech{i}_rep{j}" for i in range(26, 30) for j in range(1, 4)]
s_ko_7 = [f"180808_AM_Ech{i}_rep{j}" for i in range(30, 34) for j in range(1, 4)]
s_re_7 = [f"180808_AM_Ech{i}_rep{j}" for i in range(34, 37) for j in range(1, 4)] + ["180808_AM_Ech25_rep1", "180808_AM_Ech25_rep2", "180808_AM_Ech25_rep3"]
run_pipeline(s_wt_7, s_ko_7, s_re_7, "7W")

print("Python Script Executed Successfully")
