import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re
from matplotlib.lines import Line2D

# --- CONFIGURATION ---
METADATA_PATH = r"D:\Stage-CRBS\Stage_analyses\Transcriptomique\metadata.txt"
INPUT_DIR = r"D:\Stage-CRBS\Stage_analyses\Transcriptomique"
OUTPUT_DIR = r"D:\Stage-CRBS\Stage_analyses\Résultats\Transcriptomique"

COHORT_FILES = {
    'A': "GSE160079_Raw_gene_counts_matrix_Cohort_MTM1-a.txt",
    'B': "GSE160081_Raw_gene_counts_matrix_Cohort_MTM1-b.txt",
    'C': "GSE160083_Raw_gene_counts_matrix_Cohort_MTM1-c.txt",
    'D': "GSE160077_Raw_gene_counts_matrix_Cohort_BIN1.txt",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt",
    'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

def normalize_name(name):
    if pd.isna(name): return ""
    return re.sub(r'[^A-Z0-9]', '', str(name).upper())

def load_data_and_pca():
    print("Chargement des données...")
    meta = pd.read_csv(METADATA_PATH, sep='\t')
    if 'Outlier' in meta.columns:
        meta = meta[meta['Outlier'] != 'Yes']
    
    meta_map = {normalize_name(n): n for n in meta['Sample_name']}
    all_counts = {}

    for cid, fname in COHORT_FILES.items():
        path = os.path.join(INPUT_DIR, fname)
        if not os.path.exists(path): continue
        
        if fname.endswith('.xlsx'):
            df = pd.read_excel(path, skiprows=1, index_col=0)
        else:
            df = pd.read_csv(path, sep=None, engine='python', index_col=0)
        
        df = df.select_dtypes(include=[np.number])
        for col in df.columns:
            norm = normalize_name(col)
            if norm in meta_map:
                all_counts[meta_map[norm]] = df[col]

    counts = pd.DataFrame(all_counts).fillna(0)
    final_meta = meta[meta['Sample_name'].isin(counts.columns)].set_index('Sample_name')
    counts = counts[final_meta.index]
    
    cpm = (counts / counts.sum()) * 1e6
    log_counts = np.log2(cpm + 1)
    
    top_genes = log_counts.var(axis=1).sort_values(ascending=False).head(5000).index
    scaled = StandardScaler().fit_transform(log_counts.loc[top_genes].T)
    
    pca = PCA(n_components=2)
    coords = pca.fit_transform(scaled)
    
    df_pca = pd.DataFrame(coords, columns=['PC1', 'PC2'], index=final_meta.index)
    return pd.concat([df_pca, final_meta], axis=1), pca.explained_variance_ratio_

def setup_plot_style():
    sns.set_style("ticks")
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42

def plot_figure_1a(df, var_exp):
    """Figure 1A: Cohortes (Couleurs) et Groupes (Formes)"""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(11, 7))
    
    cohort_colors = {
        'A': '#000000', 'B': '#E69F00', 'C': '#56B4E9', 
        'D': '#009E73', 'E': '#F0E442', 'F': '#D55E00', 'G': '#0072B2'
    }
    group_markers = {
        'WT': 'o', 'WT treated': '*', 'Disease': 's', 'Disease treated': 'X'
    }

    for grp, marker in group_markers.items():
        subset = df[df['Group'] == grp]
        for coh, color in cohort_colors.items():
            sub_coh = subset[subset['Cohort'] == coh]
            if not sub_coh.empty:
                ax.scatter(sub_coh['PC1'], sub_coh['PC2'],
                           c=color, marker=marker, s=100, 
                           edgecolor='black', linewidth=0.5, alpha=0.8)

    ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)", fontweight='bold')
    ax.set_title("Figure 1A", loc='left', fontweight='bold', fontsize=14)

    # --- LÉGENDES VERTICALES ---
    cohort_handles = [Line2D([0], [0], marker='s', color='w', label=f'Cohort {k}', 
                             markerfacecolor=v, markersize=10) for k, v in cohort_colors.items()]
    group_handles = [Line2D([0], [0], marker=v, color='w', label=k, 
                            markerfacecolor='grey', markeredgecolor='black', markersize=10) for k, v in group_markers.items()]

    # Première légende (Cohortes)
    leg1 = ax.legend(handles=cohort_handles, title="Cohorts", loc='upper left', 
                      bbox_to_anchor=(1.02, 1), frameon=False)
    ax.add_artist(leg1)
    
    # Deuxième légende (Groupes) juste en dessous
    ax.legend(handles=group_handles, title="Groups", loc='upper left', 
              bbox_to_anchor=(1.02, 0.5), frameon=False)

    sns.despine()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1A_PCA.png"), dpi=300, bbox_inches='tight')
    plt.close()

def plot_figure_1b(df, var_exp):
    """Figure 1B: Groupes (Couleurs) et Modèles (Formes)"""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(11, 7))
    
    group_colors = {
        'WT': '#228B22', 'WT treated': '#90EE90', 'Disease': '#B22222', 'Disease treated': '#FF8C00'
    }
    model_markers = {
        'MTM1': 'o', 'BIN1': 's', 'DNM2': 'X'
    }

    for mdl, marker in model_markers.items():
        subset = df[df['Model'] == mdl]
        for grp, color in group_colors.items():
            sub_grp = subset[subset['Group'] == grp]
            if not sub_grp.empty:
                ax.scatter(sub_grp['PC1'], sub_grp['PC2'],
                           c=color, marker=marker, s=100,
                           edgecolor='black', linewidth=0.5, alpha=0.8)

    ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)", fontweight='bold')
    ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)", fontweight='bold')
    ax.set_title("Figure 1B", loc='left', fontweight='bold', fontsize=14)

    # --- LÉGENDES VERTICALES ---
    group_handles = [Line2D([0], [0], marker='s', color='w', label=k, 
                            markerfacecolor=v, markersize=10) for k, v in group_colors.items()]
    model_handles = [Line2D([0], [0], marker=v, color='w', label=k, 
                             markerfacecolor='grey', markeredgecolor='black', markersize=10) for k, v in model_markers.items()]

    leg1 = ax.legend(handles=group_handles, title="Groups", loc='upper left', 
                      bbox_to_anchor=(1.02, 1), frameon=False)
    ax.add_artist(leg1)
    
    ax.legend(handles=model_handles, title="Models", loc='upper left', 
              bbox_to_anchor=(1.02, 0.5), frameon=False)

    sns.despine()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1B_PCA.png"), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    try:
        data, var_ratio = load_data_and_pca()
        plot_figure_1a(data, var_ratio)
        plot_figure_1b(data, var_ratio)
        print(f"Terminé. Fichiers sauvegardés dans : {OUTPUT_DIR}")
    except Exception as e:
        print(f"Erreur : {e}")