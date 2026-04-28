import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

def create_volcano(df, name_col, lfc_col, pval_col, title, output_path, lfc_th=1.0, pval_th=0.05):
    # Vérification de l'existence des colonnes
    if not all(col in df.columns for col in [lfc_col, pval_col]):
        print(f"Colonnes absentes pour {title}. Ignoré.")
        return
    
    name_c = name_col if name_col in df.columns else df.columns[0]
    
    df_plot = df.copy()
    
    # Remplacement des valeurs problématiques pour P-value
    df_plot[pval_col] = pd.to_numeric(df_plot[pval_col], errors='coerce')
    df_plot[lfc_col] = pd.to_numeric(df_plot[lfc_col], errors='coerce')
    df_plot = df_plot.dropna(subset=[lfc_col, pval_col])
    if df_plot.empty: 
        return
    
    # Calcul de -log10 P-value (remplace les p=0 en valeurs infinitésimales pour éviter -inf)
    df_plot['-log10(p-value)'] = -np.log10(df_plot[pval_col].replace(0, 1e-300).astype(float))
    
    # Identifier la siginificativité pour les couleurs
    conditions = [
        (df_plot[pval_col] < pval_th) & (df_plot[lfc_col] > lfc_th),
        (df_plot[pval_col] < pval_th) & (df_plot[lfc_col] < -lfc_th)
    ]
    choices = ['Augmenté (Up)', 'Diminué (Down)']
    df_plot['Significatif'] = np.select(conditions, choices, default='Non Significatif')
    
    color_map = {
        'Augmenté (Up)': '#e74c3c', # Rouge
        'Diminué (Down)': '#3498db', # Bleu
        'Non Significatif': '#bdc3c7' # Gris clair
    }
    
    # Génération du Volcano Plot via Plotly Express
    fig = px.scatter(
        df_plot, 
        x=lfc_col, 
        y='-log10(p-value)', 
        color='Significatif', 
        color_discrete_map=color_map,
        hover_name=name_c,
        hover_data={lfc_col: ':.3f', pval_col: ':.3e', '-log10(p-value)': False, 'Significatif': False},
        title=title,
        labels={lfc_col: 'Log2 Fold Change', '-log10(p-value)': '-Log10(P-value)'},
        opacity=0.75
    )
                     
    # Lignes des Seuils
    fig.add_hline(y=-np.log10(pval_th), line_dash="dash", line_color="lightgrey", annotation_text=f" p={pval_th}", annotation_position="top left", annotation_font=dict(color="lightgrey"))
    if lfc_th > 0:
        fig.add_vline(x=lfc_th, line_dash="dash", line_color="lightgrey")
        fig.add_vline(x=-lfc_th, line_dash="dash", line_color="lightgrey")
    
    # Cosmétique
    fig.update_layout(
        template="plotly_dark", 
        title_x=0.5, 
        hovermode='closest',
        title_font=dict(size=20, color='white'),
        legend_title_text=''
    )
    
    # Traces plus esthétiques (bordures des points)
    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='darkgrey')))
    
    fig.write_html(output_path)
    print(f"  -> Sauvé : {os.path.basename(output_path)}")

def main():
    # Use current directory as fallback for base_res to ensure portability
    base_res = os.environ.get("OMICS_OUT_DIR", os.path.join(os.getcwd(), "Résultats"))
    out_dir = os.path.join(base_res, "Volcano_plots")
    
    # ==========================
    # 1. TRANSCRIPTOMIQUE
    # ==========================
    t_dir = os.path.join(out_dir, "Transcriptomique")
    os.makedirs(t_dir, exist_ok=True)
    in_t = os.path.join(base_res, "Transcriptomique")
    print("\n[Transcriptomique]")
    if os.path.exists(in_t):
        for f in os.listdir(in_t):
            if f.startswith("Deseq2_Genes_Analysis_Cohort_") and f.endswith(".csv"):
                cohort = f.split("_")[-1].replace(".csv", "")
                df = pd.read_csv(os.path.join(in_t, f), sep=";", index_col=0) # index_col=0 to keep Ensembl ID
                df['Ensembl_Gene_ID'] = df.index
                
                # Patho (Disease vs Ref) | Threshold LFC=1, Padj=0.05
                create_volcano(df, 'Ensembl_Gene_ID', 'Log2FC_Patho', 'Pvalue_Patho', 
                               f"Transcriptomique Patho (Cohorte {cohort})", 
                               os.path.join(t_dir, f"Volcano_Patho_Cohort_{cohort}.html"), 
                               lfc_th=1.0)
                
                # Rescue (Rescue vs Disease) | Threshold LFC=1, Padj=0.05
                create_volcano(df, 'Ensembl_Gene_ID', 'Log2FC_Rescue', 'Pvalue_Rescue', 
                               f"Transcriptomique Rescue (Cohorte {cohort})", 
                               os.path.join(t_dir, f"Volcano_Rescue_Cohort_{cohort}.html"), 
                               lfc_th=1.0)

    # ==========================
    # 2. PROTÉOMIQUE (PYTHON - VALIDÉ SUPRIYA)
    # ==========================
    pyp_dir = os.path.join(out_dir, "Proteomique_Python")
    os.makedirs(pyp_dir, exist_ok=True)
    sup_pyp_dir = os.path.join(base_res, "Protéomique", "Supriya", "Python")
    print("\n[Protéomique Python (Validé Supriya)]")
    
    for age in ['e18.5', '2w', '7w']:
        f_path = os.path.join(sup_pyp_dir, f"Proteomics_Analysis_Full_Python_{age}.csv")
        if os.path.exists(f_path):
            df = pd.read_csv(f_path)
            name_col = 'Unnamed: 2' if 'Unnamed: 2' in df.columns else df.columns[0]
            
            # Patho
            create_volcano(df, name_col, 'logFC_Patho', 'adj.P.Val_Patho', 
                           f"Protéomique Python Patho ({age})", 
                           os.path.join(pyp_dir, f"Volcano_Patho_{age}.html"), lfc_th=1.0)
            
            # Rescue
            if 'logFC_Rescue' in df.columns:
                create_volcano(df, name_col, 'logFC_Rescue', 'adj.P.Val_Rescue', 
                               f"Protéomique Python Rescue ({age})", 
                               os.path.join(pyp_dir, f"Volcano_Rescue_{age}.html"), lfc_th=1.0)
        else:
            print(f"  Fichier absent : {f_path}")
            
    # ==========================
    # 3. PROTÉOMIQUE (R_DEP - VALIDÉ SUPRIYA)
    # ==========================
    pr_dir = os.path.join(out_dir, "Proteomique_R")
    os.makedirs(pr_dir, exist_ok=True)
    sup_pr_dir = os.path.join(base_res, "Protéomique", "Supriya", "R")
    print("\n[Protéomique R (Validé Supriya)]")
    
    for age in ['e18.5', '2w', '7w']:
        f_path = os.path.join(sup_pr_dir, f"Proteomics_Analysis_Full_R_{age}.csv")
        if os.path.exists(f_path):
            df = pd.read_csv(f_path)
            name_col = 'Unnamed: 2' if 'Unnamed: 2' in df.columns else df.columns[0]
            
            # Patho
            create_volcano(df, name_col, 'logFC_Patho', 'adj.P.Val_Patho', 
                           f"Protéomique R Patho ({age})", 
                           os.path.join(pr_dir, f"Volcano_Patho_{age}.html"), 
                           pval_th=0.05, lfc_th=1.0)
            
            # Rescue
            if 'logFC_Rescue' in df.columns:
                create_volcano(df, name_col, 'logFC_Rescue', 'adj.P.Val_Rescue', 
                               f"Protéomique R Rescue ({age})", 
                               os.path.join(pr_dir, f"Volcano_Rescue_{age}.html"), 
                               pval_th=0.05, lfc_th=1.0)
        else:
            print(f"  Fichier absent : {f_path}")

    # ==========================
    # 4. MÉTABOLOMIQUE
    # ==========================
    m_dir = os.path.join(out_dir, "Metabolomique")
    os.makedirs(m_dir, exist_ok=True)
    in_m = os.path.join(base_res, "Métabolomique", "Metabolomics_Full_Analysis.csv")
    print("\n[Métabolomique]")
    if os.path.exists(in_m):
        df = pd.read_csv(in_m, sep=";")
        name_col = 'CHEMICAL_NAME'
        # Pour le métabolome, le script d'origine ne filtrait pas sur le LFC, on utilise LFC=0 (virtuellement désactivé sur le trace)
        create_volcano(df, name_col, 'Log2FC_Patho', 'Padj_Patho', 
                       "Métabolomique Patho (Disease vs WT)", 
                       os.path.join(m_dir, "Volcano_Patho.html"), 
                       lfc_th=0) 
        
        create_volcano(df, name_col, 'Log2FC_Rescue', 'Padj_Rescue', 
                       "Métabolomique Rescue (Rescue vs Disease)", 
                       os.path.join(m_dir, "Volcano_Rescue.html"), 
                       lfc_th=0)

if __name__ == '__main__':
    print("=== DÉMARRAGE DE LA GÉNÉRATION DES VOLCANO PLOTS ===")
    main()
    print("\n=== TERMINÉ ! Vos magnifiques graphiques interactifs vous attendent ! ===")
