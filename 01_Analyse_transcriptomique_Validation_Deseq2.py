import os
import pandas as pd
import numpy as np
import warnings
import sys
import re
from datetime import datetime

# --- VÉRIFICATION DES BIBLIOTHÈQUES ---
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    print("\nERREUR : PyDESeq2 n'est pas détecté. Lancez : pip install pydeseq2")
    sys.exit(1)

try:
    from pptx import Presentation
    from pptx.util import Inches, Pt
    from pptx.enum.text import PP_ALIGN
    from pptx.dml.color import RGBColor
except ImportError:
    print("\nERREUR : python-pptx n'est pas détecté. Lancez : pip install python-pptx")
    sys.exit(1)

# --- CONFIGURATION DES CHEMINS PORTABLES ---
ROOT_DIR = os.getcwd()
INPUT_DIR = os.environ.get("OMICS_IN_DIR", os.path.join(ROOT_DIR, "Transcriptomique"))
METADATA_PATH = os.environ.get("OMICS_META_PATH", os.path.join(INPUT_DIR, "metadata.txt"))
OUTPUT_DIR = os.environ.get("OMICS_OUT_DIR", os.path.join(ROOT_DIR, "Résultats", "Transcriptomique"))

warnings.filterwarnings("ignore")

COHORT_MAPPING = {
    'A': "GSE160079_Raw_gene_counts_matrix_Cohort_MTM1-a.txt",
    'B': "GSE160081_Raw_gene_counts_matrix_Cohort_MTM1-b.txt",
    'C': "GSE160083_Raw_gene_counts_matrix_Cohort_MTM1-c.txt",
    'D': "GSE160077_Raw_gene_counts_matrix_Cohort_BIN1.txt",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt",
    'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

def normalize_name(name):
    return re.sub(r'[^a-zA-Z0-9]', '', str(name)).lower()

def load_metadata():
    if not os.path.exists(METADATA_PATH):
        raise FileNotFoundError(f"Metadata non trouvé: {METADATA_PATH}")
    df = pd.read_csv(METADATA_PATH, sep="\t")
    df['Sample_name'] = df['Sample_name'].str.strip()
    return df

def calculate_cpm_means(counts_df, metadata_df):
    cpms = counts_df.div(counts_df.sum(axis=0), axis=1) * 1e6
    means = {}
    for group in ['WT', 'WT treated', 'Disease', 'Rescue']:
        samples = metadata_df[metadata_df['Group'] == group]['Sample_name'].tolist()
        valid_samples = [s for s in samples if s in cpms.columns]
        col_name = f'Moyenne_CPM_{group.replace(" ", "_")}'
        if valid_samples:
            means[col_name] = cpms[valid_samples].mean(axis=1)
        else:
            means[col_name] = 0.0
    return pd.DataFrame(means, index=counts_df.index)

def run_deseq2_logic(counts_df, metadata_df, contrast_num, contrast_den):
    sub_meta = metadata_df[metadata_df['Group'].isin([contrast_num, contrast_den])].copy()
    sub_meta.set_index('Sample_name', inplace=True)
    available = [c for c in sub_meta.index if c in counts_df.columns]
    
    if len(available) < 2: return None
    
    sub_meta = sub_meta.loc[available]
    sub_counts = counts_df[available].T
    sub_counts = sub_counts.loc[:, (sub_counts.sum(axis=0) > 0)]
    
    if sub_counts.empty or len(sub_meta['Group'].unique()) < 2:
        return None

    try:
        dds = DeseqDataSet(counts=sub_counts, metadata=sub_meta, design_factors="Group", quiet=True)
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["Group", contrast_num, contrast_den], quiet=True)
        stat_res.summary() 
        if hasattr(stat_res, 'results_df'):
            return stat_res.results_df
        return stat_res.results 
    except Exception as e:
        print(f"   [Info] Échec DESeq2 ({contrast_num} vs {contrast_den}): {e}")
        return None

def create_global_pptx(all_stats):
    prs = Presentation()
    title_slide_layout = prs.slide_layouts[0]
    slide = prs.slides.add_slide(title_slide_layout)
    slide.shapes.title.text = "Rapport d'Analyse Transcriptomique"
    slide.placeholders[1].text = f"Bilan Multi-Cohortes : DESeq2 & Efficacité Rescue\nDate : {datetime.now().strftime('%d/%m/%Y')}"

    blank_slide_layout = prs.slide_layouts[5]
    for stat in all_stats:
        slide = prs.slides.add_slide(blank_slide_layout)
        slide.shapes.title.text = f"Synthèse Cohorte {stat['id']}"

        data = [
            ["Indicateur", "Valeur"],
            ["Échantillons WT / WT Treated", f"{stat['n_wt']} / {stat['n_wtt']}"],
            ["Échantillons Disease / Rescue", f"{stat['n_dis']} / {stat['n_res']}"],
            ["Référence Choisie", stat['ref_grp']],
            ["Signature Patho (|L2FC|>1 & P<0.05)", str(stat['n_patho_signif'])],
            ["Ratio Patho / Total Gènes", f"{stat['ratio_patho']:.2%}"],
            ["Efficacité Rescue (Patho + Rescue Signif.)", str(stat['n_rescue_signif'])],
            ["Ratio Rescue / Patho", f"{stat['ratio_rescue']:.1%}"]
        ]

        rows, cols = 8, 2
        left, top, width, height = Inches(1), Inches(1.5), Inches(8), Inches(4.5)
        table = slide.shapes.add_table(rows, cols, left, top, width, height).table

        for r in range(rows):
            for c in range(cols):
                cell = table.cell(r, c)
                cell.text = data[r][c]
                paragraph = cell.text_frame.paragraphs[0]
                paragraph.font.size = Pt(18)
                paragraph.alignment = PP_ALIGN.CENTER
                if r == 0:
                    cell.fill.solid()
                    # Correction ici : fore_color est plus robuste que foreground_color
                    cell.fill.fore_color.rgb = RGBColor(47, 84, 150)
                    paragraph.font.color.rgb = RGBColor(255, 255, 255)
                    paragraph.font.bold = True
                elif r % 2 == 0:
                    cell.fill.solid()
                    cell.fill.fore_color.rgb = RGBColor(235, 241, 250)

    pptx_path = os.path.join(OUTPUT_DIR, "Deseq2_Rapport_Analyse_Transcriptomique.pptx")
    prs.save(pptx_path)
    print(f"\n> Rapport PowerPoint généré : {pptx_path}")

def process_cohort(cohort_id, metadata_all):
    file_name = COHORT_MAPPING.get(cohort_id)
    if not file_name: return None
    file_path = os.path.join(INPUT_DIR, file_name)
    if not os.path.exists(file_path): return None

    print(f"--- Analyse Cohorte {cohort_id} ---")
    if file_name.endswith('.xlsx'):
        raw_counts = pd.read_excel(file_path, skiprows=1)
    else:
        raw_counts = pd.read_csv(file_path, sep="\t")

    raw_counts.columns = [str(c).strip() for c in raw_counts.columns]
    gene_col = raw_counts.columns[0]
    raw_counts.rename(columns={gene_col: 'Ensembl_Gene_ID'}, inplace=True)
    
    cohort_meta = metadata_all[metadata_all['Cohort'] == cohort_id].copy()
    mapping = {c: c for c in raw_counts.columns if c in cohort_meta['Sample_name'].values}
    
    if not mapping:
        norm_available = {normalize_name(c): c for c in raw_counts.columns if c != 'Ensembl_Gene_ID'}
        for sample in cohort_meta['Sample_name']:
            if normalize_name(sample) in norm_available:
                mapping[norm_available[normalize_name(sample)]] = sample

    raw_counts.rename(columns=mapping, inplace=True)
    raw_counts = raw_counts[['Ensembl_Gene_ID'] + list(mapping.values())]
    raw_counts.set_index('Ensembl_Gene_ID', inplace=True)
    raw_counts = raw_counts.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    cpm_df = calculate_cpm_means(raw_counts, cohort_meta)

    # --- LOGIQUE DE RÉFÉRENCE (Seuil 1%) ---
    ref_grp = 'WT'
    total_genes = len(raw_counts)
    seuil_1_pourcent = int(total_genes * 0.01)

    if 'WT treated' in cohort_meta['Group'].unique():
        res_calib = run_deseq2_logic(raw_counts, cohort_meta, 'WT treated', 'WT')
        if res_calib is not None:
            n_diff = (res_calib['padj'] < 0.05).sum()
            if n_diff < seuil_1_pourcent:
                ref_grp = 'WT treated'
                print(f"   [Info] Utilisation de 'WT treated' comme référence ({n_diff} gènes diff. < {seuil_1_pourcent} [1%])")
            else:
                print(f"   [Info] Maintien de 'WT' comme référence ({n_diff} gènes diff. >= {seuil_1_pourcent} [1%])")

    res_patho = run_deseq2_logic(raw_counts, cohort_meta, 'Disease', ref_grp)
    res_rescue = None
    if 'Rescue' in cohort_meta['Group'].unique():
        res_rescue = run_deseq2_logic(raw_counts, cohort_meta, 'Rescue', 'Disease')

    final = pd.DataFrame(index=raw_counts.index)
    final['Moyenne_CPM_WT'] = cpm_df.get('Moyenne_CPM_WT', 0.0)
    final['Moyenne_CPM_WT_treated'] = cpm_df.get('Moyenne_CPM_WT_treated', 0.0)
    final['Moyenne_CPM_Disease'] = cpm_df.get('Moyenne_CPM_Disease', 0.0)
    final['Moyenne_CPM_Rescue'] = cpm_df.get('Moyenne_CPM_Rescue', 0.0)

    final['Pvalue_Patho'] = res_patho['padj'].reindex(final.index).fillna(1.0) if res_patho is not None else 1.0
    final['Pvalue_Rescue'] = res_rescue['padj'].reindex(final.index).fillna(1.0) if res_rescue is not None else 1.0
    final['Ref_Utilisee'] = ref_grp
    final['Log2FC_Patho'] = res_patho['log2FoldChange'].reindex(final.index).fillna(0.0) if res_patho is not None else 0.0
    final['Log2FC_Rescue'] = res_rescue['log2FoldChange'].reindex(final.index).fillna(0.0) if res_rescue is not None else 0.0

    is_patho = (final['Pvalue_Patho'] < 0.05) & (final['Log2FC_Patho'].abs() > 1)
    final['Significatif_Patho'] = is_patho.map({True: 'OUI', False: 'NON'})
    
    opposed = np.sign(final['Log2FC_Patho']) != np.sign(final['Log2FC_Rescue'])
    is_rescue = (final['Pvalue_Rescue'] < 0.05) & is_patho & opposed
    final['Significatif_Rescue'] = is_rescue.map({True: 'OUI', False: 'NON'})

    final.to_csv(os.path.join(OUTPUT_DIR, f"Deseq2_Genes_Analysis_Cohort_{cohort_id}.csv"), sep=";")
    signature_patho = final[final['Significatif_Patho'] == 'OUI']
    signature_patho.to_csv(os.path.join(OUTPUT_DIR, f"Deseq2_Genes_Signature_Patho_Cohort_{cohort_id}.csv"), sep=";")
    efficacite_rescue = final[final['Significatif_Rescue'] == 'OUI']
    efficacite_rescue.to_csv(os.path.join(OUTPUT_DIR, f"Deseq2_Genes_Efficacite_Rescue_Cohort_{cohort_id}.csv"), sep=";")
    
    print(f"   [OK] 3 fichiers CSV générés pour la Cohorte {cohort_id}")

    return {
        'id': cohort_id,
        'n_wt': len(cohort_meta[cohort_meta['Group'] == 'WT']),
        'n_wtt': len(cohort_meta[cohort_meta['Group'] == 'WT treated']),
        'n_dis': len(cohort_meta[cohort_meta['Group'] == 'Disease']),
        'n_res': len(cohort_meta[cohort_meta['Group'] == 'Rescue']),
        'ref_grp': ref_grp,
        'n_patho_signif': is_patho.sum(),
        'ratio_patho': is_patho.sum() / len(final) if len(final) > 0 else 0,
        'n_rescue_signif': is_rescue.sum(),
        'ratio_rescue': is_rescue.sum() / is_patho.sum() if is_patho.sum() > 0 else 0
    }

def main():
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    all_cohort_stats = []
    meta_data = load_metadata()
    for cid in sorted(COHORT_MAPPING.keys()):
        stats = process_cohort(cid, meta_data)
        if stats: all_cohort_stats.append(stats)
    
    if all_cohort_stats:
        create_global_pptx(all_cohort_stats)

if __name__ == "__main__":
    main()