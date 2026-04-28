import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import shutil
import gc
import re
from datetime import datetime
import plotly.express as px

# ==========================================
# DESIGN PREMIUM DARK BLUE
# ==========================================
st.set_page_config(page_title="OmicsFlow Studio", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
    [data-testid="stAppViewContainer"] { background: linear-gradient(135deg, #0d1117 0%, #161b22 100%); color: #e6edf3; }
    [data-testid="stSidebar"] { background-color: #0d1117; border-right: 1px solid #30363d; }
    h1, h2, h3 { color: #58a6ff !important; font-family: 'Inter', sans-serif; }
    .stButton>button { background: linear-gradient(90deg, #1f6feb 0%, #111b27 100%); color: white; border: 1px solid #388bfd; border-radius: 8px; font-weight: 600; width: 100%; }
    .green-btn > div > button { background: linear-gradient(90deg, #238636 0%, #2ea043 100%) !important; color: white !important; border: none !important; }
    .targets-summary { background: rgba(22, 27, 34, 0.5); padding: 20px; border-radius: 12px; border: 1px solid #30363d; }
    .target-val { font-size: 2.5rem; font-weight: 800; color: #58a6ff; margin-bottom: 10px; }
    .target-label { color: #8b949e; font-size: 0.9rem; }
    .metric-box { text-align: left; }
    .metric-label { font-size: 0.85rem; color: #8b949e; }
    .metric-value { font-size: 2.2rem; font-weight: 800; color: #58a6ff; }
    .grayed-out { opacity: 0.3; filter: grayscale(100%); pointer-events: none; }
</style>
""", unsafe_allow_html=True)

# --- DIRECTORIES ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]: os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

COHORT_MAPPING = {
    'A': "counts_mtm1_cohort.csv", 'B': "counts_bin1_cohort.csv", 'C': "counts_dnm2_cohort.csv",
    'D': "GSE160078_Raw_gene_counts_matrix_Cohort_MTM1_Updated-07-01-2024.xlsx",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt", 'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

# --- PARSING FUNCTIONS ---
def get_metrics_rna(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    with open(filepath, "r", encoding="utf-8") as f: content = f.read()
    hits = re.search(r"Intersection_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    missing = re.search(r"Missing_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    extras = re.search(r"Extras_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    h = int(hits.group(1)) if hits else 0
    m = int(missing.group(1)) if missing else 0
    rate = f"{(h / (h + m) * 100):.1f}" if (h + m) > 0 else "0"
    return str(h), str(m), (extras.group(1) if extras else "0"), rate

def get_metrics_proteo(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    with open(filepath, "r", encoding="utf-8") as f: content = f.read()
    ref = re.search(r"Gènes significatifs dans la référence \(Supriya\) : (\d+)", content)
    loc = re.search(r"Gènes significatifs dans notre analyse .* : (\d+)", content)
    inter = re.search(r"Intersection \(Gènes communs\) : (\d+)", content)
    pct = re.search(r"---> Pourcentage d'intersection \(sur réf.\) : ([\d.]+)%", content)
    h = inter.group(1) if inter else "0"
    r_total = int(ref.group(1)) if ref else 0
    l_total = int(loc.group(1)) if loc else 0
    return h, str(r_total - int(h)), str(l_total - int(h)), (pct.group(1) if pct else "0")

# ==========================================
# SIDEBAR
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    menu = st.radio("Main Navigation", ["📊 Dashboard", "🧪 Transcriptomics Pipeline", "🧪 Proteomics Pipeline", "🧪 Metabolomics Pipeline", "🔍 Reference Validations"])
    st.divider()
    if st.button("📦 DOWNLOAD ZIP"):
        zip_p = os.path.join(BASE_DIR, "Results.zip")
        shutil.make_archive(zip_p.replace(".zip", ""), 'zip', SESSION_DIR)
        with open(zip_p, "rb") as f: st.download_button("Download ZIP", f, file_name="OmicsFlow_Results.zip")
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD
# ==========================================
if menu == "📊 Dashboard":
    st.title("📊 Multi-Omics Studio")
    col_search, col_view = st.columns([3, 1])
    with col_search: st.text_input("🔍 Search Genes/Metabolites...", placeholder="Ex: MTM1, BIN1...")
    with col_view: view = st.selectbox("Switch Omics View", ["Transcriptomique", "Protéomique", "Métabolomique"])
    
    cat_dir = os.path.join(SESSION_DIR, view)
    files = [f for f in (os.listdir(cat_dir) if os.path.exists(cat_dir) else []) if f.endswith(".csv")]
    
    c_left, c_right = st.columns([3, 1])
    sig_patho, sig_rescue = 0, 0
    
    with c_left:
        if files:
            f_sel = st.selectbox("Select Cohort File", files)
            df = pd.read_csv(os.path.join(cat_dir, f_sel))
            lfc = [c for c in df.columns if "Log2FC" in c or "logFC" in c][0]
            pval = [c for c in df.columns if "Pval" in c or "adj.P" in c][0]
            df['-log10p'] = -np.log10(pd.to_numeric(df[pval], errors='coerce').replace(0, 1e-300))
            sig_patho = len(df[df[pval] < 0.05])
            sig_rescue = len(df[(df[pval] < 0.05) & (df[lfc].abs() > 1)])
            fig = px.scatter(df, x=lfc, y='-log10p', color=(df[pval]<0.05), template="plotly_dark")
            st.plotly_chart(fig, use_container_width=True)
        else: st.warning("No data found.")
            
    with c_right:
        st.markdown("<div class='targets-summary'>", unsafe_allow_html=True)
        st.subheader("Targets Summary")
        st.markdown(f"<div class='target-label'>Significatifs (Patho)</div><div class='target-val'>{sig_patho if sig_patho else '-'}</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='target-label'>Cibles Rescue Validées</div><div class='target-val'>{sig_rescue if sig_rescue else '-'}</div>", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

# ==========================================
# TRANSCRIPTOMICS PIPELINE
# ==========================================
elif "Transcriptomics" in menu:
    st.title("🧪 RNA-Seq Engine (DESeq2)")
    cat_dir = os.path.join(SESSION_DIR, "Transcriptomique")
    c1, c2 = st.columns(2)
    with c1:
        meta = st.file_uploader("Upload Metadata (metadata.txt)")
        if meta: 
            with open(os.path.join(cat_dir, "metadata.txt"), "wb") as f: f.write(meta.getbuffer())
    with c2:
        cid = st.selectbox("Select Cohort", ["A","B","C","D","E","F","G"])
        counts = st.file_uploader(f"Upload Counts ({cid})")
        if counts: 
            with open(os.path.join(cat_dir, COHORT_MAPPING[cid]), "wb") as f: f.write(counts.getbuffer())
    if st.button("🚀 EXECUTE DESEQ2"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# PROTEOMICS PIPELINE
# ==========================================
elif "Proteomics Pipeline" in menu:
    st.title("🧪 Proteomics Engine")
    cat_dir = os.path.join(SESSION_DIR, "Protéomique")
    pm = st.file_uploader("metadata.tsv")
    pg = st.file_uploader("proteinGroups.tsv")
    if pm: 
        with open(os.path.join(cat_dir, "metadata.tsv"), "wb") as f: f.write(pm.getbuffer())
    if pg: 
        with open(os.path.join(cat_dir, "proteinGroups.tsv"), "wb") as f: f.write(pg.getbuffer())
    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("🚀 RUN PYTHON PIPELINE (Supriya)"):
            subprocess.run([sys.executable, os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})
    with col_b:
        if st.button("🚀 RUN R PIPELINE (Supriya)"):
            subprocess.run(["Rscript", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya_R.R")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# METABOLOMICS PIPELINE
# ==========================================
elif "Metabolomics Pipeline" in menu:
    st.title("🧪 Metabolomics Engine")
    cat_dir = os.path.join(SESSION_DIR, "Métabolomique")
    mf = st.file_uploader("Upload Metabolomics Data")
    if mf: 
        with open(os.path.join(cat_dir, mf.name), "wb") as f: f.write(mf.getbuffer())
    if st.button("🚀 RUN ANALYSIS"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "03_Analyse_métabolomique_Validation.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# REFERENCE VALIDATIONS
# ==========================================
elif menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    val_case = st.selectbox("Select Validation Case", ["Transcriptomics (Cohort G)", "Proteomics (2w)", "Proteomics (7w)"])
    ref_f = st.file_uploader(f"Upload Reference for {val_case}", type=["xlsx"])
    if ref_f: 
        with open(os.path.join(v_dir, "benchmark.xlsx"), "wb") as f: f.write(ref_f.getbuffer())
    
    st.divider()
    
    # Section RNA
    st.subheader("A. Transcriptomics (Cohort G)")
    st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
    if st.button("Run RNA Validation"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py")], env={"OMICS_REF_FILE": os.path.join(v_dir, "benchmark.xlsx"), "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique"), "OMICS_OUT_DIR": v_dir})
    st.markdown("</div>", unsafe_allow_html=True)
    
    rna_h, rna_m, rna_e, rna_r = get_metrics_rna(os.path.join(v_dir, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.txt"))
    cols_rna = st.columns(4)
    cols_rna[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{rna_h}</div></div>", unsafe_allow_html=True)
    cols_rna[3].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{rna_r}%</div></div>", unsafe_allow_html=True)

    st.divider()
    
    # Section Proteo
    st.subheader("B. Proteomics (Table S8)")
    st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
    if st.button("Run Protein Validation"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py")], env={"OMICS_REF_DIR": v_dir, "OMICS_OUT_DIR": os.path.join(SESSION_DIR, "Protéomique")})
    st.markdown("</div>", unsafe_allow_html=True)
    
    p2_h, p2_m, p2_e, p2_r = get_metrics_proteo(os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_2w.txt"))
    p7_h, p7_m, p7_e, p7_r = get_metrics_proteo(os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_7w.txt"))
    
    cp = st.columns(2)
    with cp[0]:
        st.write("**Case: 2 Weeks**")
        cx = st.columns(2)
        cx[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p2_h}</div></div>", unsafe_allow_html=True)
        cx[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p2_r}%</div></div>", unsafe_allow_html=True)
    with cp[1]:
        st.write("**Case: 7 Weeks**")
        cx = st.columns(2)
        cx[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p7_h}</div></div>", unsafe_allow_html=True)
        cx[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p7_r}%</div></div>", unsafe_allow_html=True)
