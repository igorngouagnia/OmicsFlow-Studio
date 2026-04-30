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
    [data-testid="stAppViewContainer"] { background-color: #0d1117; color: #c9d1d9; }
    [data-testid="stSidebar"] { background-color: #0d1117; border-right: 1px solid #30363d; }
    h1, h2 { color: #ffffff !important; }
    h3 { color: #8b949e !important; }
    .stButton>button { background-color: #238636 !important; color: white !important; border-radius: 6px !important; padding: 12px !important; font-weight: 600 !important; width: 100% !important; }
    .success-box { background-color: rgba(35, 134, 54, 0.1); border: 1px solid #238636; padding: 10px; border-radius: 6px; color: #3fb950; margin-bottom: 10px; }
    .metric-box { text-align: left; margin-top: 10px; }
    .metric-label { font-size: 0.85rem; color: #8b949e; }
    .metric-value { font-size: 2.2rem; font-weight: 800; color: #58a6ff; }
    .grayed-out { opacity: 0.2; filter: grayscale(100%); pointer-events: none; }
</style>
""", unsafe_allow_html=True)

# --- DIRECTORIES ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]: os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

COHORT_MAPPING = {'A':"counts_mtm1_cohort.csv",'B':"counts_bin1_cohort.csv",'C':"counts_dnm2_cohort.csv",'D':"GSE160078_Raw_gene_counts_matrix_Cohort_MTM1_Updated-07-01-2024.xlsx",'E':"GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",'F':"GSE282489_raw_counts_bin1_cohort.txt",'G':"GSE282489_raw_counts_dnm2_cohort.txt"}

# --- PARSING ---
def get_metrics_rna(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    try:
        with open(filepath, "r", encoding="utf-8") as f: content = f.read()
        h = re.search(r"Gènes Communs \(Intersection\) : (\d+)", content)
        m = re.search(r"Gènes Manquants \(Papier seul\) : (\d+)", content)
        e = re.search(r"Gènes Extras \(Analyse seule\) : (\d+)", content)
        r = re.search(r"SCORE DE REPRODUCTIBILITÉ : ([\d.]+)%", content)
        
        return (h.group(1) if h else "0"), (m.group(1) if m else "0"), (e.group(1) if e else "0"), (r.group(1) if r else "0")
    except: return "-", "-", "-", "-"

def get_metrics_proteo(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    try:
        with open(filepath, "r", encoding="utf-8") as f: content = f.read()
        it = re.search(r"Intersection \(Gènes communs\) : (\d+)", content)
        rt = re.search(r"---> Pourcentage d'intersection \(sur réf.\) : ([\d.]+)%", content)
        return (it.group(1) if it else "0"), "-", "-", (rt.group(1) if rt else "0")
    except: return "-", "-", "-", "-"

# ==========================================
# SIDEBAR
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    menu = st.radio("Main Navigation", ["📊 Dashboard", "🧪 Transcriptomics Pipeline", "🧪 Proteomics Pipeline", "🧪 Metabolomics Pipeline", "🔍 Reference Validations"])
    st.divider()
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD
# ==========================================
if menu == "📊 Dashboard":
    st.title("📊 Multi-Omics Studio")
    view = st.selectbox("Switch Omics View", ["Transcriptomique", "Protéomique", "Métabolomique"])
    cat_dir = os.path.join(SESSION_DIR, view)
    if view == "Transcriptomique": cat_dir = os.path.join(cat_dir, "Deseq2")
    
    files = [f for f in (os.listdir(cat_dir) if os.path.exists(cat_dir) else []) if f.endswith(".csv")]
    c_l, c_r = st.columns([3, 1])
    with c_l:
        if files:
            f_sel = st.selectbox("Select Cohort File", files)
            df = pd.read_csv(os.path.join(cat_dir, f_sel))
            lfc = [c for c in df.columns if "Log2FC" in c or "logFC" in c][0]
            pval = [c for c in df.columns if "Pval" in c or "adj.P" in c][0]
            df['-log10p'] = -np.log10(pd.to_numeric(df[pval], errors='coerce').replace(0, 1e-300))
            st.plotly_chart(px.scatter(df, x=lfc, y='-log10p', color=(df[pval]<0.05), template="plotly_dark"), use_container_width=True)
    with c_r:
        st.subheader("Targets Summary")
        st.write("Dynamic calculation from selected file.")

# ==========================================
# TRANSCRIPTOMICS (STEP 1 & 2)
# ==========================================
elif menu == "🧪 Transcriptomics Pipeline":
    st.title("RNA-Seq Analysis (DESeq2)")
    cat_dir = os.path.join(SESSION_DIR, "Transcriptomique")
    st.subheader("Step 1: Resource Loading")
    c1, c2 = st.columns(2)
    with c1:
        m = st.file_uploader("Upload metadata.txt", key="rna_m")
        if m: 
            with open(os.path.join(cat_dir, "metadata.txt"), "wb") as f: f.write(m.getbuffer())
            st.markdown("<div class='success-box'>Metadata Uploaded.</div>", unsafe_allow_html=True)
    with c2:
        cid = st.selectbox("Cohort ID:", ["A","B","C","D","E","F","G"])
        c = st.file_uploader(f"Upload counts ({cid})", key="rna_c")
        if c: 
            with open(os.path.join(cat_dir, COHORT_MAPPING[cid]), "wb") as f: f.write(c.getbuffer())
    st.subheader("Step 2: Analysis Execution")
    if st.button("🚀 EXECUTE DESEQ2 WORKFLOW"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# PROTEOMICS (STEP 1 & 2)
# ==========================================
elif menu == "🧪 Proteomics Pipeline":
    st.title("Proteomics Analysis")
    cat_dir = os.path.join(SESSION_DIR, "Protéomique")
    st.subheader("Step 1: Resource Loading")
    c1, c2 = st.columns(2)
    with c1:
        pm = st.file_uploader("Upload metadata.tsv", key="p_m")
        if pm: 
            with open(os.path.join(cat_dir, "metadata.tsv"), "wb") as f: f.write(pm.getbuffer())
    with c2:
        pg = st.file_uploader("Upload proteinGroups.tsv", key="p_g")
        if pg: 
            with open(os.path.join(cat_dir, "proteinGroups.tsv"), "wb") as f: f.write(pg.getbuffer())
    st.subheader("Step 2: Analysis Execution")
    ca, cb = st.columns(2)
    with ca:
        if st.button("🚀 RUN PYTHON PIPELINE"):
            subprocess.run([sys.executable, os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})
    with cb:
        if st.button("🚀 RUN R PIPELINE"):
            subprocess.run(["Rscript", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya_R.R")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# METABOLOMICS (STEP 1 & 2)
# ==========================================
elif menu == "🧪 Metabolomics Pipeline":
    st.title("🧪 Metabolomics Analysis")
    cat_dir = os.path.join(SESSION_DIR, "Métabolomique")
    st.subheader("Step 1: Resource Loading")
    mf = st.file_uploader("Upload Metabolomics Data (CSV/XLSX)")
    if mf: 
        with open(os.path.join(cat_dir, mf.name), "wb") as f: f.write(mf.getbuffer())
        st.success("Data Uploaded.")
    st.subheader("Step 2: Analysis Execution")
    if st.button("🚀 EXECUTE METABOLOMICS WORKFLOW"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "03_Analyse_métabolomique_Validation.py")], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

# ==========================================
# REFERENCE VALIDATIONS
# ==========================================
elif menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    val_case = st.selectbox("Select Validation Case", ["Transcriptomics (Cohort G)", "Proteomics (2w)", "Proteomics (7w)"])
    st.caption("Expected: brain-2021-02002-File011.xlsx (RNA), results_2w.xlsx or results_7w.xlsx (Proteo)")
    ref_f = st.file_uploader("Upload Reference File", type=["xlsx"])
    mode = "none"
    if ref_f:
        if "brain-2021-02002-File011" in ref_f.name: mode = "rna"
        elif "results_2w" in ref_f.name: mode = "p2w"
        elif "results_7w" in ref_f.name: mode = "p7w"
        with open(os.path.join(v_dir, ref_f.name), "wb") as f: f.write(ref_f.getbuffer())
    st.divider()
    st.markdown(f"<div class='{'active' if mode in ['rna','none'] else 'grayed-out'}'>", unsafe_allow_html=True)
    st.subheader("A. Transcriptomics (Cohort G)")
    if st.button("Run RNA Validation"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py")], env={"OMICS_REF_FILE": os.path.join(v_dir, ref_f.name), "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique"), "OMICS_OUT_DIR": v_dir})
    rh, rm, re, rr = get_metrics_rna(os.path.join(v_dir, "Deseq2", "01_Analyse_transcriptomique_Comparaison_G_deseq2.txt"))
    c_rna = st.columns(4)
    c_rna[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{rh}</div></div>", unsafe_allow_html=True)
    c_rna[3].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{rr}%</div></div>", unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)
    st.divider()
    st.markdown(f"<div class='{'active' if mode in ['p2w','p7w','none'] else 'grayed-out'}'>", unsafe_allow_html=True)
    st.subheader("B. Proteomics (Table S8)")
    if st.button("Run Protein Validation"):
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py")], env={"OMICS_REF_DIR": v_dir, "OMICS_OUT_DIR": os.path.join(SESSION_DIR, "Protéomique")})
    p2h, p2m, p2e, p2r = get_metrics_proteo(os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_2w.txt"))
    p7h, p7m, p7e, p7r = get_metrics_proteo(os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_7w.txt"))
    cp = st.columns(2)
    with cp[0]:
        st.markdown(f"<div class='{'active' if mode in ['p2w','none'] else 'grayed-out'}'>", unsafe_allow_html=True)
        st.write("**Case: 2 Weeks**")
        cx = st.columns(2)
        cx[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p2h}</div></div>", unsafe_allow_html=True)
        cx[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p2r}%</div></div>", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
    with cp[1]:
        st.markdown(f"<div class='{'active' if mode in ['p7w','none'] else 'grayed-out'}'>", unsafe_allow_html=True)
        st.write("**Case: 7 Weeks**")
        cx = st.columns(2)
        cx[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p7h}</div></div>", unsafe_allow_html=True)
        cx[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p7r}%</div></div>", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)
