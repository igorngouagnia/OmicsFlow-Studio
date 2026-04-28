import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import shutil
import gc
from datetime import datetime

# ==========================================
# DESIGN PREMIUM DARK BLUE RESTAURÉ
# ==========================================
st.set_page_config(page_title="OmicsFlow Studio", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
    [data-testid="stAppViewContainer"] { background: linear-gradient(135deg, #0d1117 0%, #161b22 100%); color: #e6edf3; }
    [data-testid="stSidebar"] { background-color: #0d1117; border-right: 1px solid #30363d; }
    h1, h2, h3 { color: #58a6ff !important; font-family: 'Inter', sans-serif; }
    .stButton>button {
        background: linear-gradient(90deg, #1f6feb 0%, #111b27 100%);
        color: white; border: 1px solid #388bfd; border-radius: 8px; font-weight: 600; width: 100%;
    }
    .metric-card { background: rgba(22, 27, 34, 0.8); border: 1px solid #30363d; border-radius: 12px; padding: 20px; text-align: center; }
    .metric-value { font-size: 2rem; font-weight: 800; color: #58a6ff; }
    .green-btn > div > button {
        background: linear-gradient(90deg, #238636 0%, #2ea043 100%) !important;
        color: white !important; border: none !important;
    }
</style>
""", unsafe_allow_html=True)

# --- DIRECTORIES ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")

for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]:
    os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

COHORT_MAPPING = {
    'A': "counts_mtm1_cohort.csv", 'B': "counts_bin1_cohort.csv", 'C': "counts_dnm2_cohort.csv",
    'D': "GSE160078_Raw_gene_counts_matrix_Cohort_MTM1_Updated-07-01-2024.xlsx",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt", 'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

# --- RUNNER ---
def run_pipeline_safe(name, script, category, env_vars=None):
    gc.collect()
    env = os.environ.copy()
    env["OMICS_OUT_DIR"] = os.path.join(SESSION_DIR, category)
    if env_vars: env.update(env_vars)
    with st.spinner(f"🚀 Running {name} Engine..."):
        cmd = ["Rscript", script] if script.endswith(".R") else [sys.executable, script]
        res = subprocess.run(cmd, env=env, capture_output=True, text=True)
        gc.collect()
        if res.returncode == 0:
            st.success(f"✅ {name} Complete.")
            return True
        else:
            st.error(f"❌ {name} Error.")
            with st.expander("Logs"): st.code(res.stderr)
            return False

# ==========================================
# SIDEBAR
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    menu = st.radio("Main Navigation", [
        "📊 Dashboard", 
        "🧪 Transcriptomics Pipeline", 
        "🧪 Proteomics Pipeline", 
        "🧪 Metabolomics Pipeline",
        "🔍 Reference Validations"
    ])
    if st.button("📦 DOWNLOAD ZIP"):
        zip_p = os.path.join(BASE_DIR, "Results.zip")
        shutil.make_archive(zip_p.replace(".zip", ""), 'zip', SESSION_DIR)
        with open(zip_p, "rb") as f: st.download_button("Click to Download ZIP", f, file_name="OmicsFlow_Results.zip")
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD
# ==========================================
if menu == "📊 Dashboard":
    st.title("📊 Analysis Dashboard")
    import plotly.express as px
    view = st.selectbox("View Layer", ["Transcriptomique", "Protéomique", "Métabolomique"])
    cat_dir = os.path.join(SESSION_DIR, view)
    files = [f for f in (os.listdir(cat_dir) if os.path.exists(cat_dir) else []) if "Analysis" in f]
    if files:
        f_sel = st.selectbox("File", files)
        df = pd.read_csv(os.path.join(cat_dir, f_sel))
        lfc = f"Log2FC_{perspective}" if "Log2FC" in df.columns[0] or view=="Transcriptomique" else f"logFC_{perspective}"
        pval = [c for c in df.columns if "Pval" in c or "adj.P" in c][0]
        df['-log10p'] = -np.log10(pd.to_numeric(df[pval], errors='coerce').replace(0, 1e-300))
        fig = px.scatter(df, x=lfc, y='-log10p', color=(df[pval]<0.05), template="plotly_dark")
        st.plotly_chart(fig, use_container_width=True)

# ==========================================
# PIPELINES
# ==========================================
elif "Transcriptomics" in menu:
    st.title("RNA-Seq Engine (DESeq2)")
    cat_dir = os.path.join(SESSION_DIR, "Transcriptomique")
    c1, c2 = st.columns(2)
    with c1:
        meta = st.file_uploader("Metadata (metadata.txt)")
        if meta:
            with open(os.path.join(cat_dir, "metadata.txt"), "wb") as f: f.write(meta.getbuffer())
    with c2:
        cid = st.selectbox("Cohort", ["A","B","C","D","E","F","G"])
        counts = st.file_uploader(f"Counts ({cid})")
        if counts:
            with open(os.path.join(cat_dir, COHORT_MAPPING[cid]), "wb") as f: f.write(counts.getbuffer())
    if st.button("🚀 EXECUTE DESEQ2"):
        run_pipeline_safe("RNA", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py"), "Transcriptomique", {"OMICS_IN_DIR": cat_dir})

elif "Proteomics" in menu:
    st.title("Proteomics Engine")
    cat_dir = os.path.join(SESSION_DIR, "Protéomique")
    p_m = st.file_uploader("metadata.tsv")
    p_g = st.file_uploader("proteinGroups.tsv")
    if p_m:
        with open(os.path.join(cat_dir, "metadata.tsv"), "wb") as f: f.write(p_m.getbuffer())
    if p_g:
        with open(os.path.join(cat_dir, "proteinGroups.tsv"), "wb") as f: f.write(p_g.getbuffer())
    
    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("🚀 RUN PYTHON PIPELINE (Supriya)"):
            run_pipeline_safe("Proteo Python", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya.py"), "Protéomique", {"OMICS_IN_DIR": cat_dir})
    with col_b:
        if st.button("🚀 RUN R PIPELINE (Supriya)"):
            run_pipeline_safe("Proteo R", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya_R.R"), "Protéomique", {"OMICS_IN_DIR": cat_dir})

elif "Metabolomics" in menu:
    st.title("Metabolomics Engine")
    cat_dir = os.path.join(SESSION_DIR, "Métabolomique")
    m_file = st.file_uploader("Upload Metabolomics Data")
    if m_file:
        with open(os.path.join(cat_dir, m_file.name), "wb") as f: f.write(m_file.getbuffer())
    if st.button("🚀 RUN METABOLOMICS ANALYSIS"):
        run_pipeline_safe("Metabo", os.path.join(BASE_DIR, "03_Analyse_métabolomique_Validation.py"), "Métabolomique", {"OMICS_IN_DIR": cat_dir})

# ==========================================
# VALIDATIONS (RNA & PROTEO)
# ==========================================
elif menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    ref = st.file_uploader("Reference Benchmark (Excel)")
    if ref:
        ref_p = os.path.join(v_dir, ref.name)
        with open(ref_p, "wb") as f: f.write(ref.getbuffer())
        
        st.divider()
        st.subheader("A. Transcriptomics (Cohort G)")
        st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
        if st.button("Run RNA Validation"):
             run_pipeline_safe("RNA Val", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique")})
        st.markdown("</div>", unsafe_allow_html=True)
        
        st.divider()
        st.subheader("B. Proteomics (Supriya 2w/7w)")
        st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
        if st.button("Run Protein Validation"):
             run_pipeline_safe("Proteo Val", os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Protéomique")})
        st.markdown("</div>", unsafe_allow_html=True)
        
        # Display results below
        for r in [f for f in os.listdir(v_dir) if f.endswith(".txt")]:
            with st.expander(f"Report: {r}"):
                with open(os.path.join(v_dir, r), "r") as f: st.text(f.read())
