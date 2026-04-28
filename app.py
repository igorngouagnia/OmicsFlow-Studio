import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import os
import subprocess
import sys
import re
import shutil
from datetime import datetime
import gc

# ==========================================
# MEMORY-OPTIMIZED CONFIG v3.1
# ==========================================
st.set_page_config(page_title="OmicsFlow Studio", layout="wide", initial_sidebar_state="expanded")

# --- UTILS ---
@st.cache_data(ttl=300, max_entries=5)
def load_omics_lite(file_path, cols_to_keep):
    if not os.path.exists(file_path): return None
    try:
        # Optimization: only read required columns
        df = pd.read_csv(file_path, sep=";", usecols=lambda x: x in cols_to_keep or x == "index")
        if len(df.columns) < 2: 
            df = pd.read_csv(file_path, sep=",", usecols=lambda x: x in cols_to_keep or x == "index")
        gc.collect()
        return df
    except:
        try:
            df = pd.read_csv(file_path, sep=";") if ";" in open(file_path).read(100) else pd.read_csv(file_path, sep=",")
            return df[df.columns.intersection(cols_to_keep)]
        except: return None

def parse_validation_summary(file_path):
    if not os.path.exists(file_path): return None
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            count = re.search(r"Nombre de (?:gènes|protéines) (?:significatifs|identifiés) : (\d+)", content)
            rate = re.search(r"Taux de (?:recouvrement|matching) : ([\d\.]+)", content)
            return {"count": count.group(1) if count else "N/A", "rate": rate.group(1) if rate else "0"}
    except: return None

# --- STATE ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")

COHORT_MAPPING = {
    'A': "counts_mtm1_cohort.csv", 'B': "counts_bin1_cohort.csv", 'C': "counts_dnm2_cohort.csv",
    'D': "GSE160078_Raw_gene_counts_matrix_Cohort_MTM1_Updated-07-01-2024.xlsx",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt", 'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]:
    os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

st.markdown("""
<style>
    [data-testid="stAppViewContainer"] { background: linear-gradient(135deg, #0e1117 0%, #161b22 100%); color: #e6edf3; }
    [data-testid="stSidebar"] { background-color: #0d1117; border-right: 1px solid #30363d; }
    .metric-card { background: rgba(22, 27, 34, 0.7); border: 1px solid #30363d; border-radius: 10px; padding: 15px; text-align: center; }
    .metric-value { font-size: 2.2rem; font-weight: 800; color: #58a6ff; }
    .stButton>button { background: linear-gradient(90deg, #238636 0%, #2ea043 100%); color: white; border: none; }
</style>
""", unsafe_allow_html=True)

def run_script(name, path, category, env_vars=None):
    env = os.environ.copy()
    env["OMICS_OUT_DIR"] = os.path.join(SESSION_DIR, category)
    if env_vars: env.update(env_vars)
    cmd = ["Rscript", path] if path.endswith(".R") else [sys.executable, path]
    with st.spinner(f"🚀 Running {name}..."):
        res = subprocess.run(cmd, env=env, capture_output=True, text=True)
        gc.collect()
        if res.returncode == 0:
            st.success(f"✅ {name} Success")
            return True
        else:
            st.error(f"❌ {name} Error")
            with st.expander("Show Logs"): st.code(res.stderr)
            return False

# ==========================================
# SIDEBAR
# ==========================================
with st.sidebar:
    st.title("OmicsFlow Studio")
    st.caption(f"Workspace: {st.session_state['session_id']}")
    st.divider()
    menu = st.radio("Navigation", ["📊 Dashboard", "🧪 Pipelines", "🔍 Validations"])
    st.divider()
    if st.button("📦 ZIP ENTIRE SESSION"):
        zip_p = os.path.join(BASE_DIR, "Session_Results.zip")
        shutil.make_archive(zip_p.replace(".zip", ""), 'zip', SESSION_DIR)
        with open(zip_p, "rb") as f: st.download_button("Download ZIP", f, file_name="OmicsFlow_Results.zip")
    st.divider()
    perspective = st.selectbox("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD
# ==========================================
if menu == "📊 Dashboard":
    st.title("📊 Multi-Omics Studio")
    view = st.selectbox("Omics Layer", ["Transcriptomique", "Protéomique", "Métabolomique"])
    df_plot = None
    
    if view == "Transcriptomique":
        t_dir = os.path.join(SESSION_DIR, "Transcriptomique")
        files = [f for f in (os.listdir(t_dir) if os.path.exists(t_dir) else []) if "Analysis" in f]
        if files:
            f_sel = st.selectbox("Select Cohort", files)
            lfc, pval = f"Log2FC_{perspective}", f"Pvalue_{perspective}"
            df_plot = load_omics_lite(os.path.join(t_dir, f_sel), [lfc, pval])
            name_col = df_plot.columns[0] if df_plot is not None else ""

    elif view == "Protéomique":
        p_dir = os.path.join(SESSION_DIR, "Protéomique")
        files = [f for f in (os.listdir(p_dir) if os.path.exists(p_dir) else []) if "Analysis" in f]
        if files:
            f_sel = st.selectbox("Select Results", files)
            lfc, pval = f"logFC_{perspective}", f"adj.P.Val_{perspective}"
            df_plot = load_omics_lite(os.path.join(p_dir, f_sel), ["Unnamed: 2", lfc, pval])
            name_col = "Unnamed: 2"

    if df_plot is not None:
        c1, c2 = st.columns([3, 1])
        with c1:
            df_v = df_plot.copy()
            df_v['-log10p'] = -np.log10(pd.to_numeric(df_v[pval], errors='coerce').replace(0, 1e-300))
            df_v['Sig'] = 'None'
            df_v.loc[(df_v[lfc] > 0.5) & (df_v[pval] < 0.05), 'Sig'] = 'Up'
            df_v.loc[(df_v[lfc] < -0.5) & (df_v[pval] < 0.05), 'Sig'] = 'Down'
            
            fig = px.scatter(df_v, x=lfc, y='-log10p', color='Sig', 
                             color_discrete_map={'Up':'#ff4b4b', 'Down':'#1f77b4', 'None':'#4a4e59'},
                             hover_name=name_col, template="plotly_dark", height=600)
            st.plotly_chart(fig, use_container_width=True)
            del df_v; gc.collect()
            
        with c2:
            st.subheader("Summary")
            st.markdown(f"<div class='metric-card'><div class='metric-label'>Significatifs</div><div class='metric-value'>{len(df_plot[df_plot[pval]<0.05])}</div></div>", unsafe_allow_html=True)
            st.divider()
            for f in os.listdir(os.path.join(SESSION_DIR, view)):
                with open(os.path.join(SESSION_DIR, view, f), "rb") as fp:
                    st.download_button(f"📥 {f}", fp, file_name=f, key=f"dash_{f}")

# ==========================================
# PIPELINES
# ==========================================
elif menu == "🧪 Pipelines":
    st.title("🧪 Analysis Pipelines")
    layer = st.selectbox("Select Engine", ["Transcriptomique", "Protéomique"])
    cat_dir = os.path.join(SESSION_DIR, layer)
    
    if layer == "Transcriptomique":
        st.subheader("Step 1: Load Data")
        c1, c2 = st.columns(2)
        with c1:
            meta = st.file_uploader("Metadata (metadata.txt)")
            if meta:
                with open(os.path.join(cat_dir, "metadata.txt"), "wb") as f:
                    f.write(meta.getbuffer())
        with c2:
            cid = st.selectbox("Cohort", ["A","B","C","D","E","F","G"])
            counts = st.file_uploader(f"Counts for Cohort {cid}")
            if counts:
                with open(os.path.join(cat_dir, COHORT_MAPPING[cid]), "wb") as f:
                    f.write(counts.getbuffer())
        
        st.subheader("Step 2: Execution")
        if st.button("🚀 EXECUTE DESEQ2 WORKFLOW"):
            run_script("RNA", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py"), layer, {"OMICS_IN_DIR": cat_dir})

    elif layer == "Protéomique":
        st.subheader("Step 1: Load Data")
        c1, c2 = st.columns(2)
        with c1:
            p_m = st.file_uploader("Metadata.tsv")
            if p_m:
                with open(os.path.join(cat_dir, "metadata.tsv"), "wb") as f:
                    f.write(p_m.getbuffer())
        with c2:
            p_g = st.file_uploader("proteinGroups.tsv")
            if p_g:
                with open(os.path.join(cat_dir, "proteinGroups.tsv"), "wb") as f:
                    f.write(p_g.getbuffer())
        
        st.subheader("Step 2: Execution")
        if st.button("🚀 RUN PROTEOMICS PIPELINE"):
            run_script("Proteo Python", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya.py"), layer, {"OMICS_IN_DIR": cat_dir})
            run_script("Proteo R", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya_R.R"), layer, {"OMICS_IN_DIR": cat_dir})

# ==========================================
# VALIDATIONS
# ==========================================
elif menu == "🔍 Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    ref = st.file_uploader("Reference Benchmark File (Excel)")
    
    if ref:
        ref_p = os.path.join(v_dir, ref.name)
        with open(ref_p, "wb") as f: f.write(ref.getbuffer())
        
        st.divider()
        st.subheader("A. Transcriptomics (Cohort G)")
        if st.button("Run RNA Validation"):
            run_script("RNA Val", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique")})
        
        st.divider()
        st.subheader("B. Proteomics (Table S8)")
        if st.button("Run Proteo Validation"):
            run_script("Proteo Val", os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Protéomique")})
        
        st.divider()
        for r in [f for f in os.listdir(v_dir) if f.endswith(".txt")]:
            met = parse_validation_summary(os.path.join(v_dir, r))
            if met:
                st.write(f"**Report: {r}**")
                col1, col2 = st.columns(2)
                col1.metric("Common Hits", met['count'])
                col2.metric("Matching Rate", f"{met['rate']}%")
                with st.expander("Details"):
                    with open(os.path.join(v_dir, r), "r") as f: st.text(f.read())
