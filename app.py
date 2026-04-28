import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import shutil
import re
from datetime import datetime

# ==========================================
# RESTORATION DU DESIGN "PREMIUM DARK BLUE"
# ==========================================
st.set_page_config(page_title="OmicsFlow Studio", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
    /* Global Background */
    [data-testid="stAppViewContainer"] {
        background: linear-gradient(135deg, #0d1117 0%, #161b22 100%);
        color: #e6edf3;
    }
    [data-testid="stSidebar"] {
        background-color: #0d1117;
        border-right: 1px solid #30363d;
    }
    
    /* Headers */
    h1, h2, h3 {
        color: #58a6ff !important;
        font-family: 'Inter', sans-serif;
        font-weight: 700;
    }

    /* Buttons - Back to the Blue Gradient requested */
    .stButton>button {
        background: linear-gradient(90deg, #1f6feb 0%, #111b27 100%);
        color: white;
        border: 1px solid #388bfd;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-weight: 600;
        width: 100%;
        transition: all 0.3s;
    }
    .stButton>button:hover {
        border-color: #58a6ff;
        box-shadow: 0 0 15px rgba(88,166,255,0.4);
    }

    /* Cards & Containers */
    .metric-card {
        background: rgba(22, 27, 34, 0.8);
        border: 1px solid #30363d;
        border-radius: 12px;
        padding: 25px;
        text-align: center;
    }
    .metric-value { font-size: 2.5rem; font-weight: 800; color: #58a6ff; }
    .metric-label { font-size: 0.9rem; color: #8b949e; text-transform: uppercase; }

    /* Custom classes from original screenshots */
    .status-box {
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
        font-weight: 500;
    }
</style>
""", unsafe_allow_html=True)

# --- CONFIGURATION DES CHEMINS ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state:
    st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")

SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]:
    os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

COHORT_MAPPING = {
    'A': "counts_mtm1_cohort.csv", 'B': "counts_bin1_cohort.csv", 'C': "counts_dnm2_cohort.csv",
    'D': "GSE160078_Raw_gene_counts_matrix_Cohort_MTM1_Updated-07-01-2024.xlsx",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt", 'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

# --- FONCTIONS COEUR ---
def run_script(name, path, category, env_vars=None):
    env = os.environ.copy()
    env["OMICS_OUT_DIR"] = os.path.join(SESSION_DIR, category)
    if env_vars: env.update(env_vars)
    cmd = ["Rscript", path] if path.endswith(".R") else [sys.executable, path]
    
    with st.spinner(f"🚀 Running {name}..."):
        res = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if res.returncode == 0:
            st.success(f"✅ {name} Completed Successfully.")
            return True
        else:
            st.error(f"❌ Error in {name}")
            with st.expander("Show Logs"):
                st.code(res.stderr)
            return False

# ==========================================
# SIDEBAR (Identique aux captures)
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    st.caption(f"Workspace: {st.session_state['session_id']}")
    st.divider()
    
    menu = st.radio("Main Navigation", [
        "📊 Global Analytics Dashboard",
        "🧪 Transcriptomics Pipeline",
        "🧪 Proteomics Pipeline",
        "🧪 Metabolomics Pipeline",
        "🔍 Reference Validations"
    ])
    
    st.divider()
    # Le bouton vert comme dans la capture
    if st.button("📦 DOWNLOAD ENTIRE SESSION (ZIP)", key="dl_zip_btn"):
        zip_path = os.path.join(BASE_DIR, "OmicsFlow_Results.zip")
        shutil.make_archive(zip_path.replace(".zip", ""), 'zip', SESSION_DIR)
        with open(zip_path, "rb") as f:
            st.download_button("Click to Download ZIP", f, file_name=f"OmicsFlow_{st.session_state['session_id']}.zip")

    st.divider()
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD
# ==========================================
if menu == "📊 Global Analytics Dashboard":
    st.title("📊 Multi-Omics Studio")
    st.caption("Central visualization center for your results.")
    
    import plotly.express as px # Import local pour économiser la RAM
    
    layer = st.selectbox("Switch Omics View", ["Transcriptomique", "Protéomique"])
    cat_dir = os.path.join(SESSION_DIR, layer)
    
    files = [f for f in (os.listdir(cat_dir) if os.path.exists(cat_dir) else []) if "Analysis" in f]
    if files:
        f_sel = st.selectbox("Select Results File", files)
        df = pd.read_csv(os.path.join(cat_dir, f_sel), sep=";" if layer=="Transcriptomique" else ",")
        
        lfc = f"Log2FC_{perspective}" if layer=="Transcriptomique" else f"logFC_{perspective}"
        pval = f"Pvalue_{perspective}" if layer=="Transcriptomique" else f"adj.P.Val_{perspective}"
        name = df.columns[0] if layer=="Transcriptomique" else "Unnamed: 2"
        
        c1, c2 = st.columns([3, 1])
        with c1:
            df['-log10p'] = -np.log10(pd.to_numeric(df[pval], errors='coerce').replace(0, 1e-300))
            fig = px.scatter(df, x=lfc, y='-log10p', color=(df[pval] < 0.05),
                             color_discrete_map={True: '#ff4b4b', False: '#4a4e59'},
                             hover_name=name, template="plotly_dark", height=600)
            st.plotly_chart(fig, use_container_width=True)
        
        with c2:
            st.subheader(f"Results ({perspective})")
            sig_count = len(df[df[pval] < 0.05])
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Significatifs</div>
                    <div class='metric-value'>{sig_count}</div>
                </div>
            """, unsafe_allow_html=True)
            
            st.divider()
            st.subheader("Available Files")
            for f in os.listdir(cat_dir):
                with open(os.path.join(cat_dir, f), "rb") as fp:
                    st.download_button(f"📥 {f}", fp, file_name=f, key=f"dash_{f}")

# ==========================================
# TRANSCRIPTOMICS (Pipeline)
# ==========================================
elif menu == "🧪 Transcriptomics Pipeline":
    st.title("RNA-Seq Analysis (DESeq2)")
    cat_dir = os.path.join(SESSION_DIR, "Transcriptomique")
    
    st.subheader("Step 1: Resource Loading")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Metadata Management**")
        meta = st.file_uploader("Upload metadata.txt")
        if meta:
            with open(os.path.join(cat_dir, "metadata.txt"), "wb") as f:
                f.write(meta.getbuffer())
            st.success("Metadata Uploaded.")
            
    with col2:
        st.markdown("**Count Matrix Selection**")
        cid = st.selectbox("Assign to Cohort ID:", ["A","B","C","D","E","F","G"])
        counts = st.file_uploader(f"Upload counts for Cohort {cid}")
        if counts:
            target = COHORT_MAPPING.get(cid, counts.name)
            with open(os.path.join(cat_dir, target), "wb") as f:
                f.write(counts.getbuffer())
            st.success(f"Counts for Cohort {cid} Active.")

    st.divider()
    st.subheader("Step 2: Analysis Execution")
    if st.button("🚀 EXECUTE DESEQ2 WORKFLOW"):
        run_script("RNA", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py"), "Transcriptomique", {"OMICS_IN_DIR": cat_dir})

# ==========================================
# PROTEOMICS (Pipeline)
# ==========================================
elif menu == "🧪 Proteomics Pipeline":
    st.title("Proteomics Differential Expression")
    cat_dir = os.path.join(SESSION_DIR, "Protéomique")
    
    col1, col2 = st.columns(2)
    with col1:
        p_m = st.file_uploader("Upload metadata.tsv")
    with col2:
        p_g = st.file_uploader("Upload proteinGroups.tsv")
        
    if p_m:
        with open(os.path.join(cat_dir, "metadata.tsv"), "wb") as f: f.write(p_m.getbuffer())
    if p_g:
        with open(os.path.join(cat_dir, "proteinGroups.tsv"), "wb") as f: f.write(p_g.getbuffer())
        
    if st.button("🚀 RUN PROTEOMICS PIPELINE"):
        run_script("Proteo Python", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya.py"), "Protéomique", {"OMICS_IN_DIR": cat_dir})
        run_script("Proteo R", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_Supriya_R.R"), "Protéomique", {"OMICS_IN_DIR": cat_dir})

# ==========================================
# VALIDATIONS
# ==========================================
elif menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    
    ref = st.file_uploader("Upload Reference Excel File (mmc2.xlsx / S8)")
    if ref:
        ref_p = os.path.join(v_dir, ref.name)
        with open(ref_p, "wb") as f: f.write(ref.getbuffer())
        st.success("Reference File Active.")
        
        st.divider()
        st.subheader("A. Transcriptomics (Cohort G)")
        if st.button("Run RNA Validation"):
            run_script("RNA Val", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique")})
            
        st.divider()
        st.subheader("B. Proteomics (Table S8)")
        if st.button("Run Protein Validation"):
            run_script("Proteo Val", os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py"), "Validations", {"OMICS_REF_FILE": ref_p, "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Protéomique")})
        
        st.divider()
        for r in [f for f in os.listdir(v_dir) if f.endswith(".txt")]:
            with st.expander(f"View Report: {r}"):
                with open(os.path.join(v_dir, r), "r") as f: st.text(f.read())
