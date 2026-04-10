import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import os
import subprocess
import shutil
import tempfile
import sys
from datetime import datetime
import gc

# ==========================================
# CONFIGURATION & PREMIUM DESIGN (TRANSCRIPTOFLOW CLONE)
# ==========================================
st.set_page_config(page_title="OmicsFlow Studio", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
    /* Global Background & Text */
    .stApp { background-color: #0d1117; color: #c9d1d9; }
    [data-testid="stSidebar"] { background-color: #161b22; border-right: 1px solid #30363d; }
    
    /* Metrics / Cards Styling */
    div[data-testid="metric-container"] {
        background-color: #161b22;
        border: 1px solid #30363d;
        border-radius: 12px;
        padding: 20px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.5);
    }
    label[data-testid="stMetricLabel"] { color: #8b949e !important; font-size: 14px !important; }
    div[data-testid="stMetricValue"] { color: #58a6ff !important; font-weight: bold !important; }
    
    /* Inputs & Selectors */
    .stSelectbox>div>div, .stNumberInput>div>div, .stFileUploader {
        background-color: #21262d !important;
        color: #c9d1d9 !important;
        border-radius: 8px !important;
    }
    
    h1, h2, h3 { color: #58a6ff; font-weight: 600; }
    
    /* Buttons Styling */
    .stButton>button {
        background: linear-gradient(135deg, #238636 0%, #2ea043 100%);
        color: white;
        font-weight: bold;
        border-radius: 8px;
        border: none;
        padding: 10px 24px;
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(46, 160, 67, 0.4);
        background: linear-gradient(135deg, #2ea043 0%, #3fb950 100%);
    }
    
    /* Warning/Info Tweak */
    .stAlert { background-color: rgba(22, 27, 34, 0.8) !important; border: 1px solid #30363d !important; }
</style>
""", unsafe_allow_html=True)

# Paths - Fully Portable for Cloud (Relative to App Root)
BASE_DIR = os.getcwd() 
APP_DIR = BASE_DIR # The app now lives in the repo root for cloud deployment

COHORT_MAPPING = {
    'A': "GSE160079_Raw_gene_counts_matrix_Cohort_MTM1-a.txt",
    'B': "GSE160081_Raw_gene_counts_matrix_Cohort_MTM1-b.txt",
    'C': "GSE160083_Raw_gene_counts_matrix_Cohort_MTM1-c.txt",
    'D': "GSE160077_Raw_gene_counts_matrix_Cohort_BIN1.txt",
    'E': "GSE160078_Raw_gene_counts_matrix_Cohort_DNM2_Updated-07-01-2024.xlsx",
    'F': "GSE282489_raw_counts_bin1_cohort.txt",
    'G': "GSE282489_raw_counts_dnm2_cohort.txt"
}

# Session Management
session_defaults = {
    'session_id': datetime.now().strftime("%Y%m%d_%H%M%S"),
    'transcripto_done': False,
    'proteo_py_done': False,
    'proteo_r_done': False,
    'metabolo_done': False,
    'val_rna_done': False,
    'val_prot_done': False
}
for key, val in session_defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# For cloud deployment, we store sessions in a temporary subdirectory of the app
SESSION_DIR = os.path.join(APP_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
SUB_DIRS = ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]
for sub in SUB_DIRS:
    os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

# Run Engine
def run_script(script_name, script_path, category, extra_env=None):
    env = os.environ.copy()
    out_dir = os.path.join(SESSION_DIR, category)
    env["OMICS_OUT_DIR"] = out_dir
    if extra_env: env.update(extra_env)
        
    cmd = ["Rscript", script_path] if script_path.endswith(".R") else [sys.executable, script_path]
    
    with st.spinner(f"🚀 Processing {script_name}..."):
        try:
            res = subprocess.run(cmd, env=env, capture_output=True, text=True)
            if res.returncode == 0:
                st.success(f"{script_name} Completed.")
                # Auto-Volcano
                env_v = {"OMICS_OUT_DIR": SESSION_DIR}
                subprocess.run([sys.executable, os.path.join(BASE_DIR, "04_Volcano_Plots_Generator.py")], env=env_v)
                return True
            else:
                st.error(f"Error in {script_name}")
                with st.expander("Logs (Full Output)"):
                    st.markdown("**Standard Error:**")
                    st.code(res.stderr if res.stderr else "No error output.")
                    st.markdown("**Standard Output:**")
                    st.code(res.stdout if res.stdout else "No standard output.")
                return False
        except Exception as e:
            st.error(str(e))
            return False

# ==========================================
# SIDEBAR
# ==========================================
st.sidebar.markdown("<h2 style='text-align: center; color: #58a6ff;'>🧬 OmicsFlow</h2>", unsafe_allow_html=True)
st.sidebar.caption(f"<div style='text-align: center;'>Workspace: {st.session_state['session_id']}</div><br>", unsafe_allow_html=True)

menu = st.sidebar.radio("Main Navigation", [
    "📊 Global Analytics Dashboard",
    "🧪 Transcriptomics Pipeline", 
    "🧪 Proteomics Pipeline",
    "🧪 Metabolomics Pipeline",
    "🔍 Reference Validations"
])

st.sidebar.markdown("---")
st.sidebar.info("Application built for research reproduction. All results are isolated in the session folder.")

# --- CACHED DATA LOADERS ---
@st.cache_data(ttl=3600)
def load_omics_data(file_path):
    if not os.path.exists(file_path): return None
    return pd.read_csv(file_path, sep=";")

@st.cache_data(ttl=3600)
def load_validation_data(file_path):
    if not os.path.exists(file_path): return None
    return pd.read_csv(file_path, sep=";")

# ==========================================
# DASHBOARD GLOBAL
# ==========================================
if menu == "📊 Global Analytics Dashboard":
    st.title("📊 Multi-Omics Studio")
    st.markdown("Central visualization center for your results.")
    
    c_srch, c_sel = st.columns([2, 1])
    with c_srch: search_q = st.text_input("🔍 Search Genes/Metabolites...", placeholder="Ex: MTM1, BIN1...")
    with c_sel: layer = st.selectbox("Switch Omics View", ["Transcriptomique", "Protéomique", "Métabolomique"])

    persp_idx = 0 if menu == "🔍 Reference Validations" else (1 if st.session_state.get('perspective') == "Rescue" else 0)
    perspective = st.sidebar.radio("Perspective", ["Patho", "Rescue"], index=persp_idx, help="Switch between Disease effect (Patho) and Treatment effect (Rescue)")
    st.session_state['perspective'] = perspective
    
    col_plot, col_stats = st.columns([2.5, 1])
    
    df_plot = None
    if layer == "Transcriptomique":
        t_dir = os.path.join(SESSION_DIR, "Transcriptomique")
        csv_dir_files = os.listdir(t_dir) if os.path.exists(t_dir) else []
        csvs = [f for f in csv_dir_files if f.startswith("Deseq2_Genes_Analysis")] 
        if csvs:
            c = st.sidebar.selectbox("Select Cohort File", csvs)
            df_plot = load_omics_data(os.path.join(t_dir, c))
            if df_plot is not None:
                name_col = df_plot.index.name if df_plot.index.name else df_plot.columns[0]
                lfc_col, pval_col = f'Log2FC_{perspective}', f'Pvalue_{perspective}'
                lfc_th = 1.0

    elif layer == "Protéomique":
        p_dir = os.path.join(SESSION_DIR, "Protéomique")
        path_py = os.path.join(p_dir, "Proteomics_Analysis_Full_MTM1a.csv")
        path_r = os.path.join(p_dir, "Proteomics_Analysis_Full_MTM1.csv")
        opts = []
        if os.path.exists(path_py): opts.append("Python Engine (Welch)")
        if os.path.exists(path_r): opts.append("R Engine (DEP/Limma)")
        
        if opts:
            engine = st.sidebar.selectbox("Engine Output", opts)
            age = st.sidebar.selectbox("Target Age", ["E18.5", "2w", "7w"])
            used_path = path_py if "Python" in engine else path_r
            df_plot = load_omics_data(used_path)
            if df_plot is not None:
                name_col = 'Gene.names' if 'Gene.names' in df_plot.columns else 'Gene names'
                lfc_col = f'Log2FC_{perspective}_{age}'
                pval_col = f'Padj_{perspective}_{age}' if f'Padj_{perspective}_{age}' in df_plot.columns else f'Pvalue_{perspective}_{age}'
                lfc_th = 1.0

    elif layer == "Métabolomique":
        m_dir = os.path.join(SESSION_DIR, "Métabolomique")
        m_path = os.path.join(m_dir, "Metabolomics_Full_Analysis.csv")
        if os.path.exists(m_path):
            df_plot = load_omics_data(m_path)
            if df_plot is not None:
                name_col = 'CHEMICAL_NAME'
                lfc_col, pval_col = f'Log2FC_{perspective}', f'Padj_{perspective}'
                lfc_th = 0.0

    if df_plot is not None and not df_plot.empty and lfc_col in df_plot.columns:
        df_volc = df_plot.copy()
        if df_volc.index.name: df_volc = df_volc.reset_index()
            
        df_volc['-log10(p-value)'] = -np.log10(pd.to_numeric(df_volc[pval_col], errors='coerce').replace(0, 1e-300))
        conds = [(df_volc[pval_col] < 0.05) & (df_volc[lfc_col] > lfc_th), (df_volc[pval_col] < 0.05) & (df_volc[lfc_col] < -lfc_th)]
        df_volc['Sig'] = np.select(conds, ['Up-Regulated', 'Down-Regulated'], 'Not Significant')
        if search_query := search_q: df_volc = df_volc[df_volc[name_col].str.contains(search_query, case=False, na=False)]
        
        with col_plot:
            st.markdown(f"#### View: {perspective} ({layer})")
            # Minimal subset for plotting to save memory
            df_plot_mini = df_volc[[name_col, lfc_col, pval_col, '-log10(p-value)', 'Sig']].copy()
            fig = px.scatter(df_plot_mini, x=lfc_col, y='-log10(p-value)', color='Sig',
                             color_discrete_map={'Up-Regulated': '#ff4b4b', 'Down-Regulated': '#1f77b4', 'Not Significant': '#4a4e59'},
                             hover_name=name_col)
            fig.update_layout(template="plotly_dark", plot_bgcolor='#0d1117', paper_bgcolor='#0d1117', margin=dict(t=10, l=0, r=0, b=0))
            st.plotly_chart(fig, use_container_width=True)
            del df_plot_mini
            gc.collect()
            
        with col_stats:
            st.markdown(f"### Perspective: {perspective}")
            st.metric("Total Entities", len(df_plot), help="Total number of genes/proteins/metabolites in this dataset")
            
            n_sig = len(df_volc[df_volc['Sig'] != 'Not Significant'])
            st.metric("Significant Hits", n_sig, help=f"Entities with P < 0.05 and |Log2FC| > {lfc_th} for {perspective}")
    else:
        st.info("⚠️ No results found.")

# ==========================================
# TRANSCRIPTOMICS
# ==========================================
elif menu == "🧪 Transcriptomics Pipeline":
    st.title("RNA-Seq Analysis (DESeq2)")
    t_in = os.path.join(SESSION_DIR, "Transcriptomique")
    
    st.subheader("Step 1: Resource Loading")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("#### Metadata Management")
        meta = st.file_uploader("Upload metadata.txt", type=["txt", "tsv"])
        if meta:
            with open(os.path.join(t_in, "metadata.txt"), "wb") as f: f.write(meta.getbuffer())
            st.success("Metadata Uploaded.")
    with c2:
        st.markdown("#### Count Matrix Selection")
        coh = st.selectbox("Assign to Cohort ID:", list(COHORT_MAPPING.keys()))
        counts = st.file_uploader(f"Upload counts for Cohort {coh}", type=["txt", "csv", "tsv", "xlsx"])
        if counts:
            with open(os.path.join(t_in, COHORT_MAPPING[coh]), "wb") as f: f.write(counts.getbuffer())
            st.success(f"File Assigned to {COHORT_MAPPING[coh]}.")

    st.subheader("Step 2: Analysis Execution")
    ready = os.path.exists(os.path.join(t_in, "metadata.txt")) and any(f in COHORT_MAPPING.values() for f in os.listdir(t_in))
    if ready:
        if st.button("🚀 EXECUTE DESEQ2 WORKFLOW", use_container_width=True):
            if run_script("Transcriptomics", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Validation_Deseq2.py"), "Transcriptomique", {"OMICS_IN_DIR": t_in}):
                st.session_state['transcripto_done'] = True
    else:
        st.error("🔒 Pipeline Locked: Metadata AND at least one Cohort count file are required.")

    if st.session_state['transcripto_done']:
        st.info("✅ Latest analysis results are ready in the Dashboard.")

# ==========================================
# PROTEOMICS
# ==========================================
elif menu == "🧪 Proteomics Pipeline":
    st.title("Proteomics Analysis (MaxQuant)")
    p_in = os.path.join(SESSION_DIR, "Protéomique")
    
    st.subheader("Step 1: Resource Loading")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("#### Experimental Design")
        p_meta = st.file_uploader("Upload metadata.tsv", type=["tsv", "txt"])
        if p_meta:
            with open(os.path.join(p_in, "metadata.tsv"), "wb") as f: f.write(p_meta.getbuffer())
            st.success("Proteomics Metadata Uploaded.")
    with c2:
        st.markdown("#### Protein Intensity Table")
        p_grp = st.file_uploader("Upload proteinGroups.tsv", type=["tsv", "txt"])
        if p_grp:
            with open(os.path.join(p_in, "proteinGroups.tsv"), "wb") as f: f.write(p_grp.getbuffer())
            st.success("proteinGroups Uploaded.")

    st.subheader("Step 2: Statistics Engine Selection")
    if os.path.exists(os.path.join(p_in, "metadata.tsv")) and os.path.exists(os.path.join(p_in, "proteinGroups.tsv")):
        bt1, bt2 = st.columns(2)
        with bt1:
            if st.button("🐍 RUN PYTHON (Welch/Pingouin)", use_container_width=True):
                if run_script("Proteo Python", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation.py"), "Protéomique", {"OMICS_IN_DIR": p_in}):
                    st.session_state['proteo_py_done'] = True
        with bt2:
            if st.button("🧊 RUN R (DEP/Limma)", use_container_width=True):
                if run_script("Proteo R", os.path.join(BASE_DIR, "02_Analyse_protéomique_Validation_R.R"), "Protéomique", {"OMICS_IN_DIR": p_in}):
                    st.session_state['proteo_r_done'] = True
    else:
        st.error("🔒 Pipeline Locked: Both metadata.tsv and proteinGroups.tsv are required.")

    if st.session_state['proteo_py_done'] or st.session_state['proteo_r_done']:
        st.info("✅ Proteomics results are stored in the session folder.")

# ==========================================
# METABOLOMICS
# ==========================================
elif menu == "🧪 Metabolomics Pipeline":
    st.title("Metabolomics Analysis (IGBM)")
    m_in = os.path.join(SESSION_DIR, "Métabolomique")
    
    st.subheader("Step 1: Unified Data Upload")
    m_xlsx = st.file_uploader("Upload IGBM Excel Table", type=["xlsx"])
    if m_xlsx:
        with open(os.path.join(m_in, "IGBM-01-21VW MUSCLE DATA TABLES.XLSX"), "wb") as f: f.write(m_xlsx.getbuffer())
        st.success("Metabolomics Data Uploaded.")
        
    st.subheader("Step 2: Analysis Execution")
    if os.path.exists(os.path.join(m_in, "IGBM-01-21VW MUSCLE DATA TABLES.XLSX")):
        if st.button("🧬 EXECUTE METABOLOMICS ENGINE", use_container_width=True):
            if run_script("Metabolomics", os.path.join(BASE_DIR, "03_Analyse_métabolomique_Validation.py"), "Métabolomique", {"OMICS_IN_DIR": m_in}):
                st.session_state['metabolo_done'] = True
    else:
        st.error("🔒 Pipeline Locked: IGBM Excel file is required.")

    if st.session_state['metabolo_done']:
        st.info("✅ Metabolomics analysis complete.")

# ==========================================
# VALIDATIONS (EXPLICIT PATHO/RESCUE)
# ==========================================
elif menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    st.markdown("Compare results with published findings in the original paper.")
    
    v_in = os.path.join(SESSION_DIR, "Validations")
    ref = st.file_uploader("📁 Upload Reference Excel File (mmc2.xlsx / S8)", type=["xlsx"])
    ref_p = None
    if ref:
        ref_p = os.path.join(v_in, ref.name)
        with open(ref_p, "wb") as f: f.write(ref.getbuffer())
        st.success("Reference File Active.")

    st.markdown("---")
    
    st.subheader("A. Transcriptomics (Cohort G)")
    if st.button("Run RNA Validation"):
        if not ref_p: st.warning("Upload Reference First.")
        else:
            if run_script("ARN Validation", os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py"), "Validations", {"OMICS_REF_FILE": ref_p}):
                st.session_state['val_rna_done'] = True
    
    v_files = os.listdir(v_in) if os.path.exists(v_in) else []
    tr_hits = [f for f in v_files if f.startswith("Intersection_RNA_Patho")]
    if tr_hits:
        df_tr = load_validation_data(os.path.join(v_in, tr_hits[0]))
        miss_tr = [f for f in v_files if f.startswith("Missing_RNA_Patho")]
        extra_tr = [f for f in v_files if f.startswith("Extras_RNA_Patho")]
        
        n_common = len(df_tr) if df_tr is not None else 0
        n_miss = len(load_validation_data(os.path.join(v_in, miss_tr[0]))) if miss_tr else 0
        n_extra = len(load_validation_data(os.path.join(v_in, extra_tr[0]))) if extra_tr else 0
        
        total_ref = n_common + n_miss
        total_ana = n_common + n_extra
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Common Hits", n_common, help="Found in both paper and your analysis")
        c2.metric("Missing (Paper)", n_miss, help="In paper but not in your analysis")
        c3.metric("Extras (Analysis)", n_extra, help="Novel hits in your analysis only")
        c4.metric("Novelty Rate", f"{(n_extra/total_ana*100):.1f}%" if total_ana > 0 else "0%")
        
        st.markdown(f"**Sensitivity (Recall):** `{(n_common/total_ref*100):.1f}%` of published genes recovered.")
        with st.expander("Show RNA Overlap"): st.dataframe(df_tr)

    st.markdown("---")
    st.subheader("B. Proteomics (Table S8)")
    if st.button("Run Protein Validation"):
        if not ref_p: st.warning("Upload Reference First.")
        else:
            if run_script("Proteo Validation", os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_papier_Cohorte_A.py"), "Validations", {"OMICS_REF_FILE": ref_p}):
                st.session_state['val_prot_done'] = True
    
    p_hits = [f for f in v_files if f.startswith("Intersection_Protein_Patho")]
    if p_hits:
        # Use latest selection or E18.5/2w/7w fallback
        df_i = load_validation_data(os.path.join(v_in, p_hits[-1]))
        miss_p = [f for f in v_files if f.startswith("Missing_Protein_Patho")]
        extra_p = [f for f in v_files if f.startswith("Extras_Protein_Patho")]
        
        n_c = len(df_i) if df_i is not None else 0
        n_m = len(load_validation_data(os.path.join(v_in, miss_p[-1]))) if miss_p else 0
        n_e = len(load_validation_data(os.path.join(v_in, extra_p[-1]))) if extra_p else 0
        
        total_p_ref = n_c + n_m
        total_p_ana = n_c + n_e
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Common Hits", n_c)
        c2.metric("Missing (Paper)", n_m)
        c3.metric("Extras (Analysis)", n_e)
        c4.metric("Novelty Rate", f"{(n_e/total_p_ana*100):.1f}%" if total_p_ana > 0 else "0%")
        
        st.markdown(f"**Sensitivity (Recall):** `{(n_c/total_p_ref*100):.1f}%` of published proteins recovered.")
        with st.expander("Show Protein Overlap"): st.dataframe(df_i)

# Sidebar Actions (Download)
st.sidebar.markdown("---")
if st.sidebar.button("📦 DOWNLOAD ENTIRE SESSION (ZIP)"):
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_path = os.path.join(tmp_dir, f"OmicsFlow_Session_{st.session_state['session_id']}")
        shutil.make_archive(zip_path, 'zip', SESSION_DIR)
        with open(f"{zip_path}.zip", "rb") as f:
            st.sidebar.download_button("Download ZIP", f, file_name=f"OmicsFlow_{st.session_state['session_id']}.zip")
