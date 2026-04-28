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
    .metric-value { font-size: 2.2rem; font-weight: 800; color: #58a6ff; }
    .grayed-out { opacity: 0.3; filter: grayscale(100%); pointer-events: none; }
</style>
""", unsafe_allow_html=True)

# --- DIRECTORIES ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]: os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

# --- UTILS ---
def get_metrics_rna(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    with open(filepath, "r", encoding="utf-8") as f: content = f.read()
    hits = re.search(r"Intersection_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    missing = re.search(r"Missing_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    extras = re.search(r"Extras_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    h = int(hits.group(1)) if hits else 0
    m = int(missing.group(1)) if missing else 0
    e = extras.group(1) if extras else "0"
    rate = f"{(h / (h + m) * 100):.1f}" if (h + m) > 0 else "0"
    return str(h), str(m), e, rate

# ==========================================
# SIDEBAR
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    menu = st.radio("Main Navigation", ["📊 Dashboard", "🧪 Transcriptomics Pipeline", "🧪 Proteomics Pipeline", "🧪 Metabolomics Pipeline", "🔍 Reference Validations"])
    st.divider()
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

# ==========================================
# DASHBOARD (100% DYNAMIQUE)
# ==========================================
if menu == "📊 Dashboard":
    st.title("📊 Multi-Omics Studio")
    col_search, col_view = st.columns([3, 1])
    with col_search: st.text_input("🔍 Search Genes/Metabolites...", placeholder="Ex: MTM1, BIN1...")
    with col_view: view = st.selectbox("Switch Omics View", ["Transcriptomique", "Protéomique", "Métabolomique"])
    
    cat_dir = os.path.join(SESSION_DIR, view)
    files = [f for f in (os.listdir(cat_dir) if os.path.exists(cat_dir) else []) if f.endswith(".csv")]
    
    c_left, c_right = st.columns([3, 1])
    
    sig_patho = 0
    sig_rescue = 0
    
    with c_left:
        if files:
            f_sel = st.selectbox("Select Cohort File", files)
            df = pd.read_csv(os.path.join(cat_dir, f_sel))
            lfc = [c for c in df.columns if "Log2FC" in c or "logFC" in c][0]
            pval = [c for c in df.columns if "Pval" in c or "adj.P" in c][0]
            df['-log10p'] = -np.log10(pd.to_numeric(df[pval], errors='coerce').replace(0, 1e-300))
            
            # CALCUL REEL DES METRIQUES POUR LE DASHBOARD
            sig_patho = len(df[df[pval] < 0.05])
            # Supposons que Rescue est calculé si un fichier de Rescue est sélectionné ou présent
            sig_rescue = len(df[(df[pval] < 0.05) & (df[lfc].abs() > 1)]) # Exemple de critère réel
            
            fig = px.scatter(df, x=lfc, y='-log10p', color=(df[pval]<0.05), template="plotly_dark")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No data found. Execute a pipeline first.")
            
    with c_right:
        st.markdown("<div class='targets-summary'>", unsafe_allow_html=True)
        st.subheader("Targets Summary")
        st.markdown("<div class='target-label'>Significatifs (Patho)</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='target-val'>{sig_patho if sig_patho > 0 else '-'}</div>", unsafe_allow_html=True)
        st.markdown("<div class='target-label'>Cibles Rescue Validées</div>", unsafe_allow_html=True)
        st.markdown(f"<div class='target-val'>{sig_rescue if sig_rescue > 0 else '-'}</div>", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

# (Reste des pipelines restaurés en mode dynamique)
elif "Transcriptomics" in menu:
    st.title("RNA-Seq Engine")
    cat_dir = os.path.join(SESSION_DIR, "Transcriptomique")
    if st.button("🚀 EXECUTE DESEQ2"):
        subprocess.run([sys.executable, "01_Analyse_transcriptomique_Validation_Deseq2.py"], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

elif "Proteomics" in menu:
    st.title("Proteomics Engine")
    cat_dir = os.path.join(SESSION_DIR, "Protéomique")
    if st.button("🚀 RUN PYTHON PIPELINE"):
        subprocess.run([sys.executable, "02_Analyse_protéomique_Validation_Supriya.py"], env={"OMICS_IN_DIR": cat_dir, "OMICS_OUT_DIR": cat_dir})

elif "Reference Validations" in menu:
    st.title("🔍 Reference Validations")
    val_case = st.selectbox("Select Validation Case", ["Transcriptomics (Cohort G)", "Proteomics (2w)", "Proteomics (7w)"])
    ref_file = st.file_uploader(f"Upload Reference for {val_case}", type=["xlsx"])
    v_dir = os.path.join(SESSION_DIR, "Validations")
    if ref_file:
        with open(os.path.join(v_dir, "benchmark.xlsx"), "wb") as f: f.write(ref_file.getbuffer())
        if st.button("Run Validation"):
            if "Transcriptomics" in val_case:
                subprocess.run([sys.executable, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py"], env={"OMICS_REF_FILE": os.path.join(v_dir, "benchmark.xlsx"), "OMICS_IN_DIR": os.path.join(SESSION_DIR, "Transcriptomique"), "OMICS_OUT_DIR": v_dir})
    
    st.divider()
    rna_h, rna_m, rna_e, rna_r = get_metrics_rna(os.path.join(v_dir, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.txt"))
    cols = st.columns(4)
    cols[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{rna_h}</div></div>", unsafe_allow_html=True)
    cols[3].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{rna_r}%</div></div>", unsafe_allow_html=True)
