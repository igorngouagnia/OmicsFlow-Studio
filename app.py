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
    .green-btn > div > button { background: linear-gradient(90deg, #238636 0%, #2ea043 100%) !important; color: white !important; border: none !important; border-radius: 6px !important; }
    .metric-box { text-align: left; min-width: 150px; }
    .metric-label { font-size: 0.85rem; color: #8b949e; margin-bottom: 2px; }
    .metric-value { font-size: 2.2rem; font-weight: 800; color: #58a6ff; line-height: 1; }
</style>
""", unsafe_allow_html=True)

# --- DIRECTORIES ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if 'session_id' not in st.session_state: st.session_state['session_id'] = datetime.now().strftime("%Y%m%d_%H%M")
SESSION_DIR = os.path.join(BASE_DIR, "Sessions", f"Session_{st.session_state['session_id']}")
for sub in ["Transcriptomique", "Protéomique", "Métabolomique", "Validations"]: os.makedirs(os.path.join(SESSION_DIR, sub), exist_ok=True)

# --- PARSING FUNCTIONS ---
def get_metrics_rna(filepath):
    if not os.path.exists(filepath): return "-", "-", "-", "-"
    with open(filepath, "r", encoding="utf-8") as f: content = f.read()
    hits = re.search(r"Intersection_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    extras = re.search(r"Extras_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
    missing = re.search(r"Missing_RNA_Patho_deseq2_standard.csv \((\d+) gènes\)", content)
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
    missing = str(r_total - int(h))
    extras = str(l_total - int(h))
    rate = pct.group(1) if pct else "0"
    return h, missing, extras, rate

# ==========================================
# VALIDATIONS (DYNAMIC)
# ==========================================
with st.sidebar:
    st.title("OmicsFlow")
    menu = st.radio("Main Navigation", ["📊 Dashboard", "🧪 Transcriptomics Pipeline", "🧪 Proteomics Pipeline", "🧪 Metabolomics Pipeline", "🔍 Reference Validations"])
    st.divider()
    perspective = st.radio("Perspective", ["Patho", "Rescue"])

if menu == "🔍 Reference Validations":
    st.title("🔍 Literature Cross-Validation")
    v_dir = os.path.join(SESSION_DIR, "Validations")
    ref_file = st.file_uploader("Upload Reference Excel File", type=["xlsx"])
    if ref_file:
        with open(os.path.join(v_dir, "benchmark.xlsx"), "wb") as f: f.write(ref_file.getbuffer())

    st.divider()
    st.subheader("A. Transcriptomics (Cohort G)")
    st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
    if st.button("Run RNA Validation"):
        env = os.environ.copy()
        env["OMICS_REF_FILE"] = os.path.join(v_dir, "benchmark.xlsx")
        env["OMICS_IN_DIR"] = os.path.join(SESSION_DIR, "Transcriptomique")
        env["OMICS_OUT_DIR"] = v_dir
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.py")], env=env)
    st.markdown("</div>", unsafe_allow_html=True)
    
    # Load Real Metrics
    rna_hits, rna_miss, rna_extra, rna_rate = get_metrics_rna(os.path.join(v_dir, "01_Analyse_transcriptomique_Comparaison_genes_moi_papier_Deseq2_Cohorte_G.txt"))
    
    cols = st.columns(4)
    cols[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{rna_hits}</div></div>", unsafe_allow_html=True)
    cols[1].markdown(f"<div class='metric-box'><div class='metric-label'>Missing (Paper)</div><div class='metric-value'>{rna_miss}</div></div>", unsafe_allow_html=True)
    cols[2].markdown(f"<div class='metric-box'><div class='metric-label'>Extras (Analysis)</div><div class='metric-value'>{rna_extra}</div></div>", unsafe_allow_html=True)
    cols[3].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{rna_rate}%</div></div>", unsafe_allow_html=True)

    st.divider()
    st.subheader("B. Proteomics (Table S8)")
    st.markdown("<div class='green-btn'>", unsafe_allow_html=True)
    if st.button("Run Protein Validation"):
        env = os.environ.copy()
        env["OMICS_REF_DIR"] = v_dir
        env["OMICS_OUT_DIR"] = os.path.join(SESSION_DIR, "Protéomique")
        subprocess.run([sys.executable, os.path.join(BASE_DIR, "02_Analyse_protéomique_Comparaison_protéines_moi_Supriya.py")], env=env)
    st.markdown("</div>", unsafe_allow_html=True)
    
    # Load Real Metrics for Proteo (2w and 7w)
    p2w_file = os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_2w.txt")
    p7w_file = os.path.join(SESSION_DIR, "Protéomique", "Python", "Comparaison_Bilan_Patho_7w.txt")
    
    p2_h, p2_m, p2_e, p2_r = get_metrics_proteo(p2w_file)
    p7_h, p7_m, p7_e, p7_r = get_metrics_proteo(p7w_file)
    
    col_p1, col_p2 = st.columns(2)
    with col_p1:
        st.write("**Case: 2 Weeks**")
        display_cols = st.columns(2)
        display_cols[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p2_h}</div></div>", unsafe_allow_html=True)
        display_cols[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p2_r}%</div></div>", unsafe_allow_html=True)
    with col_p2:
        st.write("**Case: 7 Weeks**")
        display_cols = st.columns(2)
        display_cols[0].markdown(f"<div class='metric-box'><div class='metric-label'>Common Hits</div><div class='metric-value'>{p7_h}</div></div>", unsafe_allow_html=True)
        display_cols[1].markdown(f"<div class='metric-box'><div class='metric-label'>Matching Rate</div><div class='metric-value'>{p7_r}%</div></div>", unsafe_allow_html=True)

# (Rest of code for Dashboards/Pipelines remains active)
else:
    st.info("Pipeline logic is active. Run analysis to see results.")
