# 🧬 OmicsFlow Studio

Intervace interactive et portable pour l'analyse multi-omique (Transcriptomique, Protéomique, Métabolomique).

## 🚀 Déploiement Cloud (Streamlit Community Cloud)

1. Créez un dépôt sur **GitHub** et uploadez le contenu de ce dossier (`OmicsFlow_App`).
2. Rendez-vous sur [share.streamlit.io](https://share.streamlit.io/).
3. Connectez votre compte GitHub.
4. Sélectionnez votre dépôt et le fichier `app.py`.
5. Dans les paramètres de déploiement, assurez-vous que les secrets ou fichiers nécessaires (ex: références BioC) sont gérés si besoin.

## 🛠️ Installation Locale

```bash
pip install -r requirements.txt
streamlit run app.py
```

## 📂 Structure
- `app.py` : Entrée principale (Dashboard).
- `01_...` : Pipeline Transcriptomique.
- `02_...` : Pipeline Protéomique (Py & R).
- `03_...` : Pipeline Métabolomique.
- `04_...` : Générateur de Volcano Plots.
