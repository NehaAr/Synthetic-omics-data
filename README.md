

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1MYcSHwfEBKkm-AtIisCxr13qDZBFFxiY?usp=sharing)
[![Gradio](https://img.shields.io/badge/Interface-Gradio-orange.svg)](https://gradio.app)

---

# SynthProteomics

**A probabilistic, fuzzy rule-based framework for synthetic clinical and proteomic data simulation in endometrial cancer research.**

SynthProteomics is a no-code, browser-based tool (built with [Gradio](https://gradio.app)) that generates realistic, fully synthetic clinical profiles and protein abundance data for endometrial cancer. It is designed for bioinformaticians, machine learning researchers, and clinical researchers who need labeled, biologically plausible data without access to real patient cohorts or mass spectrometry facilities.

No programming experience is required to use the tool — everything runs through an interactive web interface.

---

## Table of contents

- [What this tool does](#what-this-tool-does)
- [Quick start (Google Colab)](#quick-start-google-colab)
- [Quick start (local install)](#quick-start-local-install)
- [Optional: adding real HPA expression data](#optional-adding-real-hpa-expression-data)
- [Using the interface](#using-the-interface)
  - [Tab 1 — Clinical Data](#tab-1--clinical-data)
  - [Tab 2 — Protein Abundance Data](#tab-2--protein-abundance-data)
  - [Tab 3 — Distribution Plots](#tab-3--distribution-plots)
  - [Tab 4 — Use Case: Stage Classifier](#tab-4--use-case-stage-classifier)
- [How the simulation works](#how-the-simulation-works)
- [Output files](#output-files)
- [Known limitations](#known-limitations)
- [Citation](#citation)
- [License](#license)

---

## What this tool does

SynthProteomics generates two linked synthetic datasets:

1. **Clinical data** — synthetic patient records (age, BMI, menopausal status, tumor grade/type, FIGO-style stage, histological subtype, nulliparity, treatment, survival outcome), generated using a probabilistic rule engine that encodes known clinical dependencies (e.g. higher age + nulliparity + higher BMI increases the probability of a Type 2 tumor; tumor grade shapes the stage distribution; stage and grade jointly shape survival outcome).

2. **Protein abundance data** — synthetic log2 fold-change (log2FC) values for a user-specified list of proteins, generated for the patients created in step 1. Each value is sampled based on a literature-curated regulation direction (UP/DOWN/LOW), the protein's prognostic significance, and (where relevant) the patient's clinical stage/grade/menopausal status.

The two datasets share patient identifiers implicitly through generation order, so they can be joined for downstream multi-omics analysis.

---

## Quick start (Google Colab)

This is the easiest way to run SynthProteomics — no local setup required.

1. Open a new [Google Colab](https://colab.research.google.com/) notebook.
2. Copy the full script from `synthproteomics.py` (or this repo's main script) into a cell.
3. Run the cell. The script will automatically install `gradio` and `scikit-learn`.
4. A public Gradio link will be printed at the bottom of the output — click it to open the app in a new tab.

```python
# Cell 1 — just paste and run the whole script
# (pip installs happen automatically inside the script)
```

> **Tip:** the script looks for an optional file at `/content/NOT_new_rna.tsv` (see [Optional: adding real HPA expression data](#optional-adding-real-hpa-expression-data) below). If you don't upload this file, the tool still runs — it falls back to a neutral baseline and prints a warning.

---

## Quick start (local install)

```bash
git clone https://github.com/NehaAr/Synthetic-omics-data.git
cd Synthetic-omics-data
pip install -r requirements.txt
python synthproteomics.py
```

This will start a local Gradio server (by default at `http://127.0.0.1:7860`). Open that address in your browser.

**Requirements:**
- Python 3.9+
- `gradio`
- `pandas`
- `numpy`
- `matplotlib`
- `scikit-learn` (only required for the Stage Classifier tab)

> **Note:** the script currently reads the optional HPA file from a hardcoded path (`/content/NOT_new_rna.tsv`), which is a Google Colab convention. If running locally, either create a `/content/` directory with the file in it, or edit the `pd.read_csv(...)` path near the top of the script to point at your own file location.

> **Note:** CSV downloads are written to `/tmp/`, which exists by default on Linux/macOS. Windows users should edit the `download_clinical_csv()` / `download_protein_csv()` functions to use a Windows-compatible temp path (e.g. `tempfile.gettempdir()`).

---

## Optional: adding real HPA expression data

The tool can incorporate real baseline tissue expression data from the [Human Protein Atlas](https://www.proteinatlas.org/) to refine its simulation. This is **optional** — the tool works without it.

To enable it:

1. Go to the Human Protein Atlas and download the normal tissue RNA expression + cancer prognostics export (TSV format) for endometrial tissue.
2. Rename the file `NOT_new_rna.tsv`.
3. In Colab: upload it to `/content/NOT_new_rna.tsv` (use the Colab file browser, or `files.upload()`).
4. In a local install: either place it at `/content/NOT_new_rna.tsv` or edit the file path in the script.

The file must contain at least these columns:
- `Gene`
- `Tissue RNA - endometrium 1 [nTPM]`
- `Cancer prognostics - Uterine Corpus Endometrial Carcinoma (TCGA)`

If the file is missing, the tool prints a warning and falls back to a neutral ("moderate" expression, "no" prognostic significance) baseline for every protein not otherwise found in the literature-curated dictionary described below.

---

## Using the interface

The app opens with four tabs.

### Tab 1 — Clinical Data

Generates the synthetic patient cohort.

| Control | Description |
|---|---|
| **Number of Patients** | Cohort size to generate (minimum 10) |
| **Cancer Subtype Filter** | Restrict output to one histological subtype (`Endometrioid`, `Serous`, `Clear Cell`, `Mucinous`, `Undifferentiated`), or `All` |
| **Stage Filter** | Restrict output to one FIGO stage group (`Stage1`–`Stage4`), or `All` |
| **Select Columns to Display** | Choose which fields appear in the results table |

Click **Generate Clinical Data** to run. Click **Download as CSV** to save the full cohort (`clinical_data.csv`) regardless of which columns are displayed.

**Fields generated:** `Patient_ID`, `Ages`, `Ethnicity`, `Menopause`, `Grade`, `Tumor_type`, `Histological_Subtype`, `Stage`, `Nulliparity`, `BMI`, `Myometrial_mm`, `Treatment`, `Survival_Outcome`.

> You must generate clinical data in this tab before using Tab 2 or Tab 4, since both depend on the patient cohort created here.

### Tab 2 — Protein Abundance Data

Generates synthetic log2FC protein abundance values for the current patient cohort.

| Control | Description |
|---|---|
| **Number of Patients** | How many patients (from the cohort generated in Tab 1) to generate protein data for |
| **Protein List** | Comma-separated gene symbols, e.g. `ANXA2, PKM2, ERBB2, EGFR, MMP9` |

Click **Generate Abundance Data** to view results as JSON, or **Download as CSV** to export a patient × protein matrix.

Values are log2 fold-changes clipped to the range **[-3, 3]**. A protein not found in either the literature-curated dictionary or the optional HPA file will still receive a value, sampled from a near-zero neutral distribution — it is not silently dropped.

### Tab 3 — Distribution Plots

Click **Generate Plots** to visualize the distributions of all generated clinical variables (age, BMI, grade, stage, tumor type, menopause, nulliparity, treatment, survival outcome) as histograms/bar charts. Requires clinical data to have been generated in Tab 1 first.

### Tab 4 — Use Case: Stage Classifier

Demonstrates a downstream machine learning use case: trains a Random Forest classifier (`scikit-learn`, 100 estimators, 80/20 train-test split, fixed random seed) on the synthetic clinical cohort to predict **early stage** (Stage 1–2) vs. **late stage** (Stage 3–4) disease from `Ages`, `BMI`, `Menopause`, `Nulliparity`, `Tumor_type`, and `Grade`.

Click **Run Random Forest Classifier** to view accuracy and a full classification report. Requires clinical data to have been generated in Tab 1 first.

---

## How the simulation works

**Clinical data generation** follows a cascade of probabilistic rules:
- Age ~ Normal(μ=60, σ=15), clipped to 30–85
- BMI ~ Normal(μ=25, σ=4), clipped to 21–40
- Menopausal status is deterministic on age (post-menopausal above 51)
- Nulliparity probability increases with age
- A composite "high-risk" profile (age ≥ 60, nulliparous, BMI ≥ 25) raises the probability of a Type 2 tumor
- Tumor type shapes the grade distribution; grade shapes the stage distribution; stage and grade jointly shape the simulated survival outcome

**Protein abundance generation** assigns each (patient, protein) pair a log2FC value via a small decision layer that considers:
1. The protein's literature-curated regulation direction (`UP`, `DOWN`, `UP/DOWN`, or `LOW`) from an internal, hand-curated dictionary
2. Its prognostic significance (from the optional HPA/TCGA file, if provided)
3. Its baseline tissue expression level (from the same optional file)
4. Any clinical-factor specificity noted in the literature (e.g. a protein only reported as dysregulated in Grade 3 or postmenopausal patients) — the simulator checks whether the current patient actually matches that factor before applying the stronger perturbation

Each of these combinations maps to one of six Gaussian samplers (strong/weak, up/down, low, or neutral), all clipped to stay within the biologically motivated [-3, 3] log2FC range.

---

## Output files

| File | Produced by | Contents |
|---|---|---|
| `clinical_data.csv` | Automatically saved whenever Tab 1 is run | Full synthetic clinical cohort |
| `SynthProteomics_clinical.csv` | "Download as CSV" button, Tab 1 | Same as above, via the download button |
| `SynthProteomics_protein.csv` | "Download as CSV" button, Tab 2 | Patient × protein log2FC matrix |

---

## Known limitations

- The literature-curated protein regulation dictionary currently covers a fixed, hand-compiled list of genes reported in endometrial cancer literature. Proteins outside this list (and outside the optional HPA file) are assigned neutral/near-zero values rather than being excluded.
- Proteins are sampled independently of one another; the tool does not currently model pathway-level co-expression or protein–protein correlation structure.
- All clinical and proteomic values are fully synthetic and are intended for method development, teaching, and pipeline testing — not for drawing new biological conclusions about endometrial cancer.
- File paths for the optional HPA input and CSV downloads assume a Google Colab environment (`/content/`, `/tmp/`) by default; see [Quick start (local install)](#quick-start-local-install) for adjustments needed to run outside Colab.

---

## Citation

If you use SynthProteomics in your research, please cite:
```
@software{Neha_Arora_SynthProteomics,
author = {{Neha Arora}},
title = {{SynthProteomics}},
version = {1.0.0}
}
```
---

## License

MIT License. See [`LICENSE`](LICENSE) for details.

## Author

**Neha Arora** — PhD candidate, Bioinformatics, University of Newcastle, Callaghan, NSW, Australia.
Repository: [https://github.com/NehaAr/Synthetic-omics-data](https://github.com/NehaAr/Synthetic-omics-data)
