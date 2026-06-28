# SynthProteomics

> A Probabilistic Fuzzy Rule-Based Framework for Synthetic Clinical and Proteomics Data Simulation in Endometrial Cancer Research

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/)
[![Gradio](https://img.shields.io/badge/Interface-Gradio-orange.svg)](https://gradio.app)

---

## Overview

**SynthProteomics** is an open-source Python framework for generating realistic, de-identified synthetic clinical and proteomic datasets for endometrial cancer research. It addresses a critical bottleneck in AI-driven oncology: the scarcity of large, well-annotated multi-modal datasets due to privacy laws, data silos, and the high cost of mass spectrometry.

The framework jointly synthesises:
- **Clinical profiles** — using a probabilistic rule-based engine calibrated to Australian cancer registry data (AIHW 2023)
- **Proteomic expression data** — validated against the Human Protein Atlas (HPA) and published endometrial cancer literature

A no-code **Gradio interface** makes the tool accessible to clinical researchers, bioinformaticians, and ML practitioners without programming expertise.

---

## Key Features

| Feature | SynthProteomics | Existing Tools |
|---|---|---|
| Clinical + proteomic joint synthesis | ✅ Yes | ❌ Rarely |
| Validated against Human Protein Atlas | ✅ Yes | ❌ No |
| No-code Gradio interface | ✅ Yes | ❌ No |
| Endometrial-cancer-specific rules | ✅ Yes | ❌ No |
| Open source (MIT) | ✅ Yes | ⚠️ Varies |

---

## How It Works

### Clinical Data Simulation

A probabilistic rule-based engine generates patient attributes with clinically validated interdependencies:

- **Age** — Normal distribution (μ=60, σ=15, clipped 30–85 years)
- **BMI** — Normal distribution (μ=25, σ=4, clipped 21–40)
- **Ethnicity** — 8 categories weighted to Australian demographics (ABS 2021)
- **FIGO tumour stage** — 10 substages (IA through IVB), calibrated to AIHW registry data
- **Histological subtype** — endometrioid, serous, clear cell, mucinous, undifferentiated
- **5-year survival outcome**

### Proteomic Data Simulation

Protein expression profiles are derived from the Human Protein Atlas and modulated using a **Mamdani-inspired fuzzy rule-based system**:

- Oncoproteins (e.g., HER2) are overexpressed; tumour suppressors (e.g., p53, PTEN) are downregulated
- Statistically significant proteins (HPA p≤0.05, TCGA) receive stronger perturbation: log₂FC ≈ ±2 (σ=0.4)
- Non-prognostic proteins receive: log₂FC ≈ ±1 (σ=0.3)
- A fuzzy layer links clinical variables (tumour grade, stage) to proteome perturbation magnitude

### Algorithm Summary

```
Input:  N patients, protein list P, HPA database, regulation dictionary
Output: clinical_data.csv, protein_abundance.csv

Part 1 — Clinical Generation:
  For each patient: sample Age, BMI → derive Menopause, Nulliparity →
  assign TumorType → Grade → FIGO Stage → Treatment

Part 2 — Protein Abundance Generation:
  For each patient × protein:
    Retrieve regulation direction + HPA significance
    Apply fuzzy rule → sample log2FC from N(μ, σ)
    Clip to [-3, 3]
```

Full pseudocode is in [`Algorithm 1`](docs/algorithm.md).

---

## Installation

```bash
# Clone the repository
git clone https://github.com/NehaAr/Synthetic-omics-data.git
cd Synthetic-omics-data

# Install dependencies
pip install -r requirements.txt
```

**Dependencies:** `gradio`, `numpy`, `pandas`, `scikit-learn`, `matplotlib`

### Run in Google Colab (no installation required)

Open the notebook directly: [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/)

---

## Usage

### Launch the Gradio Interface

```bash
python app.py
```

The browser-based interface exposes four parameters:

| Parameter | Description |
|---|---|
| Cohort size | Number of synthetic patients to generate |
| Histological subtype | Filter: endometrioid, serous, clear cell, or all |
| FIGO stage distribution | Filter: Stage 1–4 or all |
| Proteomic panel | Comma-separated gene symbols (e.g. `TP53, PTEN, HER2`) |

Outputs are downloaded as two CSV files:
- `clinical_data.csv`
- `protein_abundance.csv`

Both are immediately compatible with R/Bioconductor and Python scikit-learn workflows.

---

## Validation & Use Cases

### Use Case 1 — FIGO Stage Classifier

Trained a Random Forest (100 estimators) on n=1,000 synthetic patients for early vs. late FIGO stage prediction:

- **82% accuracy** on synthetic held-out test set (80/20 split)
- **74% accuracy** when transferred to real TCGA-UCEC cohort (n=47)

### Use Case 2 — Differential Expression Benchmarking

Benchmarked limma, DESeq2, and Wilcoxon rank-sum on synthetic cohorts (n=200; 100 Stage I, 100 Stage III):

- **87% of seeded differentially expressed proteins** were correctly recovered by limma (FDR < 0.05)

### Use Case 3 — Data Augmentation

Augmented a small real-world dataset (n=47, TCGA-UCEC) with 500 synthetic patients:

| | Real data only | Augmented |
|---|---|---|
| Classification accuracy | 61% | **74%** |
| Cross-validation variance | ±0.12 | **±0.06** |

### Distribution Validation

Synthetic clinical variables were validated against AIHW 2023 cancer registry statistics using Kolmogorov-Smirnov tests — **p > 0.05** for age, BMI, and stage distributions. Proteomic ranges were benchmarked against HPA expression quintiles for matched endometrial tissue.



## Target Users

- **Bioinformaticians** building and evaluating multi-omics analysis pipelines
- **ML researchers** needing labelled cancer datasets for model training and benchmarking
- **Clinical researchers** investigating endometrial cancer biomarkers without access to patient cohorts

---

## Citation

If you use SynthProteomics in your research, please cite:

```bibtex
@article{neha2025synthproteomics,
  title   = {SynthProteomics: A Probabilistic Fuzzy Rule-Based Framework for
             Synthetic Clinical and Proteomics Data Simulation in Endometrial Cancer Research},
  author  = {Neha},
  year    = {2025},
  note    = {University of Newcastle, Callaghan, NSW, Australia},
  url     = {https://github.com/NehaAr/Synthetic-omics-data}
}
```

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## Contact

**Neha**  
PhD Candidate, Bioinformatics  
Department of Biomedical Sciences and Pharmacy  
University of Newcastle, Callaghan, NSW 2308, Australia  
📧 neha10@uon.edu.au

---

## Acknowledgements

- [Human Protein Atlas](https://www.proteinatlas.org/) for tissue-level expression reference data
- [TCGA-UCEC](https://portal.gdc.cancer.gov/) for real-world validation cohort
- [Australian Institute of Health and Welfare](https://www.aihw.gov.au/) for cancer registry statistics
