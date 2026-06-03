# SynthProteomics: A Probabilistic Fuzzy Rule-Based Framework for Synthetic Clinical and Proteomics Data Simulation in Endometrial Cancer Research

SynthProteomics is introduced as an openly accessible, Gradio-based Python framework for generating the realistic, non-identifiable synthetic datasets across the clinical and proteomic domains. Clinical profiles are simulated using a probabilistic rule-based engine that models the interdependencies among patient factors, including obesity, age, and cancer stage. 

---

## Features

- **Rule-based simulation**: Generate protein expression values guided by biologically informed rules.
- **Clinical factor integration**: Incorporate variables such as age, sex, disease subtype, or treatment to modulate protein expression patterns.
- **Literature-informed**: Utilize known protein interactions, pathways, or co-expression patterns from published studies to increase realism.
- **Customizable datasets**: Control number of patients, proteins, and complexity of relationships.
- **Flexible output**: Export synthetic datasets in common formats (CSV, TSV) for downstream analysis.

---

## Installation

```bash
git clone https://github.com/yourusername/synthetic-proteomics-generator.git
cd synthetic-proteomics-generator
pip install -r requirements.txt

