# -*- coding: utf-8 -*-
"""
SynthProteomics: A Probabilistic Rule-Based Framework for Synthetic Clinical
and Proteomics Data Simulation in Endometrial Cancer Research

Author: Neha Arora
Affiliation: University of Newcastle, Callaghan, Australia
Copyright (c) 2026 Neha Arora — MIT License
GitHub: https://github.com/NehaAr/Synthetic-omics-data
"""
import subprocess, sys
subprocess.run([sys.executable, "-m", "pip", "install", "gradio", "--quiet"], check=False)
subprocess.run([sys.executable, "-m", "pip", "install", "scikit-learn", "--quiet"], check=False)

import ast
import re
import io
import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random

# ══════════════════════════════════════════════════════════════════════
# LITERATURE-COMPILED PROTEIN REGULATION DATABASE
#
# Columns: Gene, Regulation (UP/DOWN/UP/DOWN/LOW), Factor (stage/grade/
#          menopausal status specificity; NA = general)
# ══════════════════════════════════════════════════════════════════════
protein_abundance_dictionary = {
    "Gene": [
        'IFIT3','PARP9','SLC34A2','CYB5R1','PTPN1','DPT','SLP','ANXA2','PRDX1','CTNNB',
        'HMGB3','CLC1','EIF4A1','PRDX6','ENO1','ANXA4','EMD','KU70','GRP78','GSTP1',
        'ACTG','DIA3','ENOA','ALBU','ACTB','ACTG','KRT8','ANXA1','ENOA','TRFE','HSPB1',
        'EF-TU','IDH1','SOD1','CALR','RPSA','UAP56','PSME1','PDIA3','CAH1','IDHC','PPIA',
        'PPIB','ZNF844','ALDOA','ENO1','KRT10','ABRACL','PGAM2','FGB','ANXA3','CTNB1',
        'XPO2','CAPG','MMP9','EWSR1','PKM','NAMPT','ENOA','CATD','LDHA','SPIT1','OSTP',
        'MPO','CASP3','CADH1','TUBA1A','TIGAR','SEC11A','CENPV','TMSB4X','COL1A2',
        'S100A16','NEBL','OGN','COL1A2','S100A16','SLC9A3R1','DNAJB11','RBBP4','MYH11',
        'OGN','GNLY','MZB1','MX1','NANS','TMED9','TPPP3','HNRNPF','NOLC1','SLC4A1',
        'COL5A1','FGA','HBA1','COL1A2','SNRPC','UBE2V2','COL1A1','BCAM','PTMA','DEFA1',
        'S100A8','LTF','CAMP','AZU1','COL1A2','SEC63','LDHB','ABHD14B','LTF','SARS1',
        'ATP1B1','IARS1','PNP','SFN','ACTA2','TPR','MAP4','HBG2','PSMD11','SLC2A11',
        'SLC2A1','SRRM2','U2AF1','TMSB4X','DVL-2','HSP8','EIF4G2','F13A1','GFM1',
        'NPEPL1','SARS2','SNTB1','UBR4','USP47','WDR5','ASS1','PFAS','ckb','hk2',
        'MAPKAPK3','SERPINA1','TAGLN2','TPI1','ENO1','ANXA2','ANXA1','SRPK1','PTPN6',
        'ATP5A','TUBB','ERBB2','EGFR','ACTN4','UBE2N','PKM2','HSPA5','LMNAA/C','HRNR',
        'MDH2','STMN1','CKB','DJ-1','PRKCSH','NCL','GON7','APOA1','S100A','PKM2',
        'HSP10','EC1','EC2','PKM1','WFDC2','CLU','MUC5B','COX2','PRDX3','PRDX5','PRDX6',
        'RDX4','HNRNPA1','CTSB','CALU','CACYBP','LDHA','P38MAPK','NFKB','ERK1/2','PRDX6',
        'EIF4A1','CLIC1','CLIC4','TAGLN2','TPI1','TPI-1','HSPA8','ATF2','JUN','TAF1',
        'HNF4A','ATF7IP','ERBB2','EGFR','JPT1','CALR','RPSA','ACTB','IDH1','MLH1',
        'EPM2AIP1','SRPK1','CAPS','PRTN3','HMGA2','S100A8','LTF','CTSD','STMN1','TAGLN',
        'DES','CALD1','CNN1','CDH13','PARP9','IFIT3','DPT','SLP1','OXSR1','ASS1'
    ],
    "Regulation": [
        'UP','UP','UP','UP','UP','DOWN','UP','UP','UP','DOWN','UP','DOWN','UP','UP','UP',
        'DOWN','UP/DOWN','UP/DOWN','UP','UP','UP','DOWN','UP/DOWN','UP','DOWN','DOWN','UP',
        'UP','UP','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN',
        'UP/DOWN','UP','UP','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN',
        'UP/DOWN','UP/DOWN','UP/DOWN','UP','UP','UP','UP','DOWN','DOWN','DOWN','DOWN',
        'DOWN','UP','UP','UP','UP','UP','DOWN','DOWN','DOWN','UP','UP','UP','UP','UP',
        'UP','UP','DOWN','DOWN','DOWN','DOWN','DOWN','UP','UP','UP','UP','UP','DOWN',
        'DOWN','DOWN','DOWN','UP','UP','UP','UP','UP','DOWN','DOWN','DOWN','DOWN','DOWN',
        'UP','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN',
        'UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN',
        'UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP/DOWN','UP','UP/Down','UP/Down',
        'UP/Down','UP/Down','UP/DOWN','UP/Down','UP/Down','UP/Down','UP/Down','UP','UP',
        'UP','UP','UP','LOW','LOW','LOW','LOW','UP','UP','UP','UP','UP','LOW','LOW',
        'UP/DOWN','UP/DOWN','UP','UP','DOWN','UP','UP','UP','UP','DOWN','UP','UP','UP',
        'UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP',
        'UP','UP','DOWN','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP','UP',
        'UP','UP','UP','UP','DOWN','DOWN','UP','UP','UP','UP','UP','UP','UP','UP','DOWN',
        'DOWN','DOWN','DOWN','DOWN','UP','UP','DOWN','DOWN','DOWN','DOWN'
    ],
    "Factor": [
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','STAGE1A','STAGE2,stage3','NA','STAGE2','NA','NA','STAGE1B,stage3','NA',
        'STAGE3','NA','STAGE3','STAGE3','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','POSTMENOPAUSAL','POSTMENOPAUSAL',
        'POSTMENOPAUSAL','POSTMENOPAUSAL','POSTMENOPAUSAL','POSTMENOPAUSAL','POSTMENOPAUSAL',
        'POSTMENOPAUSAL','POSTMENOPAUSAL','NA','NA','NA','NA','NA','NA','NA','NA','GRADE2',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','Grade3','Grade3','Grade3',
        'Grade3','Grade3','NA','NA','NA','NA','NA','Grade3','Grade3','Grade3','Grade3',
        'Grade3','Grade3','Grade3','Grade3','Grade3','Grade3','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','TYPE 2','TYPE 1','NA','NA','NA','NA',
        'GRADE1,grade3,stage1,stage3','GRADE1,GRADE3','NA','TYPE1','TYPE1','TYPE1',
        'STAGE1A,STAGE1B','NA','NA','NA','NA','type 1','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','stage1B','NA','NA','NA','stage1A','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA'
    ]
}

# ══════════════════════════════════════════════════════════════════════
# GLOBAL STATE
# CHANGE: clinical_data initialised as DataFrame not list (BUG FIX 6)
# ══════════════════════════════════════════════════════════════════════
clinical_data   = pd.DataFrame()
protein_data    = []

# ══════════════════════════════════════════════════════════════════════
# HPA NORMAL TISSUE EXPRESSION (Human Protein Atlas)
# Reference: Uhlen M, et al. Science 347:1260419 (2015)
# TPM thresholds per paper Section 2.2:
#   >100 nTPM = highly expressed (UP)
#   10-100 nTPM = moderately expressed (MODERATE)
#   <10 nTPM = lowly expressed (DOWN)
# ══════════════════════════════════════════════════════════════════════
try:
    normal_tissue_expression = pd.read_csv('/content/NOT_new_rna.tsv', sep='\t')

    normal_tissue_prognostic = (
        normal_tissue_expression[
            'Cancer prognostics - Uterine Corpus Endometrial Carcinoma (TCGA)'
        ]
        .astype(str)
        .str.strip()
        .str.replace(r'[^\d.]', '', regex=True)   # remove non-numeric chars cleanly
    )
    normal_tissue_prognostic = pd.to_numeric(normal_tissue_prognostic, errors='coerce').astype(float)

    # TPM-based expression classification (matches paper Section 2.2)
    normal_tisue_regulation = normal_tissue_expression[
        'Tissue RNA - endometrium 1 [nTPM]'
    ].apply(lambda x: "up" if x > 100 else ("moderate" if 10 <= x <= 100 else "down"))

    # Prognostic significance: p≤0.05 = statistically prognostic
    gene_prognosis_indicator = normal_tissue_prognostic.apply(
        lambda x: "yes" if pd.notna(x) and x <= 0.05 else "no"
    )
    gene_list_lower = normal_tissue_expression['Gene'].str.lower().tolist()

except FileNotFoundError:
    print("WARNING: HPA TSV not found. Upload NOT_new_rna.tsv to /content/")
    normal_tissue_expression = pd.DataFrame(columns=['Gene'])
    normal_tisue_regulation  = pd.Series(dtype=str)
    gene_prognosis_indicator = pd.Series(dtype=str)
    gene_list_lower          = []


# ══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════

def parse_factors(factor_string):
    """
    Parse Factor field which uses mixed ',' and '/' delimiters.
    BUG FIX 1: original code used non-existent .splt() method.
    Now uses re.split() to correctly tokenise entries like
    'STAGE1A,STAGE1B' or 'GRADE1,grade3,stage1,stage3'.
    """
    if not factor_string or factor_string.upper() == "NA":
        return []
    parts = re.split(r'[/,]', factor_string)
    return [p.strip().lower() for p in parts if p.strip()]


def factor_matches_patient(factor_string, stage_arr, grade_arr, menopause_arr):
    """
    Check whether a protein's factor condition matches the current patient cohort.
    This implements the 'fuzzy rule-based layer' described in paper Section 2.2 —
    clinical parameters (stage/grade/menopausal status) gating proteomic perturbation.

    BUG FIX 4: original code used + operator on numpy arrays which concatenates
    as matrix not list. Now uses np.concatenate() correctly.
    """
    factors = parse_factors(factor_string)
    if not factors:
        return False
    all_values = [
        str(v).lower()
        for v in np.concatenate([
            np.asarray(stage_arr).flatten(),
            np.asarray(grade_arr).flatten(),
            np.asarray(menopause_arr).flatten()
        ])
    ]
    return any(f in all_values for f in factors)


# ── Gaussian abundance samplers ───────────────────────────────────────
# log2FC ranges per paper: strong signal ±2 to ±3, weak ±1 to ±2
# These correspond to the fuzzy output membership functions:
#   prognostic=YES  → strong perturbation (loc=±2, scale=0.4)
#   prognostic=NO   → weak perturbation   (loc=±1, scale=0.3)
#   neutral/unknown → baseline noise      (loc=0,  scale=0.5)

def sample_up_strong():
    """Upregulated, prognostically significant — log2FC ~2"""
    return float(np.clip(np.random.normal(loc=2.0, scale=0.4), 0.0, 3.0))

def sample_up_weak():
    """Upregulated, not prognostically significant — log2FC ~1"""
    return float(np.clip(np.random.normal(loc=1.0, scale=0.3), 0.0, 3.0))

def sample_down_strong():
    """Downregulated, prognostically significant — log2FC ~ -2"""
    return float(np.clip(np.random.normal(loc=-2.0, scale=0.4), -3.0, 0.0))

def sample_down_weak():
    """Downregulated, not prognostically significant — log2FC ~ -1"""
    return float(np.clip(np.random.normal(loc=-1.0, scale=0.3), -3.0, 0.0))

def sample_low():
    """
    LOW regulation — mild downregulation distinct from full DOWN.
    NEW: original code had no handler for 'LOW' regulation values
    that appear in the dictionary (e.g. SLC2A1, SLC2A11).
    """
    return float(np.clip(np.random.normal(loc=-0.5, scale=0.3), -1.5, 0.0))

def sample_neutral():
    """No clear direction or factor not matched — near-zero noise"""
    return float(np.clip(np.random.normal(loc=0.0, scale=0.5), -1.5, 1.5))


def assign_abundance(regulation, normal_reg, prognostic):
    """
    Core Mamdani-inspired assignment: maps fuzzy linguistic inputs
    (regulation direction × prognostic significance × baseline expression)
    to a Gaussian-sampled log2FC output value.

    This is the central 'fuzzy rule-based layer' referenced in paper Section 2.2.
    Linguistic variables:
      - regulation:  UP | DOWN | UP/DOWN | LOW
      - normal_reg:  up | down | moderate  (from HPA nTPM)
      - prognostic:  yes | no             (from TCGA p-value)
    """
    reg  = regulation.lower().strip()
    prog = str(prognostic).lower().strip()

    if reg == 'up':
        return sample_up_strong() if prog == 'yes' else sample_up_weak()

    elif reg == 'down':
        return sample_down_strong() if prog == 'yes' else sample_down_weak()

    elif 'up' in reg and 'down' in reg:  # catches UP/DOWN and UP/Down
        if normal_reg == 'up':
            return sample_up_strong()  if prog == 'yes' else sample_up_weak()
        elif normal_reg == 'down':
            return sample_down_strong() if prog == 'yes' else sample_down_weak()
        else:
            return sample_neutral()

    elif reg == 'low':
        # BUG FIX 7: LOW was previously unhandled → fell through to neutral
        return sample_low()

    else:
        return sample_neutral()


# ══════════════════════════════════════════════════════════════════════
# SECTION 2.1 — CLINICAL DATA SIMULATION
# Paper: "probabilistic rule-based engine encoding clinically established
# interdependencies among patient attributes"
# Primary variables per paper: age, BMI, ethnicity, tumor stage (FIGO),
# histological subtype, treatment modality, survival outcome
#
# CHANGES FROM ORIGINAL:
# + Added ethnicity (listed in paper Section 2.1 as primary variable)
# + Added survival_outcome (listed in paper Section 2.1 as primary variable)
# + Added subtype (endometrioid/serous/clear cell per paper Section 2.1)
# + Added subtype_filter and stage_filter params for Gradio Section 2.3
# + BUG FIX 2: Grade comparison corrected from 'Type1' → 'Grade1'
# + BUG FIX 6: clinical_data now pd.DataFrame not list
# ══════════════════════════════════════════════════════════════════════
def generate_person_data(num_records, selected_columns,
                         subtype_filter="All", stage_filter="All"):
    """
    Generates synthetic clinical profiles for endometrial cancer patients.

    Rules encoded (per FIGO staging and published epidemiology):
    - Age ≤ 51        → pre-menopausal (Menopause=1)
    - Age > 51        → post-menopausal (Menopause=2)
    - Age ≤ 35        → 50/50 nulliparity
    - Age > 35        → 60% nulliparity (higher parity with age)
    - Age ≥ 60 + nulliparous + BMI ≥ 25 → high-risk: 60% Type2 tumour
    - Type1 tumour    → Grade1 more probable (50/30/20)
    - Type2 tumour    → Grade3 more probable (20/30/50)
    - Grade1/2        → early stage weights (Stage1-2 more probable)
    - Grade3          → late stage weights  (Stage3-4 more probable)
    """
    global clinical_data
    num_records = int(num_records)

    np.random.seed(42)
    random.seed(42)

    # ── Age: Normal distribution, clipped to 30–85 ───────────────────
    Ages = np.clip(
        np.random.normal(loc=60, scale=15, size=num_records).astype(int),
        30, 85
    )

    # ── BMI: Normal distribution, clipped to 21–40 ───────────────────
    bmi_values = np.round(
        np.clip(np.random.normal(loc=25, scale=4, size=num_records), 21, 40),
        2
    )

    # ── Ethnicity: Added per paper Section 2.1 ───────────────────────
    # Distribution reflects Australian/Western population estimates
    #ethnicities = random.choices(
       # ['Caucasian','Asian','African','Hispanic','Indigenous Australian',
       #  'Middle Eastern','South Asian','Other'],
       # weights=[0.60, 0.12, 0.08, 0.06, 0.04, 0.04, 0.04, 0.02],
        #k=num_records
    #)
    ethnicities = random.choices(
    ['European','East Asian','South Asian','Middle Eastern',
     'African','Indigenous Australian','Other'],
    weights=[0.55, 0.15, 0.10, 0.05, 0.03, 0.03, 0.09],
    k=num_records
     )
    patient_id, Treatment, Menopause = [], [], []
    Grade, Stage, Myometrial_invasion = [], [], []
    Nulliparity, Tumor_type, Subtype  = [], [], []
    Survival_outcome = []  # Added per paper Section 2.1

    for i in range(num_records):
        patient_id.append(f"EC{i:04d}")
        Treatment.append(random.choice(["Surgery", "Chemotherapy",
                                        "Radiotherapy", "Combined"]))
        # Myometrial invasion: 2–40mm range (normal ≤5mm; >50% = Stage IB)
        Myometrial_invasion.append(round(random.uniform(2, 40), 2))

    # ── Menopause: clinical threshold age 51 ─────────────────────────
    for age in Ages:
        Menopause.append(1 if age <= 51 else 2)  # 1=pre, 2=post

    # ── Nulliparity: probability increases with age ───────────────────
    for age in Ages:
        if age <= 35:
            Nulliparity.extend(random.choices(["yes", "no"], [0.5, 0.5], k=1))
        else:
            Nulliparity.extend(random.choices(["yes", "no"], [0.6, 0.4], k=1))

    # ── Stage choices (FIGO system) ───────────────────────────────────
    stage_choices = [
        'Stage1','Stage1A','Stage1B',       # confined to uterus
        'Stage2','Stage2A','Stage2B',       # cervical spread
        'Stage3C1','Stage3C2',              # lymph node involvement
        'Stage4A','Stage4B'                 # distant metastasis
    ]
    early_weights = [0.1333]*6 + [0.05]*4   # Grade1/2: stages 1-2 more likely
    late_weights  = [0.05]*6   + [0.1333]*4 # Grade3: stages 3-4 more likely

    # ── Histological subtype per paper Section 1 ─────────────────────
    # Endometrioid: 70-80%; serous/clear cell: rare
    subtype_map = {
        'Type1': random.choices(
            ['Endometrioid','Mucinous'], [0.90, 0.10], k=1
        )[0],
        'Type2': random.choices(
            ['Serous','Clear Cell','Undifferentiated'], [0.55, 0.35, 0.10], k=1
        )[0]
    }

    for idx, age in enumerate(Ages):
        # ── Core probabilistic rules ──────────────────────────────────
        # Rule 1: High-risk profile → elevated Type2 probability
        high_risk = (age >= 60 and
                     Nulliparity[idx] == 'yes' and
                     bmi_values[idx] >= 25)

        if high_risk:
            Tumor_type.extend(random.choices(['Type1','Type2'], [0.4, 0.6], k=1))
        else:
            Tumor_type.extend(random.choices(['Type1','Type2'], [0.8, 0.2], k=1))

        # Rule 2: Tumor type → Grade distribution
        if Tumor_type[idx] == 'Type1':
            Grade.extend(random.choices(
                ['Grade1','Grade2','Grade3'], [0.5, 0.3, 0.2], k=1))
            Subtype.append(random.choices(
                ['Endometrioid','Mucinous'], [0.90, 0.10], k=1)[0])
        else:
            Grade.extend(random.choices(
                ['Grade1','Grade2','Grade3'], [0.2, 0.3, 0.5], k=1))
            Subtype.append(random.choices(
                ['Serous','Clear Cell','Undifferentiated'], [0.55, 0.35, 0.10], k=1)[0])

        # Rule 3: Grade → Stage weights
        # BUG FIX 2: was Grade[idx]=='Type1' — wrong type comparison
        if Grade[idx] in ('Grade1', 'Grade2'):
            Stage.extend(random.choices(stage_choices, weights=early_weights, k=1))
        else:
            Stage.extend(random.choices(stage_choices, weights=late_weights, k=1))

        # Rule 4: Stage and grade influence survival outcome
        # (added per paper Section 2.1 — survival outcome is a primary variable)
        if Stage[idx] in ('Stage1','Stage1A','Stage1B') and Grade[idx] == 'Grade1':
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.90, 0.10], k=1)[0])
        elif Stage[idx] in ('Stage4A','Stage4B'):
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.30, 0.70], k=1)[0])
        else:
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.65, 0.35], k=1)[0])

    clinical_data = pd.DataFrame({
        'Patient_ID':       patient_id,
        'Ages':             Ages,
        'Ethnicity':        ethnicities,      # ADDED: per paper Section 2.1
        'Menopause':        Menopause,
        'Grade':            Grade,
        'Tumor_type':       Tumor_type,
        'Histological_Subtype': Subtype,      # ADDED: per paper Section 2.1
        'Stage':            Stage,
        'Nulliparity':      Nulliparity,
        'BMI':              bmi_values,
        'Myometrial_mm':    Myometrial_invasion,
        'Treatment':        Treatment,
        'Survival_Outcome': Survival_outcome  # ADDED: per paper Section 2.1
    })

    # ── Subtype filter (per paper Section 2.3 — user specifies subtype) ──
    if subtype_filter != "All":
        clinical_data = clinical_data[
            clinical_data['Histological_Subtype'] == subtype_filter
        ].reset_index(drop=True)

    # ── Stage filter (per paper Section 2.3 — user specifies stage) ──
    if stage_filter != "All":
        clinical_data = clinical_data[
            clinical_data['Stage'].str.startswith(stage_filter)
        ].reset_index(drop=True)

    clinical_data.to_csv('clinical_data.csv', index=False)

    valid_cols = [c for c in selected_columns if c in clinical_data.columns]
    return clinical_data[valid_cols] if valid_cols else clinical_data


# ══════════════════════════════════════════════════════════════════════
# SECTION 2.2 — PROTEOMIC DATA SIMULATION
# Paper: "fuzzy rule-based layer associates clinical parameters (tumor
# grade and stage) with magnitude of proteomic perturbations, capturing
# continuous and graded nature of biological variation"
#
# SIX CHECKPOINTERS — now explicitly named and documented:
#   CP1: In literature DB + In HPA + Factor=NA   (general regulation)
#   CP2: NOT in literature  + In HPA             (HPA-only)
#   CP3: In literature      + NOT in HPA + NA    (literature-only, safe)
#   CP4: In literature      + NOT in HPA + Factor (stage/grade-gated)
#   CP5: Neither source available                 (neutral fallback)
#   CP6: In literature + In HPA + Factor          (full gated)
#
# BUG FIX 8: generate_abundance dict now initialised per patient inside
# the outer loop — original was a global dict causing cross-patient
# contamination of abundance values.
# ══════════════════════════════════════════════════════════════════════
def generate_protein_abundance_data(num_records, protein_list):
    """
    Generates log2 fold-change protein abundance values per patient.

    Output range: [-3, 3] — consistent with published proteomic studies
    as stated in paper comments and Section 2.2.

    For each protein × patient, the fuzzy rule-based layer determines
    abundance by evaluating:
      1. Regulation direction (UP/DOWN/UP/DOWN/LOW) from literature
      2. Prognostic significance (p≤0.05) from HPA/TCGA
      3. Baseline expression level (nTPM) from HPA
      4. Stage/grade/menopausal factor specificity from literature
    """
    all_outputs = []
    num_records  = int(num_records)

    stage_arr    = clinical_data['Stage'].unique()    if not clinical_data.empty else []
    grade_arr    = clinical_data['Grade'].unique()    if not clinical_data.empty else []
    menopause_arr= clinical_data['Menopause'].unique() if not clinical_data.empty else []

    for j in range(num_records):
        # BUG FIX 8: initialise fresh dict per patient — original was global
        generate_abundance = {}

        for i in protein_list:
            protein_lower = i.strip().lower()

            # ── Index lookup ──────────────────────────────────────────
            try:
                index_protein = [g.lower() for g in
                                 protein_abundance_dictionary['Gene']].index(protein_lower)
            except ValueError:
                index_protein = None

            try:
                index_normal = gene_list_lower.index(protein_lower)
            except ValueError:
                index_normal = None

            regulation   = protein_abundance_dictionary['Regulation'][index_protein] \
                           if index_protein is not None else None
            factor_str   = protein_abundance_dictionary['Factor'][index_protein] \
                           if index_protein is not None else "NA"
            normal_reg   = normal_tisue_regulation.iloc[index_normal] \
                           if index_normal is not None else "moderate"
            # BUG FIX 3: original code accessed gene_prognosis_indicator[None]
            # in CP3/CP4 when index_normal was None → RuntimeError
            prognostic   = gene_prognosis_indicator.iloc[index_normal] \
                           if index_normal is not None else "no"

            factor_is_na    = (factor_str == "NA")
            factor_matched  = (not factor_is_na and
                               factor_matches_patient(
                                   factor_str, stage_arr, grade_arr, menopause_arr))

            # ══════════════════════════════════════════════════════════
            # CP1: In literature (NA factor) + In HPA
            # General regulation — no stage/grade specificity
            # ══════════════════════════════════════════════════════════
            if (index_protein is not None and index_normal is not None
                    and factor_is_na):
                generate_abundance[i] = assign_abundance(
                    regulation, normal_reg, prognostic)

            # ══════════════════════════════════════════════════════════
            # CP2: NOT in literature + In HPA
            # Use HPA baseline direction as surrogate regulation signal
            # ══════════════════════════════════════════════════════════
            elif index_protein is None and index_normal is not None:
                inferred_reg = ("up"   if normal_reg == "up"   else
                               "down"  if normal_reg == "down" else "up/down")
                generate_abundance[i] = assign_abundance(
                    inferred_reg, normal_reg, prognostic)

            # ══════════════════════════════════════════════════════════
            # CP3: In literature (NA factor) + NOT in HPA
            # BUG FIX 3: was gene_prognosis_indicator[None] → crash
            # Now uses safe fallback: prognostic='no', normal_reg='moderate'
            # ══════════════════════════════════════════════════════════
            elif (index_protein is not None and index_normal is None
                      and factor_is_na):
                generate_abundance[i] = assign_abundance(
                    regulation, "moderate", "no")

            # ══════════════════════════════════════════════════════════
            # CP4: In literature (factor-specific) + NOT in HPA
            # Fuzzy rule: activate strong perturbation only when patient's
            # stage/grade/menopause matches the protein's factor condition

            elif (index_protein is not None and index_normal is None
                      and not factor_is_na):
                if factor_matched:
                    generate_abundance[i] = assign_abundance(
                        regulation, "moderate", "no")
                else:
                    generate_abundance[i] = sample_neutral()

            # ══════════════════════════════════════════════════════════
            # CP5: Neither literature nor HPA
            # Baseline noise — protein biology unknown
            # ══════════════════════════════════════════════════════════
            elif index_protein is None and index_normal is None:
                generate_abundance[i] = sample_neutral()

            # ══════════════════════════════════════════════════════════
            # CP6: In literature (factor-specific) + In HPA
            # Full fuzzy evaluation: regulation + prognostics + HPA baseline
            # + stage/grade factor gate
            # BUG FIX 1: was .splt() — AttributeError (12 occurrences here)
            # BUG FIX 4: was list1+list2+list3 on numpy arrays
            # ══════════════════════════════════════════════════════════
            elif (index_protein is not None and index_normal is not None
                      and not factor_is_na):
                if factor_matched:
                    generate_abundance[i] = assign_abundance(
                        regulation, normal_reg, prognostic)
                else:
                    generate_abundance[i] = sample_neutral()

            else:
                generate_abundance[i] = sample_neutral()

        all_outputs.append((f"Patient{j:04d}", generate_abundance.copy()))

    return all_outputs


# ══════════════════════════════════════════════════════════════════════
# SECTION 3 — USE CASE 1: Random Forest Classifier
# Paper: "training a Random Forest classifier for endometrial cancer
# stage prediction achieving 82% accuracy on synthetic data"
# ADDED: This function was described in the paper but missing from code
# ══════════════════════════════════════════════════════════════════════
def use_case_1_random_forest():
    """
    Demonstrates Use Case 1 from paper Section 3:
    Train a Random Forest classifier on synthetic clinical data
    to predict early vs late stage endometrial cancer.
    """
    if clinical_data.empty:
        return "Please generate clinical data first (Clinical Data tab)."

    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import LabelEncoder
        from sklearn.metrics import accuracy_score, classification_report

        df = clinical_data.copy()

        # Binary target: early (Stage1/2) vs late (Stage3/4)
        df['Stage_binary'] = df['Stage'].apply(
            lambda x: 0 if x.startswith(('Stage1','Stage2')) else 1
        )

        # Encode categorical features
        le = LabelEncoder()
        features = ['Ages','BMI','Menopause','Nulliparity','Tumor_type','Grade']
        X = df[features].copy()
        for col in ['Nulliparity','Tumor_type','Grade']:
            X[col] = le.fit_transform(X[col].astype(str))

        y = df['Stage_binary']

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42)

        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        acc = accuracy_score(y_test, y_pred)

        report = classification_report(y_test, y_pred,
                                       target_names=['Early Stage','Late Stage'])
        result = f"Accuracy: {acc*100:.1f}%\n\n{report}"
        return result

    except ImportError:
        return "scikit-learn not installed. Run: pip install scikit-learn"
    except Exception as e:
        return f"Error: {str(e)}"


# ══════════════════════════════════════════════════════════════════════
# WRAPPER FOR GRADIO PROTEIN TAB
# ══════════════════════════════════════════════════════════════════════
def wrapper(num_records, text_input):
    try:
        lst = ast.literal_eval(text_input)
        if isinstance(lst, str):   lst = [lst]
        if isinstance(lst, tuple): lst = list(lst)
        if not isinstance(lst, list): lst = [lst]
    except Exception:
        lst = [x.strip() for x in text_input.split(",") if x.strip()]
    return generate_protein_abundance_data(num_records, protein_list=lst)


# ══════════════════════════════════════════════════════════════════════
# CSV DOWNLOAD HELPER
# ADDED per paper Section 2.3: "outputs downloadable as CSV files"
# Original code had no download button in the Gradio interface
# ══════════════════════════════════════════════════════════════════════
def download_clinical_csv():
    """Export synthetic clinical data as downloadable CSV."""
    if clinical_data.empty:
        return None
    path = "/tmp/SynthProteomics_clinical.csv"
    clinical_data.to_csv(path, index=False)
    return path


def download_protein_csv(num_records, text_input):
    """Export synthetic protein abundance data as downloadable CSV."""
    results = wrapper(num_records, text_input)
    if not results:
        return None
    rows = []
    for patient_id, abundances in results:
        row = {'Patient_ID': patient_id}
        row.update(abundances)
        rows.append(row)
    df_out = pd.DataFrame(rows)
    path = "/tmp/SynthProteomics_protein.csv"
    df_out.to_csv(path, index=False)
    return path


# ══════════════════════════════════════════════════════════════════════
# DISTRIBUTION PLOTS
# BUG FIX 5: original plots() passed column name string to hist()
# e.g. axes[c,index].hist('Ages', data=clinical_data) — this doesn't
# work correctly. Now explicitly passes clinical_data[col] array.
# ══════════════════════════════════════════════════════════════════════
def plots():
    if clinical_data.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "Generate clinical data first",
                ha='center', va='center', fontsize=12)
        ax.axis('off')
        return fig

    plot_cols = ['Ages', 'BMI', 'Grade', 'Stage', 'Tumor_type',
                 'Menopause', 'Nulliparity', 'Treatment', 'Survival_Outcome']
    plot_cols = [c for c in plot_cols if c in clinical_data.columns]

    n    = len(plot_cols)
    rows = (n + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(14, rows*3), squeeze=False)
    fig.suptitle("SynthProteomics — Synthetic Clinical Variable Distributions",
                 fontsize=13, fontweight='bold', y=1.01)

    for idx, col in enumerate(plot_cols):
        r, c = divmod(idx, 2)
        ax = axes[r][c]
        # BUG FIX 5: was axes[c,index].hist(i,...,data=clinical_data)
        # which passes the column name string, not the data array
        if clinical_data[col].dtype in [np.float64, np.int64, float, int]:
            ax.hist(clinical_data[col].values, bins=15,
                    color='steelblue', edgecolor='white', alpha=0.85)
        else:
            counts = clinical_data[col].value_counts()
            ax.bar(counts.index.astype(str), counts.values,
                   color='steelblue', edgecolor='white', alpha=0.85)
            ax.tick_params(axis='x', rotation=35)
        ax.set_title(col, fontweight='bold', fontsize=10)
        ax.set_ylabel("Count", fontsize=9)
        ax.grid(axis='y', alpha=0.3)

    # Hide unused subplot panels
    for idx in range(n, rows * 2):
        r, c = divmod(idx, 2)
        axes[r][c].set_visible(False)

    plt.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════
# SECTION 2.3 — GRADIO INTERFACE
# Paper: "browser-accessible Gradio interface requiring no programming
# experience; users specify cohort size, cancer subtype, stage
# distribution, desired proteomic panel; outputs as CSV"
#
# CHANGES FROM ORIGINAL:
# + Added subtype selector (Endometrioid/Serous/Clear Cell/All)
# + Added stage filter (Stage1/Stage2/Stage3/Stage4/All)
# + Added CSV download buttons (completely missing from original)
# + Added Use Case 1 RF classifier tab (described in paper but absent)
# + Added descriptive tab labels and header markdown
# ══════════════════════════════════════════════════════════════════════
with gr.Blocks(title="SynthProteomics") as demo:

    gr.Markdown("""
    # SynthProteomics
    **Probabilistic Rule-Based Synthetic Clinical & Proteomic Data Generator**
    *Endometrial Cancer Research | University of Newcastle*
    GitHub: https://github.com/NehaAr/Synthetic-omics-data | MIT License
    """)

    # ── Tab 1: Clinical Data ──────────────────────────────────────────
    with gr.Tab("Clinical Data"):
        gr.Markdown("### Section 2.1 — Probabilistic Rule-Based Clinical Simulation")
        with gr.Row():
            n_records_c = gr.Number(label="Number of Patients", value=100, minimum=10)
            # ADDED: subtype selector per paper Section 2.3
            subtype_sel = gr.Dropdown(
                choices=["All","Endometrioid","Serous","Clear Cell",
                         "Mucinous","Undifferentiated"],
                value="All", label="Cancer Subtype Filter"
            )
            # ADDED: stage filter per paper Section 2.3
            stage_sel = gr.Dropdown(
                choices=["All","Stage1","Stage2","Stage3","Stage4"],
                value="All", label="Stage Filter"
            )
        col_selector = gr.CheckboxGroup(
            choices=["Patient_ID","Ages","Ethnicity","Menopause","Grade",
                     "Tumor_type","Histological_Subtype","Stage","Nulliparity",
                     "BMI","Myometrial_mm","Treatment","Survival_Outcome"],
            value=["Patient_ID","Ages","Grade","Stage","Tumor_type",
                   "Histological_Subtype","Survival_Outcome"],
            label="Select Columns to Display"
        )
        btn_clin   = gr.Button("Generate Clinical Data", variant="primary")
        out_clin   = gr.Dataframe(label="Synthetic Clinical Data")
        # ADDED: CSV download per paper Section 2.3
        btn_dl_clin = gr.Button("Download as CSV")
        file_clin   = gr.File(label="Download Clinical CSV")

        btn_clin.click(
            generate_person_data,
            inputs=[n_records_c, col_selector, subtype_sel, stage_sel],
            outputs=out_clin
        )
        btn_dl_clin.click(download_clinical_csv, inputs=[], outputs=file_clin)

    # ── Tab 2: Protein Abundance ──────────────────────────────────────
    with gr.Tab("Protein Abundance Data"):
        gr.Markdown("### Section 2.2 — Fuzzy Rule-Based Proteomic Simulation")
        gr.Markdown(
            "_Generate clinical data first. Protein abundances are gated by "
            "patient stage/grade/menopausal status._"
        )
        n_records_p  = gr.Number(label="Number of Patients", value=10, minimum=1)
        protein_input = gr.Textbox(
            lines=4,
            label="Protein List (comma-separated)",
            placeholder="e.g. ANXA2, PKM2, ERBB2, EGFR, MMP9"
        )
        btn_prot  = gr.Button("Generate Abundance Data", variant="primary")
        out_prot  = gr.JSON(label="Log2 Fold-Change Abundance Values [-3 to 3]")
        # ADDED: CSV download per paper Section 2.3
        btn_dl_prot = gr.Button("Download as CSV")
        file_prot   = gr.File(label="Download Protein CSV")

        btn_prot.click(wrapper, inputs=[n_records_p, protein_input], outputs=out_prot)
        btn_dl_prot.click(
            download_protein_csv,
            inputs=[n_records_p, protein_input],
            outputs=file_prot
        )

    # ── Tab 3: Plots ──────────────────────────────────────────────────
    with gr.Tab("Distribution Plots"):
        gr.Markdown("### Clinical Variable Distributions")
        btn_plot = gr.Button("Generate Plots", variant="primary")
        out_plot = gr.Plot(label="Synthetic Clinical Data Distributions")
        btn_plot.click(plots, inputs=[], outputs=out_plot)

    # ── Tab 4: Use Case 1 — RF Classifier ────────────────────────────
    # ADDED: Use Case 1 described in paper Section 3 was missing from code
    with gr.Tab("Use Case: Stage Classifier"):
        gr.Markdown("""
        ### Section 3 — Use Case 1: Random Forest Stage Classifier
        *Trains a Random Forest on synthetic data to predict early vs late stage.*
        Generate clinical data first, then run the classifier.
        """)
        btn_rf  = gr.Button("Run Random Forest Classifier", variant="primary")
        out_rf  = gr.Textbox(label="Classification Results", lines=15)
        btn_rf.click(use_case_1_random_forest, inputs=[], outputs=out_rf)

demo.launch()
