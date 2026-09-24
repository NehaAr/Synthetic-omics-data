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
subprocess.run([sys.executable, "-m", "pip", "install", "gradio", "scikit-learn", "pandas", "numpy", "matplotlib", "--quiet"], check=False)

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
        'NA','NA','NA','NA','NA','NA','TYPE 2','TYPE 1','NA','NA','NA','NA',
        'GRADE1,grade3,stage1,stage3','GRADE1,GRADE3','NA','TYPE1','TYPE1','TYPE1',
        'STAGE1A,STAGE1B','NA','NA','NA','NA','type 1','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','stage1B','NA','NA','NA','stage1A','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA',
        'NA','NA','NA','NA','NA','NA','NA'
    ]
}

# Annotated Pathways dictionary for Section 2.2 covariance blending
PATHWAY_MAP = {
    # Original Pathways
    'PI3K_AKT': ['CTNNB', 'EGFR', 'ERBB2', 'PKM2', 'MMP9', 'NFKB', 'PTEN'],
    'APOPTOSIS': ['CASP3', 'HMGB3', 'PRDX1', 'PRDX6', 'SOD1'],
    'GLYCOLYSIS': ['ENO1', 'ALDOA', 'LDHA', 'PKM', 'TPI1', 'PGAM2'],
    
    # Structural & Hormonal Pathways
    'WNT_BETA_CATENIN': ['CTNNB1', 'APC', 'AXIN1', 'GSK3B'],
    'MAPK_ERK': ['KRAS', 'BRAF', 'MAPK1', 'EGFR', 'ERBB2'],
    'DNA_REPAIR': ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'POLE', 'TP53'],
    
    # Genomic & Microenvironment Pathways
    'CHROMATIN_REMODELING': ['ARID1A', 'SMARCA4', 'ARID1B'],
    'TGF_BETA_EMT': ['TGFB1', 'SMAD2', 'SMAD3', 'SMAD4'],
    'ANGIOGENESIS': ['VEGFA', 'KDR', 'HIF1A', 'ANGPT2']
}

# ══════════════════════════════════════════════════════════════════════
# GLOBAL STATE
# ══════════════════════════════════════════════════════════════════════
clinical_data = pd.DataFrame()
protein_data = pd.DataFrame()

try:
    normal_tissue_expression = pd.read_csv('data/HPA.tsv', sep='\t')
    normal_tissue_prognostic = (
        normal_tissue_expression['Cancer prognostics - Uterine Corpus Endometrial Carcinoma (TCGA)']
        .astype(str).str.strip().str.replace(r'[^\d.]', '', regex=True)
    )
    normal_tissue_prognostic = pd.to_numeric(normal_tissue_prognostic, errors='coerce').astype(float)
    normal_tisue_regulation = normal_tissue_expression['Tissue RNA - endometrium 1 [nTPM]'].apply(
        lambda x: "up" if x > 100 else ("moderate" if 10 <= x <= 100 else "down")
    )
    gene_prognosis_indicator = normal_tissue_prognostic.apply(
        lambda x: "yes" if pd.notna(x) and x <= 0.05 else "no"
    )
    gene_list_lower = normal_tissue_expression['Gene'].str.lower().tolist()
except FileNotFoundError:
    normal_tissue_expression = pd.DataFrame(columns=['Gene'])
    normal_tisue_regulation = pd.Series(dtype=str)
    gene_prognosis_indicator = pd.Series(dtype=str)
    gene_list_lower = []

# ══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS & MAMDANI FUZZY ENGINE
# ══════════════════════════════════════════════════════════════════════

def compute_mamdani_scale(grade_str, stage_str):
    """
    Section 2.2 Mamdani Fuzzy Inference Engine:
    Maps Grade (1-3) and FIGO Stage index (1-10) to a scale factor kappa in [0.5, 1.5].
    """
    g_val = 1.0 if '1' in str(grade_str) else (2.0 if '2' in str(grade_str) else 3.0)
    stage_map = {'Stage1':1, 'Stage1A':1, 'Stage1B':2, 'Stage2':3, 'Stage2A':3, 'Stage2B':4, 
                 'Stage3C1':6, 'Stage3C2':7, 'Stage4A':9, 'Stage4B':10}
    s_val = float(stage_map.get(str(stage_str), 5))

    # Normalized inputs
    g_norm = (g_val - 1.0) / 2.0  # [0, 1]
    s_norm = (s_val - 1.0) / 9.0  # [0, 1]

    # Membership degrees
    g_high = g_norm
    g_low = 1.0 - g_norm
    s_high = s_norm
    s_low = 1.0 - s_norm

    # Rules
    r_high = min(g_high, s_high)
    r_low = min(g_low, s_low)
    r_med = 1.0 - max(r_high, r_low)

    # Defuzzification via Centroid
    kappa = (r_low * 0.7 + r_med * 1.0 + r_high * 1.3) / (r_low + r_med + r_high + 1e-6)
    return float(np.clip(kappa, 0.5, 1.5))

def parse_factors(factor_string):
    if not factor_string or str(factor_string).upper() == "NA":
        return []
    parts = re.split(r'[/,]', factor_string)
    return [p.strip().lower() for p in parts if p.strip()]

def factor_matches_patient(factor_string, stage, grade, menopause):
    factors = parse_factors(factor_string)
    if not factors:
        return False
    patient_vals = [str(stage).lower(), str(grade).lower(), str(menopause).lower()]
    return any(f in patient_vals for f in factors)

def sample_up_strong(): return float(np.clip(np.random.normal(loc=2.0, scale=0.4), 0.0, 3.0))
def sample_up_weak(): return float(np.clip(np.random.normal(loc=1.0, scale=0.3), 0.0, 3.0))
def sample_down_strong(): return float(np.clip(np.random.normal(loc=-2.0, scale=0.4), -3.0, 0.0))
def sample_down_weak(): return float(np.clip(np.random.normal(loc=-1.0, scale=0.3), -3.0, 0.0))
def sample_low(): return float(np.clip(np.random.normal(loc=-0.5, scale=0.3), -1.5, 0.0))
def sample_neutral(): return float(np.clip(np.random.normal(loc=0.0, scale=0.5), -1.5, 1.5))

def assign_abundance(regulation, normal_reg, prognostic):
    reg, prog = str(regulation).lower().strip(), str(prognostic).lower().strip()
    if reg == 'up':
        return sample_up_strong() if prog == 'yes' else sample_up_weak()
    elif reg == 'down':
        return sample_down_strong() if prog == 'yes' else sample_down_weak()
    elif 'up' in reg and 'down' in reg:
        if normal_reg == 'up': return sample_up_strong() if prog == 'yes' else sample_up_weak()
        elif normal_reg == 'down': return sample_down_strong() if prog == 'yes' else sample_down_weak()
        else: return sample_neutral()
    elif reg == 'low':
        return sample_low()
    else:
        return sample_neutral()

# ══════════════════════════════════════════════════════════════════════
# DATA GENERATION PIPELINES
# ══════════════════════════════════════════════════════════════════════

def generate_person_data(num_records, selected_columns, subtype_filter="All", stage_filter="All"):
    global clinical_data
    num_records = int(num_records)
    np.random.seed(42)
    random.seed(42)

    Ages = np.clip(np.random.normal(loc=60, scale=15, size=num_records).astype(int), 30, 85)
    bmi_values = np.round(np.clip(np.random.normal(loc=25, scale=4, size=num_records), 21, 40), 2)
    ethnicities = random.choices(['European','East Asian','South Asian','Middle Eastern','African','Indigenous Australian','Other'],
                                 weights=[0.55, 0.15, 0.10, 0.05, 0.03, 0.03, 0.09], k=num_records)
    
    patient_id, Treatment, Menopause, Grade, Stage, Myometrial_invasion = [], [], [], [], [], []
    Nulliparity, Tumor_type, Subtype, Survival_outcome = [], [], [], []

    stage_choices = ['Stage1','Stage1A','Stage1B','Stage2','Stage2A','Stage2B','Stage3C1','Stage3C2','Stage4A','Stage4B']
    early_weights = [0.1333]*6 + [0.05]*4
    late_weights = [0.05]*6 + [0.1333]*4

    for i in range(num_records):
        patient_id.append(f"EC{i:04d}")
        Treatment.append(random.choice(["Surgery", "Chemotherapy", "Radiotherapy", "Combined"]))
        Myometrial_invasion.append(round(random.uniform(2, 40), 2))
        Menopause.append(1 if Ages[i] <= 51 else 2)
        Nulliparity.append(random.choices(["yes", "no"], [0.5, 0.5] if Ages[i] <= 35 else [0.6, 0.4])[0])

        high_risk = (Ages[i] >= 60 and Nulliparity[i] == 'yes' and bmi_values[i] >= 25)
        t_type = random.choices(['Type1','Type2'], [0.4, 0.6] if high_risk else [0.8, 0.2])[0]
        Tumor_type.append(t_type)

        if t_type == 'Type1':
            g_val = random.choices(['Grade1','Grade2','Grade3'], [0.5, 0.3, 0.2])[0]
            sub = random.choices(['Endometrioid','Mucinous'], [0.90, 0.10])[0]
        else:
            g_val = random.choices(['Grade1','Grade2','Grade3'], [0.2, 0.3, 0.5])[0]
            sub = random.choices(['Serous','Clear Cell','Undifferentiated'], [0.55, 0.35, 0.10])[0]
        
        Grade.append(g_val)
        Subtype.append(sub)

        s_val = random.choices(stage_choices, weights=early_weights if g_val in ('Grade1', 'Grade2') else late_weights)[0]
        Stage.append(s_val)

        if s_val.startswith('Stage1') and g_val == 'Grade1':
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.90, 0.10])[0])
        elif s_val.startswith('Stage4'):
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.30, 0.70])[0])
        else:
            Survival_outcome.append(random.choices(['Alive','Deceased'], [0.65, 0.35])[0])

    clinical_data = pd.DataFrame({
        'Patient_ID': patient_id, 'Ages': Ages, 'Ethnicity': ethnicities, 'Menopause': Menopause,
        'Grade': Grade, 'Tumor_type': Tumor_type, 'Histological_Subtype': Subtype, 'Stage': Stage,
        'Nulliparity': Nulliparity, 'BMI': bmi_values, 'Myometrial_mm': Myometrial_invasion,
        'Treatment': Treatment, 'Survival_Outcome': Survival_outcome
    })

    if subtype_filter != "All":
        clinical_data = clinical_data[clinical_data['Histological_Subtype'] == subtype_filter].reset_index(drop=True)
    if stage_filter != "All":
        clinical_data = clinical_data[clinical_data['Stage'].str.startswith(stage_filter)].reset_index(drop=True)

    clinical_data.to_csv('clinical_data.csv', index=False)
    valid_cols = [c for c in selected_columns if c in clinical_data.columns]
    return clinical_data[valid_cols] if valid_cols else clinical_data


def generate_protein_abundance_data(num_records, protein_list, rho=0.4):
    """
    Algorithm 1: Fuzzy-Rule & Pathway Covariance Blended Simulation
    """
    global protein_data, clinical_data
    num_records = int(num_records)
    if clinical_data.empty or len(clinical_data) < num_records:
        generate_person_data(num_records, ["Patient_ID", "Stage", "Grade", "Menopause"])

    all_outputs = []
    gene_dict_lower = [g.lower() for g in protein_abundance_dictionary['Gene']]

    for j in range(num_records):
        patient_row = clinical_data.iloc[j] if j < len(clinical_data) else clinical_data.iloc[0]
        p_stage, p_grade, p_meno = patient_row['Stage'], patient_row['Grade'], patient_row['Menopause']
        
        # Step 1: Calculate Mamdani Fuzzy Scale Factor kappa
        kappa = compute_mamdani_scale(p_grade, p_stage)
        
        # Step 2: Draw pathway latent factors z_k ~ N(0, 1)
        pathway_z = {pw: np.random.normal(0, 1) for pw in PATHWAY_MAP}

        patient_abundances = {}
        for prot in protein_list:
            prot_clean = prot.strip()
            prot_lower = prot_clean.lower()

            idx_prot = gene_dict_lower.index(prot_lower) if prot_lower in gene_dict_lower else None
            idx_norm = gene_list_lower.index(prot_lower) if prot_lower in gene_list_lower else None

            reg = protein_abundance_dictionary['Regulation'][idx_prot] if idx_prot is not None else "UP/DOWN"
            factor_str = protein_abundance_dictionary['Factor'][idx_prot] if idx_prot is not None else "NA"
            norm_reg = normal_tisue_regulation.iloc[idx_norm] if idx_norm is not None else "moderate"
            prog = gene_prognosis_indicator.iloc[idx_norm] if idx_norm is not None else "no"

            factor_matched = (factor_str == "NA") or factor_matches_patient(factor_str, p_stage, p_grade, p_meno)
            base_log2fc = assign_abundance(reg, norm_reg, prog) if factor_matched else sample_neutral()

            # Step 3: Pathway Covariance Blending
            shared_pw = [pw for pw, genes in PATHWAY_MAP.items() if prot_clean in genes]
            if shared_pw:
                z_k = pathway_z[shared_pw[0]]
                epsilon = np.random.normal(0, 1)
                final_log2fc = base_log2fc + (rho * z_k + np.sqrt(1.0 - rho**2) * epsilon)
            else:
                final_log2fc = base_log2fc

            # Step 4: Scale by kappa and clip to [-3, 3]
            patient_abundances[prot_clean] = float(np.clip(kappa * final_log2fc, -3.0, 3.0))

        all_outputs.append((f"Patient{j:04d}", patient_abundances))

    # Update global dataframe
    rows = []
    for pid, ab in all_outputs:
        r = {'Patient_ID': pid}
        r.update(ab)
        rows.append(r)
    protein_data = pd.DataFrame(rows)
    return all_outputs


def use_case_1_random_forest():
    """
    Section 3 Use Case 1: Random Forest Classifier (Early vs Late Stage)
    Joint multi-modal classification combining clinical and proteomic features.
    """
    global clinical_data, protein_data
    if clinical_data.empty:
        generate_person_data(1000, ["Patient_ID","Ages","Ethnicity","Menopause","Grade",
                                   "Tumor_type","Histological_Subtype","Stage","Nulliparity",
                                   "BMI","Myometrial_mm","Treatment","Survival_Outcome"])

    default_panel = ['CTNNB', 'EGFR', 'ERBB2', 'PKM2', 'MMP9', 'CASP3', 'HMGB3', 'PRDX1', 'PTEN', 'TP53']
    if protein_data.empty or len(protein_data) != len(clinical_data):
        generate_protein_abundance_data(len(clinical_data), default_panel)

    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import LabelEncoder
        from sklearn.metrics import accuracy_score, classification_report

        df_clin = clinical_data.copy().reset_index(drop=True)
        df_prot = protein_data.copy().reset_index(drop=True)
        
        # Target: Early (Stage 1/2) vs Late (Stage 3/4)
        df_clin['Stage_binary'] = df_clin['Stage'].apply(
            lambda x: 0 if str(x).startswith(('Stage1', 'Stage2')) else 1
        )

        # Categorical Encoding
        le = LabelEncoder()
        cat_cols = ['Ethnicity', 'Menopause', 'Grade', 'Tumor_type', 'Histological_Subtype', 'Nulliparity', 'Treatment', 'Survival_Outcome']
        for col in cat_cols:
            if col in df_clin.columns:
                df_clin[col] = le.fit_transform(df_clin[col].astype(str))

        # Merge Clinical + Proteomic Features
        feature_cols_clin = ['Ages', 'Ethnicity', 'Menopause', 'Grade', 'Tumor_type', 'Histological_Subtype', 
                             'Nulliparity', 'BMI', 'Myometrial_mm', 'Treatment', 'Survival_Outcome']
        
        X_clin = df_clin[[c for c in feature_cols_clin if c in df_clin.columns]]
        X_prot = df_prot.drop(columns=['Patient_ID'], errors='ignore')
        
        X = pd.concat([X_clin, X_prot], axis=1)
        y = df_clin['Stage_binary']

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        
        acc = accuracy_score(y_test, y_pred)
        report = classification_report(y_test, y_pred, target_names=['Early Stage (I/II)', 'Late Stage (III/IV)'])
        
        return f"Random Forest Classification Accuracy: {acc*100:.1f}%\n\nClassification Report:\n{report}"

    except Exception as e:
        return f"Error executing classifier: {str(e)}"


def wrapper(num_records, text_input):
    try:
        lst = ast.literal_eval(text_input)
        if isinstance(lst, (str, tuple)): lst = list(lst) if isinstance(lst, tuple) else [lst]
        if not isinstance(lst, list): lst = [lst]
    except Exception:
        lst = [x.strip() for x in text_input.split(",") if x.strip()]
    return generate_protein_abundance_data(num_records, protein_list=lst)

def download_clinical_csv():
    if clinical_data.empty: return None
    path = "/tmp/SynthProteomics_clinical.csv"
    clinical_data.to_csv(path, index=False)
    return path

def download_protein_csv(num_records, text_input):
    results = wrapper(num_records, text_input)
    if not results: return None
    path = "/tmp/SynthProteomics_protein.csv"
    protein_data.to_csv(path, index=False)
    return path

def plots():
    if clinical_data.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "Generate clinical data first", ha='center', va='center', fontsize=12)
        ax.axis('off')
        return fig

    plot_cols = ['Ages', 'BMI', 'Grade', 'Stage', 'Tumor_type', 'Menopause', 'Nulliparity', 'Treatment', 'Survival_Outcome']
    plot_cols = [c for c in plot_cols if c in clinical_data.columns]

    n = len(plot_cols)
    rows = (n + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(14, rows*3), squeeze=False)
    fig.suptitle("SynthProteomics — Synthetic Clinical Variable Distributions", fontsize=13, fontweight='bold', y=1.01)

    for idx, col in enumerate(plot_cols):
        r, c = divmod(idx, 2)
        ax = axes[r][c]
        if clinical_data[col].dtype in [np.float64, np.int64, float, int]:
            ax.hist(clinical_data[col].values, bins=15, color='steelblue', edgecolor='white', alpha=0.85)
        else:
            counts = clinical_data[col].value_counts()
            ax.bar(counts.index.astype(str), counts.values, color='steelblue', edgecolor='white', alpha=0.85)
            ax.tick_params(axis='x', rotation=35)
        ax.set_title(col, fontweight='bold', fontsize=10)
        ax.set_ylabel("Count", fontsize=9)
        ax.grid(axis='y', alpha=0.3)

    for idx in range(n, rows * 2):
        r, c = divmod(idx, 2)
        axes[r][c].set_visible(False)

    plt.tight_layout()
    return fig

# ══════════════════════════════════════════════════════════════════════
# GRADIO INTERFACE
# ══════════════════════════════════════════════════════════════════════

with gr.Blocks(title="SynthProteomics") as demo:
    gr.Markdown("""
    # SynthProteomics
    **Probabilistic Rule-Based Synthetic Clinical & Proteomic Data Generator**
    *Endometrial Cancer Research | University of Newcastle*
    GitHub: https://github.com/NehaAr/Synthetic-omics-data | MIT License
    """)

    with gr.Tab("Clinical Data"):
        gr.Markdown("### Section 2.1 — Probabilistic Rule-Based Clinical Simulation")
        with gr.Row():
            n_records_c = gr.Number(label="Number of Patients", value=100, minimum=10)
            subtype_sel = gr.Dropdown(choices=["All","Endometrioid","Serous","Clear Cell","Mucinous","Undifferentiated"], value="All", label="Cancer Subtype Filter")
            stage_sel = gr.Dropdown(choices=["All","Stage1","Stage2","Stage3","Stage4"], value="All", label="Stage Filter")
        col_selector = gr.CheckboxGroup(
            choices=["Patient_ID","Ages","Ethnicity","Menopause","Grade","Tumor_type","Histological_Subtype","Stage","Nulliparity","BMI","Myometrial_mm","Treatment","Survival_Outcome"],
            value=["Patient_ID","Ages","Grade","Stage","Tumor_type","Histological_Subtype","Survival_Outcome"],
            label="Select Columns to Display"
        )
        btn_clin = gr.Button("Generate Clinical Data", variant="primary")
        out_clin = gr.Dataframe(label="Synthetic Clinical Data")
        btn_dl_clin = gr.Button("Download as CSV")
        file_clin = gr.File(label="Download Clinical CSV")

        btn_clin.click(generate_person_data, inputs=[n_records_c, col_selector, subtype_sel, stage_sel], outputs=out_clin)
        btn_dl_clin.click(download_clinical_csv, inputs=[], outputs=file_clin)

    with gr.Tab("Protein Abundance Data"):
        gr.Markdown("### Section 2.2 — Fuzzy Rule & Pathway Co-Regulation Proteomic Simulation")
        gr.Markdown("_Generate clinical data first. Protein abundances are modulated via Mamdani inference and pathway covariance._")
        n_records_p = gr.Number(label="Number of Patients", value=10, minimum=1)
        protein_input = gr.Textbox(lines=4, label="Protein List (comma-separated)", value="ANXA2, PKM2, ERBB2, EGFR, MMP9, CASP3, HMGB3, PRDX1")
        btn_prot = gr.Button("Generate Abundance Data", variant="primary")
        out_prot = gr.JSON(label="Log2 Fold-Change Abundance Values [-3 to 3]")
        btn_dl_prot = gr.Button("Download as CSV")
        file_prot = gr.File(label="Download Protein CSV")

        btn_prot.click(wrapper, inputs=[n_records_p, protein_input], outputs=out_prot)
        btn_dl_prot.click(download_protein_csv, inputs=[n_records_p, protein_input], outputs=file_prot)

    with gr.Tab("Distribution Plots"):
        gr.Markdown("### Clinical Variable Distributions")
        btn_plot = gr.Button("Generate Plots", variant="primary")
        out_plot = gr.Plot(label="Synthetic Clinical Data Distributions")
        btn_plot.click(plots, inputs=[], outputs=out_plot)

    with gr.Tab("Use Case: Stage Classifier"):
        gr.Markdown("### Section 3 — Use Case 1: Random Forest Classifier")
        gr.Markdown("_Train a Random Forest classifier (100 estimators) on combined clinical and proteomic features to predict early vs late stage._")
        btn_rf = gr.Button("Run Random Forest Evaluation", variant="primary")
        out_rf = gr.Textbox(lines=10, label="Classification Results")
        btn_rf.click(use_case_1_random_forest, inputs=[], outputs=out_rf)

if __name__ == "__main__":
    demo.launch()
