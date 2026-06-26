# Algorithm 1 — SynthProteomics: Probabilistic Fuzzy Rule-Based Synthetic Data Generation

## Overview

SynthProteomics generates synthetic data in two sequential parts:

1. **Clinical Data Generation** — probabilistic rule-based simulation of patient attributes
2. **Protein Abundance Generation** — fuzzy inference over HPA-validated expression profiles

---

## Inputs & Outputs

| Symbol | Description |
|---|---|
| `N` | Number of synthetic patients to generate |
| `P` | List of protein gene symbols (user-supplied) |
| `HPA` | Human Protein Atlas expression database |
| `D` | Protein regulation dictionary (207 genes; regulation direction + factor) |
| `C` | Output: synthetic clinical dataset |
| `A` | Output: protein abundance matrix (log₂ fold-change values) |

---

## Part 1 — Clinical Data Generation

```
for each patient i = 1 to N:

    1. Sample continuous variables
       Age_i     ~ Normal(μ=60, σ=15)   clipped to [30, 85]
       BMI_i     ~ Normal(μ=25, σ=4)    clipped to [21, 40]

    2. Derive binary variables from Age_i
       Menopause_i  = TRUE   if Age_i >= 50,  else FALSE
       Nulliparity_i = TRUE  if Age_i <= 35,  else FALSE  (with probabilistic weight)

    3. Assign TumorType_i
       ~ Categorical( weights = f(Age_i, BMI_i, Nulliparity_i) )
       Classes: endometrioid | serous | clear cell | mucinous | undifferentiated

    4. Assign Grade_i
       ~ Categorical( weights = f(TumorType_i) )
       Classes: Grade 1 | Grade 2 | Grade 3

    5. Assign FIGO_Stage_i
       ~ Categorical( weights = f(Grade_i) )
       Substages: IA | IB | II | IIIA | IIIB | IIIC1 | IIIC2 | IVA | IVB
       (calibrated to AIHW 2023 Australian cancer registry data)

    6. Assign Treatment_i
       ~ Categorical( weights = f(FIGO_Stage_i) )
       Classes: surgery | chemo | radiation | combined

end for

C ← { Age, BMI, Menopause, Nulliparity, TumorType, Grade, FIGO_Stage, Treatment }
```

---

## Part 2 — Protein Abundance Generation

```
for each patient j = 1 to N:
  for each protein p in P:

    1. Retrieve from regulation dictionary D:
       r_p  = regulation direction  { UP | DOWN | UP-DOWN | LOW }
       f_p  = perturbation factor

    2. Retrieve from HPA database:
       e_p  = baseline expression level in endometrial tissue
       s_p  = prognostic significance flag  { yes | no }
              (yes if HPA p ≤ 0.05 from TCGA survival analysis)

    3. Apply fuzzy rule → sample log2FC:

       ┌─────────────────────────────────────────────────────────────────┐
       │  IF r_p = UP      AND s_p = yes  → log2FC ~ Normal(+2.0, 0.4)  │
       │  IF r_p = UP      AND s_p = no   → log2FC ~ Normal(+1.0, 0.3)  │
       │  IF r_p = DOWN    AND s_p = yes  → log2FC ~ Normal(−2.0, 0.4)  │
       │  IF r_p = DOWN    AND s_p = no   → log2FC ~ Normal(−1.0, 0.3)  │
       │  IF r_p = UP-DOWN AND s_p = yes  → log2FC ~ Normal(±e_p×2, 0.4)│
       │  IF r_p = UP-DOWN AND s_p = no   → log2FC ~ Normal(±e_p×1, 0.3)│
       │  OTHERWISE                        → log2FC ~ Normal(0, 0.5)     │
       └─────────────────────────────────────────────────────────────────┘

    4. Apply clinical fuzzy layer
       Scale log2FC_p,j by tumour grade and FIGO stage of patient j
       (higher grade / later stage → amplified perturbation magnitude)

    5. Clip log2FC_p,j to [−3, +3]

  end for
end for

A ← { log2FC_p,j }  for all patients j and proteins p
```

---

## Fuzzy Rule Logic — Explained

The fuzzy layer is inspired by **Mamdani fuzzy inference** (Zadeh, 1965). Rather than a hard threshold (e.g., "if Stage III then high expression"), it maps linguistic input variables to graded, continuous output values sampled from Gaussian distributions.

```
Linguistic Inputs:
  - Regulation direction  →  { UP, DOWN, UP-DOWN, LOW }
  - HPA prognostic flag   →  { significant (p≤0.05), non-significant }
  - Tumour grade          →  { Grade 1 (low), Grade 2 (intermediate), Grade 3 (high) }
  - FIGO stage            →  { early (I–II), advanced (III–IV) }

Linguistic Output:
  - log₂ fold-change magnitude  →  { low (~±1), moderate (~±1.5), high (~±2) }

Defuzzification:
  - Output is a Gaussian-sampled continuous value, not a crisp category
  - This captures the graded, continuous nature of biological expression variation
```

---

## Design Rationale

| Design Choice | Reason |
|---|---|
| Gaussian sampling over fixed values | Models biological noise and individual variation |
| Larger σ for non-significant proteins (0.5) | Higher uncertainty when no prognostic evidence |
| Smaller σ for significant proteins (0.4) | Tighter, evidence-backed perturbation |
| log₂FC clipped to [−3, +3] | Matches the range reported in published endometrial cancer proteomic studies (Chen et al., 2021) |
| FIGO weights calibrated to AIHW 2023 | Ensures stage distribution matches real Australian population |
| HPA baseline as anchor | Grounds simulated profiles in validated tissue-level expression |

---

## Outputs

| File | Contents |
|---|---|
| `clinical_data.csv` | One row per patient; columns: Age, BMI, Menopause, Nulliparity, TumorType, Grade, FIGO_Stage, Treatment, Survival |
| `protein_abundance.csv` | One row per patient; columns = gene symbols; values = log₂ fold-change relative to normal endometrial tissue |

Both files are directly compatible with:
- Python: `pandas`, `scikit-learn`
- R: `limma`, `DESeq2`, Bioconductor pipelines

---

## References

- Zadeh, L.A. (1965). Fuzzy sets. *Information and Control*, 8(3), 338–353.
- Uhlen, M. et al. (2015). Tissue-based map of the human proteome. *Science*, 347, 1260419.
- Chen, W. et al. (2021). Privacy-preserving synthetic health data. *Artificial Intelligence in Medicine*, 117, 102082.
- Australian Institute of Health and Welfare (2023). *Cancer in Australia 2023*. AIHW, Canberra.
