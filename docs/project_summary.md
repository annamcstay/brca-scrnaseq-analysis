# BRCA scRNA-seq MSc Dissertation — Project Summary
**Student:** Anna McStay  
**Supervisor:** Adrien Kissenpfennig (GitHub: AKissenpfennig)  
**Repository:** https://github.com/annamcstay/brca-scrnaseq-analysis  
**Environment:** scrna conda env (Python 3.10), Windows laptop  
**Last updated:** Week 9 of 18

---

## Project Overview
MSc dissertation investigating immune cell heterogeneity in the breast cancer tumour microenvironment (TME) using single-cell RNA sequencing. Two complementary publicly available datasets were selected to answer different but complementary biological questions.

---

## Datasets

### GSE114725 (Azizi et al. 2018, Cell)
- **Title:** Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment
- **Cells:** 47,016 raw → 44,662 post-QC (95% retained)
- **Patients:** 8 breast cancer patients
- **Tissues:** Tumour, Blood, Normal, Lymph node (within-patient comparisons)
- **Protocol:** InDrop v2, CD45+ immune-enriched
- **Key feature:** Within-patient tissue comparisons — ideal for tumour vs blood immune cell analysis
- **Limitation:** Small cohort, no clinical subtype information

### GSE176078 (Wu et al. 2021, Nature Genetics)
- **Title:** A single-cell and spatially resolved atlas of human breast cancers
- **Cells:** 100,064 raw → 91,425 post-QC (91% retained)
- **Patients:** 26 patients across ER+, HER2+, TNBC subtypes
- **Protocol:** 10x Chromium, full TME (immune + epithelial + stromal)
- **Key feature:** Three clinical subtypes, published cell type annotations for validation
- **Limitation:** No within-patient tissue comparisons, between-patient subtype comparison

**Complementarity:** GSE114725 asks "how does the tumour alter immune cell behaviour?" (tumour vs blood, same patients). GSE176078 asks "how does breast cancer subtype shape the immune landscape?" (ER+ vs HER2+ vs TNBC, different patients).

---

## File Paths
PROJECT_DIR = Path(r"C:\Users\annam\Dissertation 2026")
RAW_DIR = PROJECT_DIR / "Data" / "Raw"
PROCESSED_DIR = PROJECT_DIR / "Data" / "Processed"
FIGURE_DIR_QC = PROJECT_DIR / "figures" / "phase1_qc_v2"
FIGURE_DIR_P2 = PROJECT_DIR / "figures" / "phase2_clustering_v2"
FIGURE_DIR_P3 = PROJECT_DIR / "figures" / "phase3"
RESULTS_DIR_P2 = PROJECT_DIR / "results" / "phase2_clustering_v2"
RESULTS_DIR_P3 = PROJECT_DIR / "results" / "phase3"
RESULTS_DIR_CORRECTED = PROJECT_DIR / "results" / "phase3" / "corrected"
RESULTS_DIR_LIANA = PROJECT_DIR / "results" / "phase3" / "corrected" / "liana"

---

## Key Processed Files

### Raw (true integer counts — use for PyDESeq2)
- `Data/Raw/GSE114725_raw.h5ad` — 332MB, 47,016 cells × 14,875 genes, integer counts (max ~1236)
- `Data/Raw/GSE176078_raw.h5ad` — 1374MB, 100,064 cells × 27,343 genes, integer counts (max ~13141)

### Phase 1 — QC processed (2000 HVGs, scaled, Harmony)
- `Data/Processed/GSE114725_phase1_v2.h5ad` — log-normalised in adata.raw, scaled in adata.X
- `Data/Processed/GSE176078_phase1_v2.h5ad` — same

### Phase 2 — Annotated (corrected cluster labels)
- `Data/Processed/GSE114725_phase2_v2_annotated.h5ad` — leiden_0.6, corrected cell_type labels
- `Data/Processed/GSE176078_phase2_v2_annotated_corrected.h5ad` — leiden_0.6, corrected labels

### Phase 2 — Sub-clustered objects
- `Data/Processed/GSE114725_tcells_subclustered.h5ad` — 31,718 T cells, 5 subtypes including mast cells
- `Data/Processed/GSE114725_macrophages_subclustered.h5ad` — 8,624 macrophages, 4 subtypes

### IMPORTANT — Misnamed file (DO NOT USE FOR PyDESeq2)
- `Data/Processed/GSE114725_phase1_v2_rawcounts.h5ad` — MISLEADING NAME — contains log-normalised data not raw counts. Safe for LIANA but NOT for PyDESeq2. Use GSE114725_raw.h5ad instead.

---

## Notebooks (in order)

| Notebook | Status | Contents |
|---|---|---|
| 00_environment_setup_checks.ipynb | Complete | Environment verification |
| 01_GSE114725_data_loading.ipynb | Complete | Data download and initial loading |
| 02_GSE176078_data_loading.ipynb | Complete | Data download and initial loading |
| 03_phase1_QC_V2.ipynb | Complete | v2 QC pipeline both datasets |
| 04_phase2_clustering_annotation_clean.ipynb | Complete | **CLEAN** clustering notebook, runs top to bottom |
| 05_phase3_DE_pathways_interactions.ipynb | SUPERSEDED | Old exploratory DE+LIANA, use 07 and 08 instead |
| 06_phase3_DE_pseudobulk.ipynb | SUPERSEDED | DE on log-normalised data — statistically invalid |
| 07_phase3_DE_pseudobulk_corrected.ipynb | **CURRENT** | Corrected DE with true raw counts + pathway enrichment |
| 08_phase3_LIANA_clean.ipynb | **CURRENT** | Corrected LIANA with updated cluster labels |

---

## Phase 1 — QC Pipeline (COMPLETE)

### Key steps:
1. V(D)J variable + constant region gene removal (prevents clonal diversity skewing clustering)
2. Mitochondrial gene removal (entirely, not filtered by percentage)
3. Per-sample Scrublet doublet detection (649 removed GSE114725, 197 GSE176078)
4. Cell filtering: n_genes > 200 (GSE114725), n_genes > 500 (GSE176078)
5. Normalise (sc.pp.normalize_total, target_sum=1e4) + log1p
6. **adata.raw = adata** ← snapshot taken here (log-normalised, full gene set, used for LIANA/CellTypist/markers)
7. HVG selection (top 2,000)
8. Scaling (sc.pp.scale, max_value=10) ← adata.X is scaled from here onward
9. PCA (50 components)
10. Harmony batch correction (by patient)

### Data layers explained:
- `adata.X` — scaled, HVG-only — use for PCA/clustering ONLY
- `adata.raw.X` — log-normalised, full gene set, unscaled — use for LIANA, CellTypist, marker genes
- `GSE114725_raw.h5ad` — true integer counts, all genes, all cells pre-QC — use for PyDESeq2 (subset to QC barcodes first)

---

## Phase 2 — Clustering and Annotation (COMPLETE)

### GSE114725 — Final cluster labels (leiden_0.6, 6 clusters)
```python
cluster_labels_1 = {
    "0": "T cells",              # CD3D/CD3E high, ribosomal dominant
    "1": "Activated T cells",    # FOS/JUN/DUSP1/NFKBIA — immediate early genes
    "2": "NK/Cytotoxic T cells", # NKG7/GNLY/PRF1/KLRD1/GZMB
    "3": "Macrophages",          # CD68/LYZ/TYROBP/AIF1 — confirmed by canonical markers
    "4": "B cells",              # CD79A/MS4A1/CD19 — n=790, genuine B cell cluster
    "5": "Monocytes/DC"          # S100A8/S100A9/CD14/FCGR3A/VCAN
}
```

**Resolution history:** Tested 0.2, 0.4, 0.6, 0.8, 1.0. Resolution 0.8 previously used but had ambiguous T cell split (resting vs naive/memory not distinguishable by canonical markers). Resolution 0.6 adopted — merges ambiguous T cell clusters into one, reveals genuine B cell population.

### GSE176078 — Final cluster labels (leiden_0.6, 26 clusters)
```python
cluster_labels_2 = {
    "0": "Endothelial cells", "1": "Endothelial cells",
    "2": "CAFs", "3": "PVL", "4": "Basal epithelial",
    "5": "B cells", "6": "Cycling cells", "7": "Plasma cells",
    "8": "Cycling epithelial", "9": "CD8 T cells",
    "10": "NK cells", "11": "T cells", "12": "Naive/memory T cells",
    "13": "Luminal epithelial", "14": "Macrophages",
    "15": "Unassigned",  # n=22, no clear lineage markers
    "16": "Unassigned",  # n=19, no clear lineage markers
    "17": "pDC",         # LILRA4/CLEC4C/IRF7 confirmed
    "18": "Luminal epithelial", "19": "Luminal epithelial",
    "20": "Epithelial (ambiguous)",  # non-specific metabolic signature
    "21": "Epithelial (ambiguous)",  # ribosomal dominant, ambiguous
    "22": "Luminal epithelial", "23": "Luminal epithelial",
    "24": "Luminal epithelial", "25": "Luminal epithelial"
}
```

**Resolution note:** Lower resolutions (0.2-0.4) tested but merge biologically distinct immune populations (e.g. CD8 T cells + NK cells + Naive/memory T cells at res 0.2) — resolution 0.6 retained.

### Validation against Wu et al. 2021
>95% concordance all major immune populations:
- NK cells: 100%, CD8 T cells: 99.9%, T cells: 99.8%, Naive/memory T cells: 99.6%
- Macrophages: 98.8%, Plasma cells: 96.1%, B cells: 96.9%, pDC: 99.7%
- CAFs: 89.3%, PVL: 94.6% (lower due to annotation granularity, not errors)
- Unassigned (C15/C16): confirmed mapping to mixed Cancer Epithelial/Myeloid — correctly relabelled

### CellTypist annotation
- Model: Immune_All_High.pkl
- Run on adata.raw.to_adata() (log-normalised, unscaled)
- GSE114725: majority_voting=True
- GSE176078: batched (10,000 cells) due to memory constraints

---

## Phase 2 Stretch Tasks

### Stretch 1 — T cell sub-clustering (GSE114725) ✓ COMPLETE
- **Input:** 31,718 cells from T cells + Activated T cells + NK/Cytotoxic T cells clusters
- **Method:** Leiden resolution 0.7 on Harmony PCA
- **Result:** 5 subtypes:
  - Resting T cells (n=10,057) — CD3D high, TCF7 high, ribosomal dominant
  - Naive/Memory T cells (n=7,307) — CD3D high, CCR7 slightly higher
  - Activated T cells (n=7,953) — FOS/JUN/DUSP1 high
  - NK/Cytotoxic T cells (n=6,200) — NKG7/GNLY/PRF1
  - **Mast cells (n=201)** — CPA3/KIT/TPSB2 — unexpected finding, captured at cluster boundary
- **CD4/CD8 separation:** Attempted but not achievable at mRNA level — CD4 and CD8A both expressed at low overlapping levels. Requires CITE-seq protein data (as used by Wu et al. 2021)
- **Saved:** `GSE114725_tcells_subclustered.h5ad`

### Stretch 2 — Macrophage sub-clustering (GSE114725) ✓ COMPLETE
- **Input:** 8,624 macrophage cells
- **Method:** Leiden resolution 0.7 on Harmony PCA
- **Result:** 4 subtypes:
  - Resting/Resident (n=2,020) — ribosomal dominant, low activation
  - LAM-like/Antigen-presenting (n=2,964) — C1QA/APOE/HLA-DRA — LAM-like and antigen-presenting not cleanly separable at this resolution, merged. Consistent with macrophage plasticity continuum (Azizi et al. 2018)
  - **Cytotoxic macrophages NKG7+ (n=1,953)** — NKG7/PRF1/GNLY — tumour-enriched (1,014 tumour vs 688 blood) — real population, not contamination
  - Monocyte-like S100A8/A9+ (n=1,687) — recently recruited
- **Note:** Previous 6-cluster solution (from earlier sessions) could not be reproduced after corrected broad clustering — 4-cluster solution is the correct current version
- **Saved:** `GSE114725_macrophages_subclustered.h5ad`

### Stretch 3 — Cell type proportions ✓ PARTIALLY COMPLETE
- GSE114725: proportions per patient calculated and visualised ✓
- GSE176078: proportions per subtype (ER+/HER2+/TNBC) calculated and visualised ✓
- **Formal statistical tests: NOT YET RUN** — scCODA had Windows dependency conflicts, propeller requires R
- **Action needed:** Implement scipy permutation test or confirm tool with Adrien

### Stretch 4 — Reference atlas comparison ✓ COMPLETE
- GSE176078 only (GSE114725 has no published per-cell annotations)
- >95% concordance across all major immune populations
- Saved: `GSE176078_our_vs_published_comparison.csv`

---

## Phase 3 — Differential Expression (CORRECTED — COMPLETE)

### CRITICAL CORRECTION
All previous DE results (notebook 06) were run on log-normalised data. PyDESeq2 requires true raw integer counts. All results have been re-run using:
- **Option 1 (Adrien's recommendation):** Subset `GSE114725_raw.h5ad` / `GSE176078_raw.h5ad` to QC-surviving barcodes
- Results saved in `results/phase3/corrected/`
- Old results moved to `results/phase3/superseded/`

### GSE114725 Broad Cluster DE (TUMOR vs BLOOD) ✓ CORRECTED
| Cell type | Significant DEGs |
|---|---|
| T cells | 692 |
| NK/Cytotoxic T cells | 333 |
| Activated T cells | 148 |
| Macrophages | 408 |
| Monocytes/DC | 253 |

**Key findings (confirmed with corrected data):**
- TIPIN↓ across ALL cell types — reduced DNA replication/proliferative capacity in TME
- ATF3↑ across ALL cell types — universal stress response
- HLA-DRA↑ — increased antigen presentation in tumour
- SPP1↑ in Activated T cells — immune evasion signal
- CCL22↑ in Activated T cells — Treg-recruiting chemokine

### GSE114725 T Cell Subcluster DE (TUMOR vs BLOOD) ✓ CORRECTED
| T cell subtype | Significant DEGs |
|---|---|
| Resting T cells | 533 |
| Naive/Memory T cells | 330 |
| Activated T cells | 324 |
| NK/Cytotoxic T cells | 250 |
| Mast cells | excluded (n=201, too small) |

**Key finding:** TIPIN↓ and NOSIP↓ confirmed across ALL 4 T cell subtypes — exhaustion/reduced proliferation signal is universal across T cell compartment, not subtype-specific. IFNG↑ specifically in Resting T cells.

### GSE114725 Macrophage Subcluster DE (TUMOR vs BLOOD) ✓ CORRECTED
| Macrophage subtype | Significant DEGs |
|---|---|
| Resting/Resident | 228 |
| LAM-like/Antigen-presenting | 172 |
| Cytotoxic macrophages (NKG7+) | 245 |
| Monocyte-like (S100A8/A9+) | 212 |

**Key findings:**
- LAM-like: FABP5↑ (lipid metabolism, direct LAM validation), GPNMB↑
- Cytotoxic macrophages: GPNMB↑ (7.09 logFC), ACP5↑ — overlap with LAM signature, consistent with plasticity continuum
- All subtypes: ATF3↑, TIPIN↓ confirmed

### GSE176078 Subtype DE ✓ CORRECTED (NULL RESULT — BIOLOGICAL FINDING)
All three subtype comparisons (TNBC vs ER+, TNBC vs HER2+, HER2+ vs ER+) across 6 immune cell types yielded 0-16 significant DEGs. Even with improved power (ER+ vs TNBC, n=11 vs n=10) results remained minimal.

**Conclusion:** Immune cell transcriptional states are conserved across breast cancer subtypes — subtype differences are **compositional not transcriptional**. Consistent with Wu et al. 2021. Pathway enrichment not run (insufficient DEGs for valid enrichment).

---

## Phase 3 — Pathway Enrichment (CORRECTED — COMPLETE)

All pathway enrichment run on corrected DE results. Method: gseapy, MSigDB Hallmark + KEGG 2021, background = genes detected per cell type.

### GSE114725 Broad Cluster — Top pathways (tumour upregulated)
- **TNF-alpha/NF-kB** — #1 hit in ALL 5 cell types (padj 1e-10 to 1e-63)
- Inflammatory Response
- Interferon Gamma Response
- IL-2/STAT5 Signalling
- IL-6/JAK/STAT3 (Macrophages only)
- **Cholesterol Homeostasis (Monocytes/DC)** — new finding with corrected data, LAM-like lipid signal

**Consistently downregulated:**
- Ribosome (KEGG) — reduced protein synthesis all cell types
- Myc Targets V1 — reduced proliferative programme

### T cell subcluster pathway enrichment
- TNF-alpha/NF-kB #1 in ALL 4 T cell subtypes
- NK/Cytotoxic uniquely: Hypoxia, Apoptosis
- Ribosome/Myc consistently downregulated

### Macrophage subcluster pathway enrichment
- **LAM-like: Cholesterol Homeostasis** — directly validates LAM macrophage identity
- Resting/Resident: EMT upregulated — tissue remodelling
- Cytotoxic: Complement upregulated

---

## Phase 3 — LIANA Cell-Cell Communication (CORRECTED — COMPLETE)

### CRITICAL CORRECTION AND RETRACTION
Previous GSE176078 headline finding (CCL19→CCR7) was retracted. It was entirely driven by a 22-cell "Monocytes/DC" cluster (C15) that failed canonical marker validation and was relabelled "Unassigned". This cluster had artificially high specificity in LIANA due to its small size and unique expression profile — not biologically meaningful.

### GSE114725 LIANA (corrected, resolution 0.6 labels)
- **Total interactions tested:** 2,356
- **Significant interactions (specificity_rank < 0.05):** 184
- **Top interactions:**
  - Monocytes/DC → B cells: CXCL8-CD79A (alarmin to B cell receptor)
  - Monocytes/DC autocrine: C1QB-LRP1, SERPINA1-LRP1, S100A8-CD68
  - Monocytes/DC → Monocytes/DC: complement signalling dominant
- **Note:** Down from 458 with old labels — reduction due to T cell cluster merging (fewer source/target combinations) and changes in relative specificity across network

### GSE176078 LIANA (corrected, immune subset 43,798 cells, no Monocytes/DC)
- **Total interactions tested:** 5,720
- **Significant interactions (specificity_rank < 0.05):** 413
- **NEW headline finding: APOE-TREM2 macrophage autocrine signalling** — directly validates LAM macrophage identity from sub-clustering
- **Top interactions:**
  - Macrophages → Plasma cells: CXCL8-SDC1, FN1-SDC1
  - Macrophages autocrine: C1QB-LRP1, S100A9-CD68, APOE-TREM2, APOE-ABCA1, PLAU-PLAUR
  - pDC → Macrophages: SERPINF1-PLXDC2
- **Saved in:** `results/phase3/corrected/liana/`

---

## Results Folder Structure
results/
├── phase2_clustering_v2/          # All Phase 2 results (markers, proportions, annotations)
└── phase3/
├── corrected/                 # ALL VALID CURRENT RESULTS
│   ├── DE/                    # Corrected DE CSVs (GSE114725 + GSE176078)
│   ├── pathways/              # Corrected pathway enrichment CSVs
│   └── liana/                 # Corrected LIANA results (both datasets)
└── superseded/                # OLD INVALID RESULTS (log-normalised DE input)

---

## Figures Folder Structure
figures/
├── phase1_qc_v2/                  # QC figures (violin plots, Harmony before/after, scrublet)
├── phase2_clustering_v2/          # UMAP figures, proportions, sub-clustering
└── phase3/
├── [corrected figures]        # All _corrected.png files — USE THESE
└── superseded/                # Old figures based on wrong DE/clustering

### Current valid figures (phase3/):
- `GSE114725_volcano_*_corrected.png` — 5 volcano plots (T cells, NK/Cytotoxic, Activated, Macrophages, Monocytes/DC)
- `GSE114725_tcell_subcluster_DE_heatmap_corrected.png`
- `GSE114725_macrophage_subcluster_DE_heatmap_corrected.png`
- `GSE114725_liana_heatmap_corrected.png`
- `GSE176078_liana_heatmap_corrected.png`
- `GSE114725_liana_dotplot_corrected.png`
- `GSE176078_liana_dotplot_corrected.png`

---

## Key Biological Findings (Corrected, Defensible)

### GSE114725 (Azizi et al.) — Tumour vs Blood
1. **TNF-alpha/NF-kB** is the dominant inflammatory pathway across ALL immune cell types in tumour — universal, not cell-type specific
2. **TIPIN downregulation** across ALL immune cell types and ALL T cell subtypes — reduced proliferative/DNA replication capacity is a universal feature of tumour-infiltrating immune cells
3. **TCF7 downregulation** in T cells — loss of stem-like/naive marker, early exhaustion signal
4. **Mast cells** identified at T cell/NK cluster boundary — unexpected finding from sub-clustering correction
5. **Cytotoxic macrophages (NKG7+)** identified as tumour-enriched macrophage subtype — new finding from corrected macrophage sub-clustering
6. **LAM-like macrophages** validated by FABP5↑, Cholesterol Homeostasis pathway, APOE-TREM2 signalling — cross-validates Wu et al. 2021 LAM1/LAM2

### GSE176078 (Wu et al.) — Subtype comparisons
1. **Compositional not transcriptional:** Immune cell transcriptional states conserved across ER+/HER2+/TNBC — subtype differences are in which cells are recruited, not how they behave
2. **TNBC:** Macrophage-rich (14.9%), plasma cell-rich (5.2%) — myeloid-dominant landscape
3. **HER2+:** T cell-rich (CD8 15.8%, naive/memory 24.9%) — most immune-infiltrated
4. **ER+:** Epithelial-dominant (39%), low immune infiltration
5. **APOE-TREM2** macrophage autocrine signalling — new LIANA headline finding, validates LAM macrophage identity

### Retracted findings
- **CCL19→CCR7** (GSE176078) — retracted. Was artefact of 22-cell "Monocytes/DC" cluster that failed canonical marker validation. Cluster relabelled Unassigned and removed from analysis.

---

## Corrections Made This Session

### Correction 1 — GSE114725 clustering resolution
- **Wrong:** Resolution 0.8 had ambiguous T cell resting/naive split (not supported by canonical markers)
- **Fixed:** Resolution 0.6 — T cells merged into one valid cluster, genuine B cells (n=790) identified

### Correction 2 — GSE176078 cluster annotations
- **Wrong:** C15 "Monocytes/DC" (n=22), C16 "Cycling myeloid" (n=19) — too small, non-specific markers. C20/C21 "Epithelial" — non-specific metabolic signatures
- **Fixed:** C15/C16 → Unassigned; C20/C21 → Epithelial (ambiguous)

### Correction 3 — T cell sub-cluster labels
- **Wrong:** Cluster 3 labelled "B cell contamination"
- **Fixed:** Relabelled "Mast cells" (CPA3=4.43, KIT=3.24, TPSB2=5.46) — genuine finding

### Correction 4 — Macrophage sub-clustering
- **Wrong:** 6-cluster solution from earlier session, including "NK/T cell contamination" and "B cell contamination" labels
- **Fixed:** 4-cluster solution — Resting/Resident, LAM-like/Antigen-presenting, Cytotoxic macrophages (NKG7+), Monocyte-like

### Correction 5 — PyDESeq2 input data (MOST CRITICAL)
- **Wrong:** All DE run on log-normalised data (max ~8, non-integer). PyDESeq2 negative binomial model requires integer counts
- **Fixed:** Used GSE114725_raw.h5ad / GSE176078_raw.h5ad (integer counts confirmed) subsetted to QC-surviving barcodes (Option 1)
- **Impact:** DEG counts changed but biological story consistent — all confirmed with valid data

### Correction 6 — LIANA cluster labels
- **Wrong:** LIANA run with old cluster labels including ambiguous T cell split and 22-cell Monocytes/DC in GSE176078
- **Fixed:** Re-run with corrected labels. CCL19-CCR7 finding retracted (artefact of 22-cell cluster). New finding: APOE-TREM2

---

## GitHub Issues Status
- #22 Fix raw counts — CLOSED (used original raw files via Option 1)
- #23 Re-cluster GSE114725 — CLOSED (res 0.6, 6 clusters confirmed)
- #24 Check GSE176078 clusters — CLOSED (corrected, Unassigned/ambiguous labels)
- #25 Re-run DE with correct raw counts — CLOSED (notebook 07, all DE corrected)
- #26 Re-run pathway enrichment — CLOSED (on corrected DE results)
- #27 Re-run LIANA if labels change — CLOSED (re-run, CCL19-CCR7 retracted)
- #28 PAGA trajectory inference — OPEN (Phase 3 stretch, next priority)

---

## Pending Tasks (Priority Order)

### Immediate
1. Composition statistical testing — scipy permutation test or propeller (R). Closes Phase 2 stretch gap. **Not yet started.**
2. Update docs: methods_clustering.md, methods_phase3.md
3. Final figure check in new chat (image limit reached in current chat)
4. Push all remaining changes

### Phase 3 Stretch
1. **PAGA trajectory inference** (Issue #28) — most feasible. T cell subclusters already saved. TCF7 exhaustion signal provides hypothesis (naive → activated → exhausted)
2. TCGA bulk deconvolution — lower priority
3. GSE176078 macrophage sub-clustering — possible but memory constrained

### Phase 4
- Sensitivity analyses (QC threshold changes, clustering resolution stability)
- Literature cross-check (Wu et al. 2021 exhaustion/LAM signatures)
- Lock final figure set (5-10 main + supplementary)
- Biological interpretation paragraphs for Results chapter

### Phase 5 (Weeks 12-18)
- Thesis writing: Methods → Results → Introduction → Discussion → Abstract

---

## Timeline
- **Week:** 9 of 18
- **Phases 0-3:** Complete (with corrections)
- **Phase 4:** Not started
- **Estimated remaining coding:** 6-10 hours (composition stats, PAGA, sensitivity analyses)
- **Status:** Ahead of schedule — Phase 3 due Week 11, complete at Week 9

---

## Standard Session Setup
```python
import scanpy as sc, pandas as pd, numpy as np, matplotlib.pyplot as plt
import liana as li, gseapy as gp, gc
from pathlib import Path
from scipy.sparse import issparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

PROJECT_DIR = Path(r"C:\Users\annam\Dissertation 2026")
PROCESSED_DIR = PROJECT_DIR / "Data" / "Processed"
RAW_DIR = PROJECT_DIR / "Data" / "Raw"
RESULTS_DIR = PROJECT_DIR / "results" / "phase3" / "corrected"
FIGURE_DIR = PROJECT_DIR / "figures" / "phase3"

cluster_labels_1 = {
    "0": "T cells", "1": "Activated T cells",
    "2": "NK/Cytotoxic T cells", "3": "Macrophages",
    "4": "B cells", "5": "Monocytes/DC"
}

cluster_labels_2 = {
    "0": "Endothelial cells", "1": "Endothelial cells",
    "2": "CAFs", "3": "PVL", "4": "Basal epithelial",
    "5": "B cells", "6": "Cycling cells", "7": "Plasma cells",
    "8": "Cycling epithelial", "9": "CD8 T cells",
    "10": "NK cells", "11": "T cells", "12": "Naive/memory T cells",
    "13": "Luminal epithelial", "14": "Macrophages",
    "15": "Unassigned", "16": "Unassigned", "17": "pDC",
    "18": "Luminal epithelial", "19": "Luminal epithelial",
    "20": "Epithelial (ambiguous)", "21": "Epithelial (ambiguous)",
    "22": "Luminal epithelial", "23": "Luminal epithelial",
    "24": "Luminal epithelial", "25": "Luminal epithelial"
}
```
