\# BRCA scRNA-seq Analysis — Project Summary

\## MSc Dissertation 2026 — Anna McStay

\## Supervisor: Adrien Kissenpfennig



\---



\## Project Overview

Integrative analysis of single-cell RNA-seq data to map immune cell heterogeneity in 

breast cancer tumour microenvironments. Two publicly available datasets are used:



\- \*\*GSE114725\*\* — immune-enriched breast cancer dataset, 8 patients, 4 tissue types 

&#x20; (TUMOR, BLOOD, LYMPHNODE, NORMAL)

\- \*\*GSE176078\*\* — full TME dataset, 26 samples, Wu et al. 2021 (Nat Genet), 

&#x20; breast cancer subtypes ER+, HER2+, TNBC



\---



\## Repository

GitHub: https://github.com/annamcstay/brca-scrnaseq-analysis

Language: Python (scrna conda environment, Python 3.10)

Key tools: Scanpy, Dask, Zarr, Harmony, Scrublet, CellTypist, LIANA, gseapy, pydeseq2



\---



\## Data Files (local only, not on GitHub)

All in `C:\\Users\\annam\\Dissertation 2026\\Data\\`



\### Raw

\- `GSE114725\_raw.h5ad` — 47,016 cells × 14,875 genes, raw counts

\- `GSE176078\_raw.h5ad` — 100,064 cells × 29,733 genes, raw counts



\### Processed

\- `GSE114725\_phase1\_v2.h5ad` — 44,662 cells × 2,000 HVGs, post-QC + Harmony

\- `GSE176078\_phase1\_v2.h5ad` — 91,425 cells × 2,000 HVGs, post-QC + Harmony

\- `GSE114725\_phase2\_v2\_clustered.h5ad` — + Leiden clusters all resolutions

\- `GSE176078\_phase2\_v2\_clustered.h5ad` — + Leiden clusters all resolutions

\- `GSE114725\_phase2\_v2\_annotated.h5ad` — + cell\_type labels + CellTypist

\- `GSE176078\_phase2\_v2\_annotated.h5ad` — + cell\_type labels + CellTypist

\- `GSE114725\_tcells\_subclustered.h5ad` — T cell sub-clustering object



\---



\## Notebooks

\- `00\_environment\_setup\_checks.ipynb`

\- `01\_GSE114725\_data\_loading.ipynb`

\- `02\_GSE176078\_data\_loading.ipynb`

\- `03\_phase1\_QC\_V2.ipynb` — \*\*current v2 QC pipeline\*\*

\- `04\_phase2\_clustering\_annotation.ipynb` — \*\*current Phase 2\*\*

\- `05\_phase3\_DE\_pathways\_interactions.ipynb` — \*\*current Phase 3\*\*



\---



\## Phase 1 — QC Pipeline (v2, COMPLETE)

Notebook: `03\_phase1\_QC\_V2.ipynb`



\### Key steps in order:

1\. Load raw h5ad as Dask float32 sparse arrays (chunk\_size=2000)

2\. `gc.collect()` after each heavy step

3\. Filter V(D)J variable + constant region genes:

&#x20;  - Prefixes: IGHV, IGLV, IGKV, TRAV, TRBV, TRGV, TRDV, IGHD, IGHJ, IGLJ, IGKJ

&#x20;  - Constant: TRAC, TRBC, TRGC, TRDC, IGHA, IGHE, IGHG, IGHM, IGLC, IGKC

&#x20;  - Removed 42 genes from GSE114725, 372 from GSE176078

4\. Remove mitochondrial genes (MT- prefix) completely

&#x20;  - Removed 12 genes from GSE114725, 13 from GSE176078

5\. Chunked QC metrics (n\_genes\_by\_counts, total\_counts) — no MT% since MT genes removed

6\. QC filtering:

&#x20;  - GSE114725: n\_genes > 200 (mean=599, no MT threshold)

&#x20;  - GSE176078: n\_genes > 500 (mean=1763, no MT threshold)

7\. Gene filtering: min\_cells=3

8\. Per-sample Scrublet doublet detection (NOT full dataset — memory constraint)

&#x20;  - GSE114725: 8 samples, 649 doublets removed → 44,662 cells

&#x20;  - GSE176078: 26 samples, 197 doublets removed → 91,425 cells

9\. Normalise to 10,000 counts + log1p

&#x20;  - Check for double normalisation before running

10\. `adata.raw = adata` — stores full gene log-normalised matrix

11\. HVG selection: top 2,000, flavor="seurat"

&#x20;   - Required for memory — Scrublet and Harmony crash on full gene matrix

&#x20;   - V(D)J genes filtered before HVG selection in v2 (fixes original contamination issue)

12\. Scale (max\_value=10) — note: densifies sparse matrix

13\. PCA (svd\_solver="arpack") — adata1 uses standard; adata2 uses dense conversion trick

14\. Harmony batch correction:

&#x20;   - GSE114725: key="patient"

&#x20;   - GSE176078: key="orig.ident"

15\. UMAP before and after Harmony (random\_state=42)

16\. Save as compressed h5ad



\### Key design decisions:

\- \*\*Why HVGs if Adrien wanted full gene?\*\* Full gene PCA crashes laptop RAM — 

&#x20; PCA requires dense matrix internally. HVGs for PCA/Harmony/clustering only. 

&#x20; Full genes preserved in adata.raw for biological analyses.

\- \*\*Why per-sample Scrublet?\*\* Doublets only form within a sample. Full dataset 

&#x20; Scrublet caused MemoryError on both datasets.

\- \*\*Why MT genes removed rather than filtered by %?\*\* Adrien's recommendation. 

&#x20; Cleaner than threshold filtering.

\- \*\*Why V(D)J constant region genes added?\*\* Adrien's recommendation in GitHub 

&#x20; Issue #9 comment.



\---



\## Phase 2 — Clustering and Annotation (COMPLETE + all stretch tasks)

Notebook: `04\_phase2\_clustering\_annotation.ipynb`



\### Clustering:

\- Neighbours built on X\_pca\_harmony (critical — not X\_pca)

\- Leiden at resolutions 0.2, 0.4, 0.6, 0.8, 1.0

\- \*\*Selected resolution: 0.8 for GSE114725 (6 clusters), 0.6 for GSE176078 (26 clusters)\*\*

\- Justification: visual inspection of UMAP plots across resolutions



\### GSE114725 cluster labels (leiden\_0.8):"0": "T cells (resting)"

"1": "T cells (naive/memory)"

"2": "NK/Cytotoxic T cells"

"3": "Activated T cells"

"4": "Macrophages"

"5": "Monocytes/DC"

\### GSE176078 cluster labels (leiden\_0.6):

"0": "Endothelial cells"

"1": "Endothelial cells"

"2": "CAFs"

"3": "PVL"

"4": "Basal epithelial"

"5": "B cells"

"6": "Cycling cells"

"7": "Plasma cells"

"8": "Cycling epithelial"

"9": "CD8 T cells"

"10": "NK cells"

"11": "T cells"

"12": "Naive/memory T cells"

"13": "Luminal epithelial"

"14": "Macrophages"

"15": "Monocytes/DC"

"16": "Cycling myeloid"

"17": "pDC"

"18": "Luminal epithelial"

"19": "Luminal epithelial"

"20": "Epithelial"

"21": "Epithelial"

"22": "Luminal epithelial"

"23": "Luminal epithelial"

"24": "Luminal epithelial"

"25": "Luminal epithelial"



\### CellTypist annotation:

\- Model: Immune\_All\_High.pkl

\- GSE114725: run on adata.raw.to\_adata(), majority\_voting=True

\- GSE176078: batched (10,000 cells/batch), majority\_voting=False (memory constraint)

\- Validated against published Wu et al. 2021 labels — >95% concordance on all 

&#x20; major immune populations

\- Result stored in adata.obs\["celltypist"] and adata.obs\["majority\_voting"]



\### Marker genes:

\- Wilcoxon rank-sum tests

\- V(D)J prefixes filtered from marker interpretation

\- Saved as CSVs in results/phase2\_clustering\_v2/



\### Stretch tasks completed:

1\. \*\*T cell sub-clustering\*\* — 31,892 T/NK cells from GSE114725, leiden 0.7 → 5 subtypes:

&#x20;  NK/Cytotoxic, Resting T cells, Naive/Memory T cells, Activated T cells, 

&#x20;  B cell contamination

2\. \*\*Cell type proportions\*\* — per patient (GSE114725) and per subtype (GSE176078)

&#x20;  Key finding: TNBC highest macrophages (14.85%), HER2+ highest CD8 T cells (15.81%)

3\. \*\*Reference atlas comparison\*\* — CellTypist vs Wu et al. 2021 published labels

&#x20;  >94% concordance for all major populations



\---



\## Phase 3 — DE, Pathways, Cell-Cell Communication (PARTIAL)

Notebook: `05\_phase3\_DE\_pathways\_interactions.ipynb`



\### LIANA (COMPLETE):

\- Tool: LIANA v1.7.3, rank\_aggregate method

\- Resource: consensus (combination of multiple LR databases)

\- Parameters: expr\_prop=0.1, min\_cells=10

\- GSE114725: 458 significant interactions (specificity\_rank < 0.05)

&#x20; Key: S100A8/A9→CD36/CD68, CXCL8→CXCR2, IL1B→IL1R2

\- GSE176078: 3,086 significant interactions (immune cells only — memory constraint)

&#x20; Key: CCL19→CCR7 (Monocytes/DC → T/B/NK cells), CXCL9→FCGR2A, CD14→TLR9

\- Results saved in results/phase3/



\### DE analysis (PENDING — awaiting Adrien's reply):

\- Question on GitHub Issue #11: Scanpy Wilcoxon vs PyDESeq2

\- PyDESeq2 needs raw integer counts — adata.raw has log-normalised not raw counts

\- Proposed: Scanpy Wilcoxon on v2 objects directly

\- Planned comparisons:

&#x20; - GSE114725: TUMOR vs BLOOD within T cells, Macrophages, NK cells

&#x20; - GSE176078: ER+ vs HER2+ vs TNBC within immune cell types



\### Pathway enrichment (PENDING — depends on DE):

\- Tool: gseapy, gene\_sets=\["MSigDB\_Hallmark\_2020", "KEGG\_2021\_Human"]

\- Background: all genes detected in relevant cell type



\---



\## Key Technical Issues and Resolutions



| Issue | Resolution |

|-------|-----------|

| harmonypy install fails on Python 3.13 Windows | Use scrna env (Python 3.10) |

| Scrublet MemoryError on full dataset | Run per-sample instead |

| CellTypist ValueError on scaled data | Use adata.raw.to\_adata() |

| Full gene PCA MemoryError | HVG selection before PCA; adata.raw preserves full genes |

| LIANA MemoryError on GSE176078 | Subset to immune cells only |

| CellTypist MemoryError on GSE176078 | Batch size 10,000 cells |

| Git detached HEAD / merge conflicts | Always: git pull --rebase before push |

| Kernel restarts losing variables | Save h5ad after every major step |



\---



\## GitHub Issues Status

\- #9 memory issues — CLOSED (resolved in v2)

\- #10 RAM issues — CLOSED (resolved in v2)

\- #11 Phase 3 DE — OPEN (awaiting Adrien reply on Wilcoxon vs PyDESeq2)

\- #12 Phase 3 DE GSE176078 subtypes — OPEN

\- #13 Phase 3 LIANA — CLOSED (complete)



\---



\## Important Paths

```python

PROJECT\_DIR = Path(r"C:\\Users\\annam\\Dissertation 2026")

RAW\_DIR = PROJECT\_DIR / "Data" / "Raw"

PROCESSED\_DIR = PROJECT\_DIR / "Data" / "Processed"

FIGURE\_DIR\_P2 = PROJECT\_DIR / "figures" / "phase2\_clustering\_v2"

FIGURE\_DIR\_P3 = PROJECT\_DIR / "figures" / "phase3"

RESULTS\_DIR\_P2 = PROJECT\_DIR / "results" / "phase2\_clustering\_v2"

RESULTS\_DIR\_P3 = PROJECT\_DIR / "results" / "phase3"

```



\---



\## Standard Setup Cell (use at start of every session)

```python

import scanpy as sc

import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

import liana

import gseapy as gp

import gc

from pathlib import Path

from scipy.sparse import issparse



PROJECT\_DIR = Path(r"C:\\Users\\annam\\Dissertation 2026")

PROCESSED\_DIR = PROJECT\_DIR / "Data" / "Processed"

FIGURE\_DIR = PROJECT\_DIR / "figures" / "phase3"

RESULTS\_DIR = PROJECT\_DIR / "results" / "phase3"



cluster\_labels\_1 = {

&#x20;   "0": "T cells (resting)", "1": "T cells (naive/memory)",

&#x20;   "2": "NK/Cytotoxic T cells", "3": "Activated T cells",

&#x20;   "4": "Macrophages", "5": "Monocytes/DC"

}



cluster\_labels\_2 = {

&#x20;   "0": "Endothelial cells", "1": "Endothelial cells",

&#x20;   "2": "CAFs", "3": "PVL", "4": "Basal epithelial",

&#x20;   "5": "B cells", "6": "Cycling cells", "7": "Plasma cells",

&#x20;   "8": "Cycling epithelial", "9": "CD8 T cells",

&#x20;   "10": "NK cells", "11": "T cells", "12": "Naive/memory T cells",

&#x20;   "13": "Luminal epithelial", "14": "Macrophages",

&#x20;   "15": "Monocytes/DC", "16": "Cycling myeloid", "17": "pDC",

&#x20;   "18": "Luminal epithelial", "19": "Luminal epithelial",

&#x20;   "20": "Epithelial", "21": "Epithelial",

&#x20;   "22": "Luminal epithelial", "23": "Luminal epithelial",

&#x20;   "24": "Luminal epithelial", "25": "Luminal epithelial"

}



IMMUNE\_TYPES = \[

&#x20;   "T cells", "CD8 T cells", "NK cells", "Naive/memory T cells",

&#x20;   "B cells", "Plasma cells", "Macrophages", "Monocytes/DC", "pDC"

]

```



\---



\## Timeline Status (Week 6 of 18)

\- Phase 0 ✓ Complete

\- Phase 1 ✓ Complete (v2 pipeline)

\- Phase 2 ✓ Complete including all 3 stretch tasks

\- Phase 3 ⚡ LIANA complete, DE pending Adrien reply

\- Phase 4 ✗ Not started (Weeks 9-11 target)

\- Phase 5 ✗ Not started (Weeks 12-18)



Grade trajectory: Merit territory, distinction possible with DE + interpretation

