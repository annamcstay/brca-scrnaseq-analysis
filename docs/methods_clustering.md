# Methods: Clustering and Cell Type Annotation (Phase 2)

**Datasets:** GSE114725 (Azizi et al. 2018, CD45+ sorted, breast cancer immune cells, tumour/blood/normal/lymph node) and GSE176078 (Wu et al. 2021, whole-tumour breast cancer scRNA-seq, ER+/HER2+/TNBC).

This document records the clustering and annotation decisions made for Phase 2, and the reasoning behind each one. Written during active analysis rather than reconstructed afterward, per standard practice for keeping an accurate record of parameter choices.

---

## 1. Dimensionality reduction

PCA was run on the 2,000 highly variable genes selected in Phase 1 (`sc.tl.pca`, `svd_solver="arpack"`), producing 50 principal components. **`n_pcs=30`** was used as input to the neighbour graph and downstream clustering for both datasets.

**Justification:** A variance-ratio elbow plot showed the informative signal flattening at approximately PC 8–12 for both datasets, meaning `n_pcs=30` retains substantially more components than the minimum needed. This was investigated rather than assumed reasonable by default:

- A sensitivity check compared clustering under `n_pcs=30` vs `n_pcs=12` using Adjusted Rand Index (ARI). Results: GSE114725 ARI = 0.615, GSE176078 ARI = 0.725 — moderate agreement, not the "PC count barely matters" result initially expected.
- To interpret this properly, `n_pcs=12` clusters were cross-tabulated against the already-validated `cell_type` labels (from `n_pcs=30`). For GSE114725, all 9 clusters mapped to a single dominant cell type at 82.5–100% concordance. For GSE176078, 16 of 19 clusters showed >80% concordance; the remaining 3 involved known-ambiguous regions (a 28-cell artefact cluster, cycling epithelial cells, and closely related T-cell subtypes) — cases where boundary reassignment is expected regardless of PC count, not evidence of lost biological signal.

**Conclusion:** `n_pcs=30` was retained. The moderate ARI reflects boundary reshuffling in genuinely ambiguous regions rather than loss of real separation, and the resulting clusters have been independently validated against canonical markers, CellTypist, and (for GSE176078) published reference labels — a stronger justification than PC count sensitivity testing alone would provide.

*(Figure: `PCA_elbow_plots.png`)*

---

## 2. Reproducibility: random seeding and single-threaded execution

Two issues were identified and fixed during Phase 1/2 development that affect how results should be interpreted if this pipeline is ever re-run:

**Harmony batch correction and PCA were initially unseeded.** `harmonypy` initialises its correction using `sklearn.KMeans`, which is stochastic by default. Two otherwise-identical runs of the pipeline produced materially different clustering (e.g. GSE114725 res 0.6 gave 6 clusters in one run, 13 in another). `random_state=0` was added to `sc.tl.pca()` and `sce.pp.harmony_integrate()` to fix this.

**Even after seeding, GSE114725's resolution counts remained unstable** (still shifting between runs). This was traced to `sc.pp.neighbors()`'s approximate nearest-neighbour search (`pynndescent`), which runs multi-threaded and is not bit-reproducible across threads even with a fixed seed. `NUMBA_NUM_THREADS=1` and `OMP_NUM_THREADS=1` were set before any scanpy import to force single-threaded execution, and `random_state=0` was added explicitly to every `neighbors()` and `leiden()` call.

**Practical consequence:** the originally-reported 6-cluster solution for GSE114725 at resolution 0.6 was not a stable, reproducible result — it was one outcome of an unseeded stochastic process. Once seeding was fixed, resolution 0.6 consistently produced 13 clusters (reproducible, but with mixed-identity populations on marker inspection — see Section 3). This is why the working resolution for GSE114725 changed from 0.6 to 0.2 partway through analysis.

---

## 3. Clustering (Leiden) and resolution selection

Leiden clustering (`flavor="igraph"`, seeded, single-threaded) was run at five resolutions (0.2, 0.4, 0.6, 0.8, 1.0) for both datasets, on the Harmony-corrected PCA embedding.

**GSE114725: resolution 0.2 used (9 clusters).** Res 0.6 (13 clusters, the originally intended resolution) was re-examined once seeding was fixed, because a canonical marker screen showed clusters with mixed identity — e.g. one cluster elevated in both macrophage markers (CD68, LYZ) and monocyte markers (S100A8) simultaneously, rather than being cleanly one or the other. A marker-panel comparison across res 0.2/0.4/0.6 showed res 0.2 gave 7 of 9 clusters with clean, single-category canonical marker signal (T cells, NK/Cytotoxic, B cells, Macrophages, Monocytes/DC, Mast cells, pDC), plus two smaller clusters requiring separate interpretation (see Section 5). Resolutions 0.4 and 0.6 fragmented these same populations without adding distinguishable biology.

**GSE176078: resolution 0.6 used (28 clusters).** This dataset's UMAP structure and cluster-level marker signal remained stable across resolutions and across the seeding fix (26→28 clusters is a modest change, not the qualitative shift seen in GSE114725), so resolution 0.6 — chosen for preserving immune population diversity without over-fragmenting the epithelial compartment — was retained.

*(Figure: `GSE114725_resolution_comparison.png`, `GSE176078_resolution_comparison.png`)*

---

## 4. Marker gene identification

Marker genes per cluster were identified using the Wilcoxon rank-sum test (`sc.tl.rank_genes_groups`, `method="wilcoxon"`), run on `adata.raw` (log-normalised, unscaled, full gene set — not the 2,000-HVG-subsetted, scaled matrix used for clustering itself). Top 10 genes per cluster saved to `results/phase2_clustering_v2/`.

---

## 5. Cell type annotation

Annotation combined three sources, per standard practice for this kind of analysis:

1. **Canonical markers** — a panel covering expected populations (T cell, NK, B cell, macrophage, monocyte/DC, mast cell, pDC, fibroblast/stromal, and dataset-specific epithelial/stromal markers for GSE176078), checked as mean expression per cluster from `adata.raw`.
2. **CellTypist** (`Immune_All_High.pkl` model) as an automated cross-check. Both datasets were run in 10,000-cell batches with `majority_voting=False` — GSE114725 initially crashed under the default full-dataset, `majority_voting=True` path due to a memory allocation failure inside CellTypist's own mean-centering step; batching resolved this and made both datasets' CellTypist methodology consistent.
3. **Published reference labels** (Wu et al. 2021, available for GSE176078 only) — used as an independent cross-tabulation check against final cluster labels.

**Every label was required to be supported by an explicit canonical marker value, not top-DE-genes alone.** This standard was applied retrospectively to catch several labels initially assigned from DE gene inspection without a full marker-panel check (notably GSE114725 clusters 1, 2, and 6, and GSE176078 clusters 12–27), all of which were subsequently re-validated.

### GSE114725 cluster labels (resolution 0.2, 9 clusters)

| Cluster | Label | Basis |
|---|---|---|
| 0 | T cells | TCF7, CD3D/E |
| 1 | CD8/Effector T cells | CD8A=1.43, CD8B=0.79 (highest) |
| 2 | Mixed/stromal-contaminated (CD8+fibroblast signal) | COL1A1=3.05/DCN=2.28 *and* CD8A=0.93/CD8B=0.36 both elevated — likely low-level stromal contamination in an otherwise immune-sorted (CD45+) dataset |
| 3 | NK/Cytotoxic T cells | NKG7/GNLY/PRF1 |
| 4 | B cells | CD79A/MS4A1/CD19 |
| 5 | Macrophages | CD68/LYZ/C1QA/C1QB |
| 6 | pDC | JCHAIN=2.90, IRF7=2.42, IRF8=2.63 |
| 7 | Monocytes/DC | S100A8/S100A9 |
| 8 | Mast cells | CPA3/KIT/TPSB2 |

**Note on the cluster-3-vs-5 macrophage/monocyte question:** at the originally-used resolution 0.6, cluster 5 ("Monocytes/DC") showed *higher* canonical macrophage marker expression (CD68, LYZ, TYROBP) than the cluster labelled "Macrophages" — a genuine inconsistency investigated extensively before the resolution correction. At resolution 0.2, this ambiguity resolved cleanly: cluster 5 (Macrophages) and cluster 7 (Monocytes/DC) are both clean, canonically-distinct populations with no overlap. The ambiguity was an artefact of over-resolved clustering at 0.6, not a real biological question.

### GSE176078 cluster labels (resolution 0.6, 28 clusters)

23 of 28 clusters showed clean, single-category canonical marker evidence (Endothelial, CAFs, PVL, Basal/Luminal epithelial, B cells, T cell subsets, Plasma cells, Macrophages, pDC). Full table in `GSE176078_cluster_annotations_v2_corrected.csv`.

**Known caveats:**
- Clusters 12, 15, 23, 24 labelled "Epithelial (ambiguous)" — ESR1/FOXA1 (luminal-defining transcription factors) absent despite keratin expression; consistent with EMT-like or otherwise non-canonical epithelial states.
- Cluster 25 labelled "Luminal epithelial" with lower confidence — ESR1/FOXA1 values sit at the borderline threshold used to distinguish luminal from ambiguous elsewhere.
- Cluster 21 labelled "Macrophages" with a caveat — highest C1QA of any cluster, but also elevated VIM/COL1A2 (fibroblast/mesenchymal signal), suggesting a possible fibrotic-TAM phenotype rather than a fully clean population.
- **Cluster 21 (original numbering pre-relabel) — n=28 cells — reclassified as "Unassigned (doublet/mixed-identity artefact)".** This cluster showed strong macrophage marker signal (C1QA=3.33, the highest in the dataset) from both our own canonical panel and CellTypist independently, but the published Wu et al. label for this cluster is 100% "Cancer Epithelial" — a three-way disagreement. Given the cluster's negligible size (0.03% of the dataset), this was resolved as a doublet or index-hopping artefact rather than a genuine cell type, consistent with existing project convention for other very small, ambiguous clusters (n<30).

### Reference atlas validation (GSE176078 only)

Final cluster labels were cross-tabulated against Wu et al. 2021's published `celltype_major` annotations. Concordance was strong across nearly all categories (B cells 96.9%, Endothelial 99.5%, PVL 94.9%, Macrophages→Myeloid 98.8%, CAFs 89.3%). The one soft spot — "Epithelial (ambiguous)" matching published labels at only 60.4%, split across Cancer Epithelial/Normal Epithelial/CAFs — corroborates rather than undermines the ambiguous label, since the published atlas shows the same underlying uncertainty. The retracted/reclassified artefact cluster (above) mapping to yet another category (Cancer Epithelial) in the published data further supports treating it as noise rather than a real population.

*(Figure: `GSE114725_annotated_umap_v2.png`, `GSE176078_annotated_umap_split.png` — split into immune and stromal/epithelial panels for readability given the 17-category legend)*

---

## 6. Sub-clustering (stretch goal)

Sub-clustering was applied to GSE114725 only, not GSE176078, as a deliberate choice: GSE114725 is CD45+ sorted with a tumour/blood/normal design, making immune fine-structure the natural focus; GSE176078 is a whole-tumour dataset better suited to the subtype-level composition and reference-validation analysis already performed on it. Adrien's brief specifies sub-clustering "a specific immune cell type of interest" (singular), which this satisfies.

### T cell sub-clustering

T cells, CD8/Effector T cells, and NK/Cytotoxic T cells (n=33,003) were extracted and re-clustered. Resolution 0.5 (5 clusters) was chosen over 0.7 (10 clusters) after the same marker-validation approach used for the parent clustering: 0.5 gave 5 clean populations including a genuine, well-defined Treg cluster; 0.7 fragmented further without adding clearly interpretable biology, plus introduced one small artefact-like cluster.

**CD4/CD8 lineage assignment:** CD4 mRNA has well-documented high dropout in droplet-based scRNA-seq. Lineage was determined primarily by CD8A/CD8B exclusion (low/absent = presumed CD4 lineage), with direct CD4 expression checked as supporting evidence rather than the sole criterion. Both lines of evidence agreed for every cluster.

| Cluster | Label | n | Basis |
|---|---|---|---|
| 0 | CD4 Naive/Resting T cells | 10,807 | TCF7/CCR7/SELL high; CD8A/B low |
| 1 | CD4 Activated T cells | 8,398 | CD69/FOS/JUN elevated; CD8A/B low |
| 2 | Activated CD8 T cells | 6,256 | CD8A=1.64, CD8B=0.92 (highest); CD4 low |
| 3 | Regulatory T cells (Tregs, CD4+) | 2,000 | FOXP3/IL2RA/CTLA4/TIGIT dominant; CD4 highest of any cluster |
| 4 | NK/Cytotoxic T cells | 5,542 | NKG7/GNLY/GZMB/PRF1 dominant; CD8-leaning |

A Treg population (n=2,000, ~6% of T cells) was not resolvable at the parent clustering level and is a genuine finding of the sub-clustering step.

*(Figure: `GSE114725_tcell_annotated_umap.png`)*

### Macrophage sub-clustering

Macrophages (n=5,914; 84% tumour-derived) were extracted and re-clustered. Resolution 0.5 (8 clusters) was chosen: resolution 0.3 missed a genuine population entirely (see below), and resolution 0.7 fragmented two clean populations (LAM-like, Monocyte-like) into weaker/stronger split pairs without adding distinguishable biology, while the smallest population's marker profile was identical at 0.5 and 0.7 — confirming stability rather than an artefact of resolution choice.

| Cluster | Label | n | Top DE genes |
|---|---|---|---|
| 0 | LAM-like macrophages | 946 | SPP1, APOE, FTL, CTSD, CSTB |
| 1 | Lipid-laden/Foam-cell macrophages | 735 | LPL, CSTB, SPP1, CD9, ACP5 |
| 2 | Antigen-presenting macrophages | 910 | HLA-DPB1, HLA-DRA, HLA-DPA1, CD74 |
| 3 | Complement-high macrophages | 1,075 | C3, C1QB, C1QA, C1QC |
| 4 | Resting/Resident macrophages | 469 | SEPP1, C1QA, RNASE1 |
| 5 | Monocyte-like macrophages | 1,265 | VCAN, FCN1, LYZ, S100A9, S100A8 |
| 6 | Non-classical monocytes (CD16+) | 423 | LST1, FCGR3A, CDKN1C, IFITM2 |
| 7 | Unassigned (contamination artefact) | 91 | see below |

**Retracted finding — cytotoxic/NKG7+ macrophage cluster.** An earlier analysis pass identified cluster 7 as a "cytotoxic NKG7+ macrophage" population based on elevated NK-marker expression relative to background (NKG7=0.54, GNLY=0.78, PRF1=0.86 vs ~0.1 elsewhere). On closer marker validation, this cluster's actual top differentially-expressed genes are `DCN, COL1A1, MGP, IGFBP7, HBB` — unambiguous fibroblast and erythrocyte markers, not appearing in any cytotoxic signature. A confirmatory check showed the stromal/RBC markers elevated 17–70× over background (COL1A1 3.36 vs 0.20; HBB 2.02 vs 0.04) — an order of magnitude stronger than the NK-marker elevation (~5–9× over background). Core macrophage identity markers (CD68, LYZ) were *lower* in this cluster than every other macrophage cluster, and QC metrics (elevated `n_genes`/`total_counts`) were consistent with multi-cell or high-ambient-RNA droplets. **Conclusion: cluster 7 represents stromal/erythrocyte contamination, not a genuine cytotoxic macrophage phenotype. The earlier "cytotoxic NKG7+ macrophage" finding is retracted.**

*(Figure: `GSE114725_macrophage_subclusters_umap.png`)*

---

## 7. Cell-type composition and statistical testing (stretch goal)

**GSE114725 — descriptive only.** Per-patient composition showed substantial heterogeneity (e.g. macrophage proportion ranged from 1.9% to 64.7% across 8 patients). A formal paired tumour-vs-blood statistical test was not pursued: only 2 of 8 patients had matched tumour and blood samples, insufficient for a paired test with meaningful power. Composition is reported descriptively.

**GSE176078 — Kruskal-Wallis by subtype (ER+/HER2+/TNBC, n=11/5/10 samples).** No cell type reached significance after FDR correction (all FDR > 0.05), despite descriptive trends consistent with published biology (e.g. macrophage infiltration: TNBC 14.7% vs ER+ 5.3%, uncorrected p=0.030). This is interpreted as a power limitation given the number of cell types tested relative to sample size per group, not evidence against a true biological difference.

**Note on method:** cell-type proportions are compositional data (sum to 100% within each sample), which standard non-parametric tests do not formally account for — an approximation used here for a first-pass analysis; purpose-built compositional methods (e.g. `scCODA`, `propeller`) would be a more rigorous approach if pursued further.

*(Figures: `GSE114725_celltype_proportions.png`, `GSE176078_celltype_proportions_subtype.png`; results tables in `results/phase2_clustering_v2/`)*

---

## 8. Post-review validation (mid-project supervision feedback)

Following supervisor review of the clustering approach, three additional checks were run.

### Tissue-pooling sensitivity analysis

**Concern raised:** GSE114725 clusters were derived by pooling all four tissue compartments
(Blood, Lymphnode, Normal, Tumour) together — could this pooling mask a population that
only exists in one compartment?

**Method:** each tissue was re-clustered completely independently (fresh neighbour graph
and Leiden clustering, same `X_pca_harmony` embedding, resolution 0.2, `random_state=0`),
then every independently-derived cluster was cross-tabulated against the existing pooled
`cell_type` labels to check whether it mapped cleanly onto one label or represented
something new and unmatched.

**Result:** every independently-derived cluster across all four tissues mapped cleanly onto
an existing pooled label (lowest match 71.2% in Blood, most >90%; no cluster showed the
"no dominant match" pattern that would indicate a hidden population).

| Tissue | Cells | Clusters found alone | Clusters found pooled | Lowest match |
|---|---|---|---|---|
| Tumour | 19,594 | 9 | 9 | 92.6% |
| Blood | 15,592 | 5 | 9 | 71.2% |
| Normal | 4,340 | 8 | 9 | 80.8% |
| Lymphnode | 5,136 | 5 | 9 | 89.0% |

Fewer clusters emerge per-tissue in Blood and Lymphnode than in the pooled analysis —
expected, not evidence of masking: these compartments genuinely lack certain tissue-resident
populations (no distinct pDC or Mast cell cluster in Blood; no distinct Macrophage or NK
cluster in Lymphnode), consistent with known immunology for those compartments rather than
a clustering artefact.

**Conclusion, confirmed by supervisor:** pooled clustering is appropriate and was not
hiding tissue-specific structure. Retained as the approach for defining `cell_type`, with
this independent per-tissue analysis kept as supporting evidence. Downstream comparative
analyses (DE, LIANA, scCODA) additionally follow a **pool-to-define, split-to-compare**
principle — see `methods_phase3.md` §1 for the corresponding DE design discussion,
particularly the exclusion of Blood from the Tumour-vs-Normal comparison on both power and
principled grounds.

*(Figure: `GSE114725_tissue_pooling_sensitivity_umap.png`)*

### Expanded marker matrix

Cluster identities re-validated against a broader canonical marker panel spanning 12
lineage categories (including several markers not in the original validation set: TRAC,
GZMK, NCAM1, KLRF1, CD79B, XBP1, DERL3, FCER1A, LILRA4, CLEC4C, LUM, TOP2A), visualised as
a dot plot (cluster × marker, mean expression and fraction of cells expressing). Every
cluster's strongest signal fell exactly on its assigned lineage's marker columns, with no
exceptions. Cluster 2 ("Mixed/stromal-contaminated") showed dots at both T-cell and
stromal-contamination columns, visually confirming the mixed-identity label rather than
suggesting a cleaner single label was missed.

*(Figure: `GSE114725_marker_matrix_dotplot.png`)*

### Per-cluster signature scoring

Each cluster scored (via `sc.tl.score_genes`) against *all* 12 marker panels, not just its
own assigned panel — testing specificity rather than presence. Every cluster's own assigned
panel scored highest specifically within that cluster (see table below), with no exceptions.

| Cluster | Assigned label | Top-scoring panel |
|---|---|---|
| 0 | T cells | Pan-T (anchor) |
| 1 | CD8/Effector T cells | CD8 T |
| 2 | Mixed/stromal-contaminated | Stromal (contam flag) |
| 3 | NK/Cytotoxic T cells | Cytotoxic (shared) |
| 4 | B cells | B cell |
| 5 | Macrophages | Macrophage (vs mono) |
| 6 | pDC | pDC (specific) |
| 7 | Monocytes/DC | Monocyte/DC (vs mac) |
| 8 | Mast cells | Mast |

Cluster 2 scoring highest on the stromal contamination-flag panel is independent, formal
confirmation of that label — not a hedge, but a specific, testable claim that held up.

*(Results: `GSE114725_signature_scores_per_cluster.csv`)*

---

## 9. Fine-resolution annotation (`cell_type_fine`)

Following supervisor review, a second, finer-resolution annotation layer (`cell_type_fine`)
was built alongside the primary `cell_type` column. `cell_type` remains the primary
resolution for all headline findings and is unaffected by anything in this section;
`cell_type_fine` is an additive, optional layer used only where explicitly noted.

### NK/cytotoxic T cell reclassification

**Concern raised:** the original "NK/Cytotoxic T cells" label (both datasets) was validated
using shared cytotoxic-function markers (NKG7, GNLY, PRF1, GZMB) — genes expressed by both
true NK cells and cytotoxic CD8 T cells, since both kill via the same granule machinery.
These markers cannot distinguish the two cell types; only T-cell-receptor genes (CD3D,
CD3E) can.

**Method:** cells in the original "NK/Cytotoxic T cells" cluster (both datasets) were
re-classified using a per-cell marker check: CD3D/CD3E-positive cells were considered
T-lineage; NCAM1/KLRF1/FCGR3A-positive cells were considered NK-lineage. This produces
three groups — True NK (T-marker negative, NK-marker positive), NKT cells (both marker
sets positive — a genuine, distinct innate-like lineage, not an artefact), and Cytotoxic
CD8 T cells (T-marker positive, NK-marker negative, confirmed via elevated CD8A).

**Results:**
- GSE114725: of 4,929 cells, 2,192 (44.5%) True NK, 748 (15.2%) NKT, 843 (17.1%)
  reclassified Cytotoxic CD8. Remaining 1,146 (23.3%) could not be confidently classified.
- GSE176078: of 4,083 cells, 803 (19.7%) True NK, 504 (12.3%) NKT, 1,206 (29.5%)
  reclassified Cytotoxic CD8, 570 (14.0%) retained at the coarse label. GSE176078 showed
  substantially more contamination than GSE114725 (~48% vs ~32% of the original cluster
  reclassified as CD8 or NKT), confirming this correction was necessary in both datasets.

**Unresolved residual cells (GSE114725, n=1,146):** investigated directly rather than
forced into a category. These cells showed (a) lower sequencing depth than the confidently
classified groups (mean 442 genes/804 counts vs 510–586 genes/933–1,052 counts elsewhere),
(b) intact core cytotoxic marker expression at rates matching the cluster overall (85.7%
NKG7-positive vs 89.6% overall), and (c) a top differentially-expressed gene signature
dominated by ribosomal protein genes — a recognised technical signature associated with
lower transcriptional complexity, not a distinct biological state. Three independent
attempts to resolve this group further (single-gene positivity at looser thresholds,
background-corrected composite scoring, raw-sum multi-gene comparison) each failed to
improve on or worsened the original split. **Conclusion: these cells are genuinely
unclassifiable at this resolution given available sequencing depth — excluded from
`cell_type_fine` (set to NaN) rather than force-assigned, while remaining fully valid
in all `cell_type`-level analyses.**

*(Results: `GSE114725_NK_cluster_3way_classification.csv`)*

### CD8 T cell sub-clustering

Following supervisor-provided marker panels (Naive, Central Memory/Effector, Tissue
Resident Memory, Exhausted states), CD8 T cells were sub-clustered for both datasets using
the same method validated for CD4/macrophage sub-clustering (resolution sweep, Harmony-
corrected, seeded).

**Critical prerequisite finding:** the CD8 population pool required combining three
sources per dataset — the original parent-level CD8 cluster, cells reclassified from the
NK/cytotoxic split (above), **and, for GSE114725, 1,090 cells (5.4% of the parent "T
cells" cluster) that were found to be CD8A/CD8B-positive despite being labelled "T cells"
at parent level.** This is a genuine, pre-existing property of the original resolution-0.2
clustering (Section 3) — not a new error — reflecting that coarse-resolution clustering
does not guarantee lineage purity within a correctly-identified broader population. It
was only discovered because CD4 sub-clustering was rebuilt from the single, non-pooled
parent "T cells" cluster (20,213 cells) rather than reusing an earlier sub-clustering
result that had — separately — pooled and reshuffled cells across three parent clusters,
an artefact of that earlier construction rather than of the underlying data.

**Resolution and validation:** resolution 0.5 selected for both datasets (5–6 clusters).
Marker-panel composite scoring proved unreliable where panels shared most of their genes
(Naive and Central Memory panels overlap almost entirely) — ambiguous clusters were
instead resolved via unbiased differential expression (Wilcoxon test, no pre-selected
panel), revealing sharper, more specific signatures than panel scoring alone.

**GSE114725 (5 states, from 9,794 pooled cells):** Activated (4,941, immediate-early/AP-1
signature: CXCR4, DUSP1, CD69, FOS), Naive/Memory (2,246, CCR7/TCF7/SELL-high), Effector
(1,668, GNLY/FGFBP2/GZMH terminal cytotoxic signature), NK-like (840, TYROBP/KLRD1/XCL1-2 —
an innate-like CD8 phenotype, a recognised but less common population), Cycling (99,
STMN1/TUBA1B/HMGB2 proliferation markers).

**GSE176078 (6 states, from 10,916 cells):** Stress-Response/Activated (2,569, heat-shock
protein signature: HSPA1A/B, DNAJB1), Resting/Memory-like (2,391), Effector (2,238,
GZMA/GZMH/NKG7), GZMK+ (1,862), Exhausted (1,232, CXCL13/TIGIT/LAG3/HAVCR2 — a strong,
specific exhaustion signature), Interferon-Response (624, ISG15/MX1/IFIT1/STAT1 — a
type-I interferon-stimulated gene signature not anticipated in the original marker panel,
representing a genuine additional finding).

*(Files: `GSE114725_cd8_subclustered_v2.h5ad`, `GSE176078` fine labels within
`GSE176078_phase2_v2_annotated_corrected_finelabels.h5ad`)*

### Viability for downstream statistical analysis

`cell_type_fine` categories were checked for adequate per-patient/per-sample
representation (≥10 cells/sample, ≥3 viable samples per comparison group — identical
threshold to the original Phase 3 viability checks) before any statistical analysis was
attempted.

- **GSE114725 (Tumour vs Normal):** 7 of 24 fine categories viable, reflecting the same
  small Normal-tissue sample size (n=4 patients) that constrained the original
  `cell_type`-level analysis.
- **GSE176078 (pairwise subtypes):** 20 of 24 fine categories viable for at least one
  comparison — substantially better supported, reflecting the larger sample base (26
  samples vs 8 patients).

`cell_type` was retained as the primary analytical resolution rather than being replaced,
for two independent reasons: (1) all previously reported findings were derived from and
validated against `cell_type`; (2) the viability results above confirm `cell_type_fine`
does not achieve universal coverage and cannot support the same breadth of statistical
testing — the two resolutions serve different analytical purposes (broad, fully-powered
comparison vs specific, sample-permitting sub-lineage resolution) rather than one
superseding the other.

*(Results: `GSE114725_finelabels_viability_check.csv`,
`GSE176078_finelabels_viability_check.csv`)*

---

## 10. Summary of key methodological decisions

1. `n_pcs=30` retained despite an empirical elbow at ~PC 8–12, justified by downstream biological validation rather than PC count alone.
2. Random seeding (`random_state=0`) and single-threaded execution applied throughout to ensure reproducible clustering — without this, the pipeline produced materially different results between runs.
3. GSE114725 resolution changed from 0.6 to 0.2 after the reproducibility fix revealed 0.6 produced mixed-identity clusters; GSE176078 remained at 0.6.
4. All cluster labels required explicit canonical marker support; several were corrected after retrospective review.
5. Two- to three-way triangulation (canonical markers, CellTypist, published reference labels) used wherever possible, catching one mislabelled artefact cluster and confirming the resolution correction.
6. One prior finding (cytotoxic NKG7+ macrophages) was investigated further and retracted after being found to represent contamination rather than genuine biology.
7. Statistical testing on cell-type proportions was applied where sample size allowed and honestly reported as non-significant after correction where that was the case, rather than relying on uncorrected p-values.
8. A supervisor-raised concern about tissue pooling was tested directly (independent per-tissue re-clustering) rather than argued from principle alone, confirming pooled clustering did not mask any tissue-specific population; two further supervisor-requested validation steps (expanded marker matrix, per-cluster signature scoring) both confirmed every existing label with no exceptions.
9. A supervisor comment on NK/cytotoxic marker ambiguity led to a full reclassification (both datasets) and CD8 sub-clustering, uncovering a genuine pre-existing lineage-purity limitation in the original resolution-0.2 clustering (~5% of one parent cluster) — corrected, documented, and folded into a new fine-resolution annotation layer (`cell_type_fine`) that supplements rather than replaces the primary annotation.
10. An honest, unresolved residual (2.3% of one dataset's NK/cytotoxic population) was excluded from fine-resolution analysis rather than force-classified, after three independent classification approaches failed to improve on the original split and root-cause investigation confirmed a genuine sequencing-depth limitation.
