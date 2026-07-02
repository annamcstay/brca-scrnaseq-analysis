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

## 8. Summary of key methodological decisions

1. `n_pcs=30` retained despite an empirical elbow at ~PC 8–12, justified by downstream biological validation rather than PC count alone.
2. Random seeding (`random_state=0`) and single-threaded execution applied throughout to ensure reproducible clustering — without this, the pipeline produced materially different results between runs.
3. GSE114725 resolution changed from 0.6 to 0.2 after the reproducibility fix revealed 0.6 produced mixed-identity clusters; GSE176078 remained at 0.6.
4. All cluster labels required explicit canonical marker support; several were corrected after retrospective review.
5. Two- to three-way triangulation (canonical markers, CellTypist, published reference labels) used wherever possible, catching one mislabelled artefact cluster and confirming the resolution correction.
6. One prior finding (cytotoxic NKG7+ macrophages) was investigated further and retracted after being found to represent contamination rather than genuine biology.
7. Statistical testing on cell-type proportions was applied where sample size allowed and honestly reported as non-significant after correction where that was the case, rather than relying on uncorrected p-values.
