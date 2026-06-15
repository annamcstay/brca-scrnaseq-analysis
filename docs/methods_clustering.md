Input Data
Phase 2 analysis was performed using two processed AnnData objects derived from the Phase 1 v2 QC pipeline. GSE114725 comprised 44,662 cells and 2,000 highly variable genes following QC filtering, doublet removal, normalisation, HVG selection, scaling, PCA, and Harmony batch correction. GSE176078 comprised 91,425 cells and 2,000 highly variable genes processed through the same pipeline. Both objects contained Harmony-corrected PCA embeddings (X_pca_harmony) and full-gene log-normalised counts stored in adata.raw (14,800 genes for GSE114725, 27,343 genes for GSE176078).
Dimensionality Reduction and Clustering
Neighbour graphs were constructed using the Harmony-corrected PCA embedding (use_rep="X_pca_harmony", k=15, 30 PCs) via sc.pp.neighbors. UMAP embeddings were computed with a fixed random seed (random_state=42) for reproducibility. Leiden clustering was performed at resolutions 0.2, 0.4, 0.6, 0.8, and 1.0 using the igraph implementation with n_iterations=2.
Resolution 0.8 was selected for GSE114725 (6 clusters) and resolution 0.6 for GSE176078 (26 clusters) following visual inspection of UMAP plots across all resolutions. The lower cluster count for GSE114725 reflects the immune-enriched nature of the dataset — immune cell populations are transcriptionally more similar than the mixed TME populations in GSE176078.
Marker Gene Identification
Cluster marker genes were identified using Wilcoxon rank-sum tests (sc.tl.rank_genes_groups). All V(D)J variable region and constant region gene prefixes (IGHV, IGLV, IGKV, TRAV, TRBV, TRGV, TRDV, IGHD, IGHJ, IGLJ, IGKJ, TRAC, TRBC, TRGC, TRDC, IGHA, IGHE, IGHG, IGHM, IGLC, IGKC) were excluded from marker gene interpretation as these genes reflect clonal diversity rather than cell type identity.
Cell Type Annotation
Cluster identities were assigned using a combination of three approaches. First, cleaned Wilcoxon marker genes were matched against canonical cell type markers from the literature. Second, automated annotation was performed using CellTypist with the Immune_All_High reference model — annotation was run directly on adata.raw.to_adata() which contains log-normalised unscaled expression values as required by CellTypist, eliminating the need for separate data objects used in the original pipeline. GSE114725 was annotated with majority voting enabled; GSE176078 was processed in batches of 10,000 cells without majority voting due to memory constraints. Third, biological judgement was applied where automated and marker-based annotations disagreed, particularly for stromal populations absent from the Immune_All_High model.
Annotation Validation
CellTypist predictions for GSE176078 were cross-referenced against published cell type labels from Wu et al. 2021. Concordance exceeded 95% for all major immune populations including T cells, B cells, NK cells, macrophages, endothelial cells, and plasma cells. Discrepancies were limited to stromal populations (CAFs, PVL) where the Immune_All_High model lacks dedicated categories, and to small clusters (Cycling myeloid, n=19; Monocytes/DC, n=22) which showed mixed signal consistent with cluster boundary artefacts.
Stretch Tasks
T cell sub-clustering: T cells, NK/Cytotoxic T cells, and Activated T cells from GSE114725 (31,892 cells) were extracted and re-clustered independently using the same Harmony-corrected PCA embedding. Leiden clustering at resolution 0.7 identified 5 sub-populations: NK/Cytotoxic T cells (NKG7, GNLY, PRF1), Resting T cells (ribosomal gene signature), Naive/Memory T cells, Activated T cells (FOS, JUN, DUSP1, NFKBIA), and a small B cell contamination cluster at the cluster boundary.
Cell type proportion analysis: Cell type proportions were quantified per patient for GSE114725 and per breast cancer subtype for GSE176078. Notable findings: BC3 and BC6 showed elevated macrophage/monocyte proportions (>30%) compared to other patients; ER+ tumours were dominated by luminal epithelial cells (39%) with lower immune infiltration; HER2+ tumours showed highest CD8 T cell (15.8%) and naive/memory T cell (24.9%) proportions; TNBC showed highest macrophage (14.9%) and cycling epithelial (6%) proportions.
Reference atlas comparison: Our GSE176078 annotations were compared against the published Wu et al. 2021 cell type labels using cross-tabulation. Agreement exceeded 94% for B cells, basal epithelial, CAFs, endothelial cells, macrophages, PVL, plasma cells, and all T cell populations. Luminal epithelial cells spanned both Cancer Epithelial (85.8%) and Normal Epithelinal (14%) in the published labels, consistent with our use of a single luminal epithelial category at this resolution.


Cluster Annotation Justifications
GSE114725 (leiden_0.8, 6 clusters)
Cluster 0 — T cells (resting)

Marker genes: RPL/RPS ribosomal gene signature dominant, low activation markers. CellTypist confirmed T cell identity. Resting designation based on absence of activation (FOS, JUN) or cytotoxic (NKG7, GNLY) markers.
Cluster 1 — T cells (naive/memory)

Marker genes: CCR7, TCF7, LEF1 — canonical naive/central memory T cell markers. CellTypist confirmed. Distinct from resting cluster by higher CCR7 expression indicating lymph node homing capacity.
Cluster 2 — NK/Cytotoxic T cells

Marker genes: NKG7, GNLY, PRF1, KLRD1, CCL5, CST7. Strong cytotoxic effector signature. CellTypist confirmed NK/cytotoxic T cell identity. KLRD1 (CD94) expression consistent with NK cell fraction.
Cluster 3 — Activated T cells

Marker genes: FOS, JUN, JUNB, DUSP1, NFKBIA, ZFP36 — immediate early gene response signature indicating recent activation. CellTypist confirmed T cell identity. Immediate early gene signature consistent with recently stimulated T cells.
Cluster 4 — Macrophages

Marker genes: LYZ, CD68, TYROBP, AIF1, CD74. Classic myeloid/macrophage markers. CellTypist confirmed myeloid identity. CD68 and TYROBP co-expression confirms macrophage rather than monocyte identity.
Cluster 5 — Monocytes/DC

Marker genes: S100A8, S100A9, VCAN, FCN1. Classical monocyte markers. CellTypist confirmed. S100A8/A9 alarmin expression and VCAN distinguish from tissue-resident macrophages.

GSE176078 (leiden_0.6, 26 clusters)
Clusters 0-1 — Endothelial cells

Marker genes: PECAM1, CD34, VWF. Canonical endothelial markers. Two clusters reflect endothelial heterogeneity consistent with Wu et al. 2021 (stalk-like vs tip-like states).
Cluster 2 — CAFs

Marker genes: COL1A1, PDGFRA, FAP, ACTA2. Cancer-associated fibroblast markers. Consistent with Wu et al. 2021 CAF populations.
Cluster 3 — PVL

Marker genes: MCAM, PDGFRB, RGS5, ACTA2. Perivascular-like cell markers. Distinct from CAFs by MCAM (CD146) expression.
Cluster 4 — Basal epithelial

Marker genes: KRT5, KRT14, TP63. Basal epithelial markers. Validated against Wu et al. 2021 Normal Epithelial category (94.8% concordance).
Cluster 5 — B cells

Marker genes: MS4A1, CD79A, CD19. Canonical B cell markers. 96.9% concordance with Wu et al. B-cells category.
Cluster 6 — Cycling cells

Marker genes: MKI67, TOP2A, PCNA. Proliferation markers. Mixed lineage cycling population — 96.7% map to T-cells in Wu et al. suggesting predominantly cycling T cells.
Cluster 7 — Plasma cells

Marker genes: JCHAIN, MZB1, IGKC. Plasma cell/plasmablast markers. 96.1% concordance with Wu et al. Plasmablasts category.
Cluster 8 — Cycling epithelial

Marker genes: MKI67 + epithelial markers. 98.5% concordance with Wu et al. Cancer Epithelial category.
Cluster 9 — CD8 T cells

Marker genes: CD8A, CD8B, GZMK. 99.9% concordance with Wu et al. T-cells category.
Cluster 10 — NK cells

Marker genes: NKG7, GNLY, KLRD1, NCAM1. 100% concordance with Wu et al. T-cells category (Wu et al. do not separate NK at major lineage level).
Cluster 11 — T cells

Marker genes: CD3D, CD3E, IL7R. 99.8% concordance with Wu et al. T-cells category.
Cluster 12 — Naive/memory T cells

Marker genes: CCR7, TCF7, SELL. 99.6% concordance with Wu et al. T-cells category.
Cluster 13 — Luminal epithelial

Marker genes: KRT8, KRT18, EPCAM. Luminal epithelial markers. 85.8% Cancer Epithelial + 14% Normal Epithelial in Wu et al. — single label spans both cancer and normal luminal cells at this resolution.
Cluster 14 — Macrophages

Marker genes: CD68, LYZ, TYROBP, AIF1. 98.8% concordance with Wu et al. Myeloid category.
Cluster 15 — Monocytes/DC

Marker genes: S100A8, S100A9, VCAN. 59.1% Myeloid + 31.8% CAFs in Wu et al. — small cluster (22 cells) likely represents a cluster boundary artefact.
Cluster 16 — Cycling myeloid

Marker genes: MKI67 + myeloid markers. Only 19 cells — likely artefactual cluster, noted as limitation.
Cluster 17 — pDC

Marker genes: IRF7, LILRA4, CLEC4C. Plasmacytoid dendritic cell markers. 99.7% concordance with Wu et al. Myeloid category.
Clusters 18-25 — Luminal epithelial / Epithelial variants

Marker genes: KRT8, KRT18, EPCAM variants. Multiple luminal and epithelial clusters reflecting tumour epithelial heterogeneity consistent with Wu et al. 2021.
