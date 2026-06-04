# Clustering and Cell Type Annotation
# Input Data
Phase 2 analysis was performed using two processed AnnData objects derived from the Phase 1 quality control pipeline. For GSE114725, the input comprised 36,607 cells and 2,000 highly variable genes following QC filtering, doublet removal, normalisation, HVG selection, scaling, PCA, and Harmony batch correction. For GSE176078, the input comprised 82,931 cells and 2,000 highly variable genes processed through the same pipeline. These objects contained Harmony-corrected PCA embeddings (X_pca_harmony) used for all downstream neighbour graph construction and clustering.
For CellTypist annotation, a separate set of full-gene QC objects was used (GSE114725_phase1_qc_fullgenes.h5ad, GSE176078_phase1_qc_fullgenes.h5ad), which retained all genes passing QC filtering (14,751 genes for GSE114725, 27,221 genes for GSE176078) but had not undergone scaling. These objects were subsetted to the post-doublet-removal cell barcodes present in the Harmony-corrected objects (removing 6 cells from GSE114725 and 102 cells from GSE176078 that had been flagged as doublets), normalised to 10,000 counts per cell and log1p transformed, and used exclusively for CellTypist input. UMAP coordinates and Leiden cluster labels were transferred from the Harmony-corrected objects into these full-gene objects prior to annotation to ensure spatial and cluster consistency across all analyses.

The top 2,000 highly variable genes were selected prior to scaling, PCA, Harmony integration, and Scrublet doublet detection. This reduces computational memory requirements for these steps while retaining the genes carrying the most biological signal. Full-gene expression matrices were retained in separate objects for downstream analyses requiring the complete gene set, including CellTypist annotation and marker gene validation against canonical signatures.

The presence of V(D)J variable region genes in the top 2,000 highly variable genes is a consequence of performing HVG selection prior to V(D)J filtering. These genes exhibit high variability across cells due to somatic recombination diversity within B and T cell clones rather than biologically meaningful cell-type differences. This is a known characteristic of datasets containing lymphocytes and has been reported in published breast cancer single-cell studies. For future analyses, V(D)J genes should be excluded from the gene list prior to HVG selection to prevent them from consuming variable gene slots that could otherwise capture informative biology.
## Dimensionality Reduction

Principal component analysis (PCA) was performed using highly variable genes.

To reduce batch effects between patients/samples, Harmony integration was applied:

* GSE114725: batch variable = `patient`
* GSE176078: batch variable = `orig.ident`

Neighbour graphs and UMAP embeddings were subsequently generated using the Harmony-corrected PCA representation.

## Leiden Clustering

Leiden clustering was performed across multiple resolutions to assess cluster stability and biological interpretability.

### GSE114725

| Resolution | Clusters |
| ---------- | -------- |
| 0.2        | 9        |
| 0.4        | 10       |
| 0.6        | 11       |
| 0.8        | 14       |
| 1.0        | 19       |

### GSE176078

| Resolution | Clusters |
| ---------- | -------- |
| 0.2        | 13       |
| 0.4        | 17       |
| 0.6        | 21       |
| 0.8        | 23       |
| 1.0        | 27       |

Resolution 0.6 was selected for both datasets following visual inspection of UMAP plots across resolutions 0.2–1.0. At resolution 0.2–0.4, clusters were too broad and likely merged biologically distinct populations. At resolution 0.8–1.0, clusters became fragmented with no clear biological basis. Resolution 0.6 produced 11 clusters for GSE114725 and 21 for GSE176078, consistent with the expected diversity of immune and stromal cell types in the breast tumour microenvironment.
## Marker Gene Identification

Cluster marker genes were identified using Scanpy's Wilcoxon rank-sum test:

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
)
```

Top marker genes were exported for manual inspection and cluster annotation.

## Cell Type Annotation

Cluster identities were assigned using:

1. Cluster-specific marker genes.
2. Canonical cell-type markers reported in the literature.
3. Biological interpretation of marker gene combinations.

Annotation remains ongoing and will be validated using automated reference-based annotation tools.

Automated Annotation with CellTypist
Automated cell type prediction was performed using CellTypist (Domínguez Conde et al., 2022) with the Immune_All_High reference model, which covers a broad range of human immune cell types derived from multi-tissue single-cell atlases.
CellTypist requires log1p normalised expression in the input matrix. As the Harmony-corrected objects had undergone scaling (zero-centering and variance normalisation), annotation was performed on the full-gene QC objects normalised to 10,000 counts per cell prior to log1p transformation. UMAP coordinates and Leiden cluster labels were transferred from the Harmony-corrected objects to ensure spatial consistency.
For GSE114725, CellTypist was run with majority voting enabled, which refines per-cell predictions by propagating the dominant label within each local neighbourhood. For GSE176078, memory constraints precluded running the full dataset in a single pass; annotation was therefore performed in sequential batches of 20,000 cells without majority voting, and per-cell predicted labels were concatenated across batches.
Marker Gene Filtering
For GSE176078, initial Wilcoxon marker gene results were dominated in several clusters by immunoglobulin and T-cell receptor variable region genes (prefixes IGHV, IGLV, IGKV, TRAV, TRBV, TRGV, TRDV). These genes reflect somatic V(D)J recombination diversity within B and T cell clones and are not informative for cell type classification. They were excluded from marker gene interpretation, after which cluster-specific markers matched canonical immune and stromal signatures cleanly.
Annotation Validation
CellTypist predictions for GSE176078 were cross-referenced against the published cell type labels provided with the dataset (Wu et al., 2021). Strong concordance was observed for all major populations including T cells, B cells, plasma cells, macrophages, and epithelial cells. The primary sources of disagreement were cancer-associated fibroblasts (CAFs) and perivascular-like cells (PVL), which were distributed across multiple CellTypist categories — consistent with the absence of dedicated stromal categories in the Immune_All_High model. These populations were annotated manually based on marker gene expression (CAFs: COL1A1, COL1A2, LUM; PVL: ACTA2, RGS5, BGN).
