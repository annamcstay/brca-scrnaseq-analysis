# Clustering and Cell Type Annotation

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
