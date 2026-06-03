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
