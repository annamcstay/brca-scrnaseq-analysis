# Quality Control and Preprocessing

## Datasets

Two publicly available breast cancer single-cell RNA-seq datasets were analysed:

* GSE114725
* GSE176078

## Cell-level Quality Control

Quality control metrics were calculated using Scanpy, including:

* Number of detected genes per cell (`n_genes_by_counts`)
* Total counts per cell (`total_counts`)
* Percentage mitochondrial reads (`pct_counts_mt`)

Cells were filtered using dataset-specific thresholds selected after visual inspection of QC distributions.

### GSE114725

Cells were retained if:

* `n_genes_by_counts > 200`
* `pct_counts_mt < 10`

### GSE176078

Cells were retained if:

* `n_genes_by_counts > 300`
* `pct_counts_mt < 10`

## Gene Filtering

Genes detected in fewer than three cells were removed using:

```python
sc.pp.filter_genes(adata, min_cells=3)
```

## Doublet Detection

Potential doublets were identified using Scrublet. Cells predicted as doublets were removed before downstream analysis.

## Normalisation

Data were normalised to 10,000 counts per cell using total-count normalisation followed by log-transformation:

```python
sc.pp.normalize_total(target_sum=1e4)
sc.pp.log1p()
```

## Highly Variable Gene Selection

Highly variable genes (HVGs) were identified using the Seurat method implemented in Scanpy:

```python
sc.pp.highly_variable_genes(
    n_top_genes=2000,
    flavor="seurat"
)
```

The top 2,000 HVGs were retained for dimensionality reduction and clustering.

## Saved Data Objects

Two types of processed AnnData objects were retained:

1. Full-gene QC-filtered objects for marker validation and annotation.
2. HVG-restricted Harmony-integrated objects for dimensionality reduction and clustering.

This approach reduced memory requirements while retaining access to canonical marker genes during annotation.

# QC Methods

Datasets were filtered based on:
- number of detected genes
- total counts
- mitochondrial percentage

QC thresholds were adjusted per dataset according to distributional characteristics.
