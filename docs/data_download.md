# Data Download and Processing

---

## GSE176078 — Wu et al. 2021 BRCA scRNA-seq Dataset

### Source
- GEO accession: GSE176078
- Dataset folder used locally:
  `GSE176078_Wu_etal_2021_BRCA_scRNASeq`

### Original Download Format
The dataset was downloaded as a processed sparse count matrix with accompanying barcode, gene, and metadata files.

Files used:
- `count_matrix_sparse.mtx`
- `count_matrix_barcodes.tsv`
- `count_matrix_genes.tsv`
- `metadata.csv`

### Processing Workflow
The sparse matrix was loaded using Scanpy:

```python
adata = sc.read_mtx("count_matrix_sparse.mtx")
```

Gene and barcode files were loaded separately using pandas.

The matrix was originally arranged as genes × cells, so it was transposed to cells × genes:

```python
adata = adata.T
```

Metadata from `metadata.csv` was aligned to the AnnData object using the cell barcode identifiers.

Metadata columns included:
- `orig.ident`
- `nCount_RNA`
- `nFeature_RNA`
- `percent.mito`
- `subtype`
- `celltype_subset`
- `celltype_minor`
- `celltype_major`

A dataset label was added:

```python
adata.obs["dataset"] = "GSE176078"
```

Gene and cell names were made unique before saving.

### Output
The resulting AnnData object contained:

```text
100064 cells × 29733 genes
```

The object was saved locally as:

```text
data/raw/GSE176078_raw.h5ad
```

### Loading Example

```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE176078_raw.h5ad")
```

---

## GSE114725 — Breast Cancer scRNA-seq Dataset

### Source
- GEO accession: GSE114725
- Downloaded from:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725

### Original Download Format
The dataset was downloaded as a corrected raw expression matrix CSV file:

```text
raw_corrected.csv
```

### Processing Workflow
The dataset was loaded into Python using pandas and converted into a sparse matrix representation to reduce memory usage.

Metadata columns included:
- `patient`
- `tissue`
- `replicate`
- `cluster`
- `cellid`

Expression columns were converted into a sparse matrix using `scipy.sparse`.

The sparse matrix and metadata were assembled into an AnnData object using Scanpy.

Cell IDs were used as observation names and duplicate identifiers were removed before saving.

### Output
The resulting AnnData object contained:

```text
47016 cells × 14875 genes
```

The object was saved locally as:

```text
data/raw/GSE114725_raw.h5ad
```

### Loading Example

```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE114725_raw.h5ad")
```
