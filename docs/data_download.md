# Data Download and Storage

## Dataset 1 — GSE176078

### Source

* GEO accession: GSE176078
* Downloaded from NCBI GEO:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

### Data Type

* Breast cancer tumour microenvironment scRNA-seq dataset

### Downloaded Files

* Processed expression matrix / AnnData-compatible files

### Local Storage

Stored locally in:

```text
data/raw/GSE176078_raw.h5ad
```

### Loading Method

Loaded into Python using Scanpy:

```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE176078_raw.h5ad")
```

### Notes

Intermediate processed datasets are stored in:

```text
data/processed/
```

including QC-filtered, clustered, and annotated objects.

---

## Dataset 2 — GSE114725

### Source

* GEO accession: GSE114725
* Downloaded from NCBI GEO:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725

### Data Type

* Breast cancer single-cell RNA-seq dataset

### Downloaded Files

* Processed expression matrix / AnnData-compatible files

### Local Storage

Stored locally in:

```text
data/raw/GSE114725_raw.h5ad
```

### Loading Method

Loaded into Python using Scanpy:

```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE114725_raw.h5ad")
```

### Notes

Processed and annotated datasets are stored in:

```text
data/processed/
```

