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


# Data Download and Processing

## GSE114725 — Breast Cancer scRNA-seq Dataset

### Source
- GEO accession: GSE114725
- Downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725

### Original Download Format
- Raw corrected expression matrix CSV file:
  `raw_corrected.csv`

### Processing Workflow
The dataset was loaded into Python using pandas and converted into a sparse matrix representation to reduce memory usage.

Metadata columns:
- patient
- tissue
- replicate
- cluster
- cellid

Expression columns were converted into a sparse matrix using scipy.sparse.

The sparse matrix and metadata were then assembled into an AnnData object using Scanpy.

### Output
The processed AnnData object was saved locally as:

data/raw/GSE114725_raw.h5ad

### Loading Example

```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE114725_raw.h5ad")
