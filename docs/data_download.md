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
