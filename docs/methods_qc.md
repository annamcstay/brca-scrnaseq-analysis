Quality Control and Preprocessing
Datasets
Two publicly available breast cancer single-cell RNA-seq datasets were analysed:

GSE114725 — an immune-enriched breast cancer dataset comprising 47,016 cells across 8 patients
GSE176078 — a full tumour microenvironment dataset comprising 100,064 cells across 26 samples (Wu et al., 2021)

Raw data were loaded from h5ad files and converted to Dask-backed sparse arrays to enable memory-efficient chunk-wise computation throughout the pipeline.
V(D)J Gene Filtering
Prior to quality control, immunoglobulin and T-cell receptor variable region genes were removed from the gene set. These genes (prefixes IGHV, IGLV, IGKV, IGHD, IGHJ, IGLJ, IGKJ, TRAV, TRBV, TRGV, TRDV) exhibit high variability across cells due to somatic V(D)J recombination within B and T cell clones rather than meaningful cell-type differences. Their inclusion in downstream analyses inflates marker gene results without adding annotation value. This filtering removed 24 genes from GSE114725 and 345 genes from GSE176078, leaving 14,851 and 29,388 genes respectively.
Cell-level Quality Control
Quality control metrics were calculated chunk-wise across the Dask array without materialising the full matrix in memory, including:

Number of detected genes per cell (n_genes_by_counts)
Total counts per cell (total_counts)
Percentage mitochondrial reads (pct_counts_mt)

Thresholds were selected after visual inspection of QC distributions for each dataset.
GSE114725 — cells were retained if:

n_genes_by_counts > 200
pct_counts_mt < 15

GSE176078 — cells were retained if:

n_genes_by_counts > 500
pct_counts_mt < 20

The mitochondrial percentage threshold was set more leniently than in standard pipelines (10%) to reflect the naturally elevated mitochondrial expression in immune cell populations and the distributional characteristics of each dataset. The gene count threshold for GSE176078 was set higher (500) given the substantially higher mean genes per cell in this dataset (mean = 1,776) compared to GSE114725 (mean = 611).
Gene Filtering
Genes detected in fewer than three cells were removed:
pythonsc.pp.filter_genes(adata, min_cells=3)
After cell and gene filtering, GSE114725 retained 44,533 cells and 14,828 genes, and GSE176078 retained 92,236 cells and 27,379 genes.
Doublet Detection
Doublets were identified using Scrublet, run independently per patient or sample rather than on the pooled dataset. This per-sample approach is biologically appropriate because doublets can only form within a sample during library preparation, not across samples. Running per sample also avoids the memory constraints of running Scrublet on the full dataset simultaneously.
GSE114725 was processed across 8 patient samples and GSE176078 across 26 samples. Doublet rates per sample ranged from 0.0% to 3.4% in GSE114725 and 0.0% to 0.9% in GSE176078. In total, 675 doublets were removed from GSE114725 and 194 from GSE176078, leaving 43,858 and 92,042 cells respectively.
Normalisation
Data were normalised to 10,000 counts per cell using total-count normalisation followed by log1p transformation:
pythonsc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
The full log-normalised gene matrix was stored in adata.raw before any further dimensionality reduction steps. This ensures that all biological analyses requiring the complete gene set — including automated cell type annotation with CellTypist and canonical marker gene validation — can be performed directly from the stored raw matrix without requiring separate data objects.
Dimensionality Reduction and Batch Correction
(To be completed pending resolution of full-gene PCA approach — see note below)
Batch correction will be performed using Harmony integration. For GSE114725 the batch variable is patient identity; for GSE176078 it is sample identifier (orig.ident). UMAP embeddings will be generated before and after Harmony correction to assess the degree of batch effect and confirm that biological structure is preserved after integration.

Pipeline note: Standard PCA on the full gene matrix exceeds available RAM on the local analysis machine due to internal matrix densification during SVD computation. The approach for this step is under discussion with the project supervisor. All upstream steps (V(D)J filtering, QC, Scrublet, normalisation) were completed on the full gene matrix using Dask-backed chunked computation.
