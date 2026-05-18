# Hardware and Memory Notes

## System Specifications

- CPU: AMD Ryzen 5 7520U with Radeon Graphics
- RAM: 8 GB installed (7.2 GB usable)
- Storage: SSD (NVMe)
- Operating System: Windows

## Dataset Loading Tests

The following scRNA-seq datasets were successfully loaded and processed locally using Scanpy:

- GSE176078
- GSE114725

Completed workflows included:
- quality control
- highly variable gene selection
- PCA
- neighbour graph construction
- UMAP
- Leiden clustering
- broad cell-type annotation
- T-cell heterogeneity analysis
- macrophage heterogeneity analysis

## Observations

The laptop was capable of handling moderately large scRNA-seq datasets locally without crashes during preprocessing, clustering, or annotation workflows.

Memory usage became elevated during larger analyses, indicating that extremely large integrated datasets or computationally intensive downstream methods may require optimisation or access to higher-memory compute resources.
