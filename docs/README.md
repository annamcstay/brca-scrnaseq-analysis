# Integrative Analysis of Single-Cell RNA-Seq Data to Map Immune Cell Heterogeneity in Breast Cancer Microenvironments

## Project Overview

This MSc project investigates immune cell heterogeneity within breast cancer tumour microenvironments using publicly available single-cell RNA sequencing (scRNA-seq) datasets.

The project focuses on identifying major tumour microenvironment populations and characterising immune subpopulations, particularly T-cell and macrophage heterogeneity, using computational single-cell analysis workflows in Python with Scanpy.

## Objectives

- Perform quality control and preprocessing of scRNA-seq datasets
- Construct tumour microenvironment cell atlases
- Identify major immune and stromal populations
- Characterise T-cell heterogeneity
- Characterise macrophage heterogeneity
- Explore immune activation and suppression states within tumours

## Datasets

Current shortlisted datasets:

- GSE176078
- GSE114725

Further dataset information is provided in `docs/datasets.md`.

## Analysis Workflow

1. Data loading
2. Quality control and preprocessing
3. Highly variable gene selection
4. PCA and dimensionality reduction
5. Clustering and UMAP visualisation
6. Broad cell-type annotation
7. Immune heterogeneity analysis
   - T-cell states
   - Macrophage states

## Tools and Environment

Analysis performed using:

- Python
- Scanpy
- AnnData
- NumPy
- Pandas
- Matplotlib
- Seaborn

Reproducible environment specified in `environment.yml`.

## Repository Structure

```text
docs/          Project documentation
figures/       Saved figures and plots
Notebooks/     Jupyter notebooks for analysis workflow
