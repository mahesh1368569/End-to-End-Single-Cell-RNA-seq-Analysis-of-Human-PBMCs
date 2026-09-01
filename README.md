# End-to-End Single-Cell RNA-seq Analysis of Human PBMCs

[![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-2C3E50)](https://satijalab.org/seurat/)
[![10x Genomics](https://img.shields.io/badge/10x%20Genomics-3'%20Gene%20Expression-C92A39)](https://www.10xgenomics.com/)
[![Cell Ranger](https://img.shields.io/badge/Cell%20Ranger-count-5B6770)](https://www.10xgenomics.com/support/software/cell-ranger)
[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-yellow.svg)](LICENSE)

## Project overview

This portfolio project demonstrates an end-to-end single-cell RNA-sequencing workflow using a publicly available human peripheral blood mononuclear cell (PBMC) dataset from 10x Genomics. Raw FASTQ reads were processed with Cell Ranger to perform alignment, barcode and UMI processing, cell calling, and gene-level quantification. The resulting filtered feature-barcode matrix was analyzed in R with Seurat for quality control, normalization, dimensionality reduction, graph-based clustering, marker detection, and immune-cell annotation.

The repository is designed to show both command-line/HPC workflow skills and biological interpretation of single-cell transcriptomic data.

