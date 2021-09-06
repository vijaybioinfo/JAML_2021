## JAML immunotherapy preferentially targets tumor-infiltrating tissue-resident memory T cells

This repository contains the scripts used to analyzed our asthma airways bulk and single-cell data.

Developer: Ciro Ramírez-Suástegui (ciro@lji.org)

Vijayanand Lab
Division of Vaccine Discovery
La Jolla Institute for Immunology
La Jolla, CA 92037, USA

Global description
---

*Demultiplexing libraries*: Cell Ranger was used to demultiplex the 10x libraries. And our in-house
mapping [pipeline](https://github.com/ndu-UCSD/LJI_RNA_SEQ_PIPELINE_V2) was used for the bulk data.

*Quality control*: An in-house [script](https://github.com/vijaybioinfo/quality_control)
was used to explore the quality of the data and select the thresholds.

*Doublets*: Scrublet was used to detect doublets in the single-cell data.

*Clustering*: Seurat was used to cluster the data.

For more specific information about the data generation and processing, please check the methods.

Code
---

### Single-cell RNA-seq integration analyses

*Download*
The main script to download the data is `download_sc_data.R`, but not all
af the studies provide an easy accessible dataset. So, there are other ways
some of the data were downloaded in `/home/ciro/simon/scripts/downloads.R`

*Cell selection/filter*
`integration/select_cells.R` contains the instructions to subset the datasets.

*Feature/gene selection*
Before we can integrate the data we need to make the features compatible
across the datasets. They are not easily matched in part because of the gene
annotation utilized in the each study. You can accomplish this by using
`integration/select_features.R`. This was explored in
`integration/gene_compatibility.R`.

*Integration with Seurat*
We do the alignment with Seurat in `integration/clustering.R`. After a
round of analysis, some clusters were filtered out with
`integration/clustering_filter.R`.

*Further analyses*
`integration/clustering_further_analyses.R` shows a few visualizations useful
for various papers.

#### SiEs05 TCR analysis
File: SiEs05_tcr_analysis.R

It also includes some SiEs012 viz.

Citation
---
Please cite the following manuscript if you are using this repository:

*Under review*

Contact
---
Please email Ciro Ramírez-Suástegui (ciro@lji.org) and Vijayanand Pandurangan (vijay@lji.org).
