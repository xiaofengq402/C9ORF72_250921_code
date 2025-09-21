Single-Cell RNA Sequencing Analysis Project
This repository contains the code for the analysis and visualization of single-cell RNA sequencing data as presented in our study.

Project Structure & Script Descriptions
File Name	Description
check.ipynb	The main pipeline for comprehensive single-cell analysis. This notebook includes data loading, quality control, normalization, clustering, and cell type annotation.
heatmap_total.ipynb	Code for generating all heatmaps presented in the manuscript. This includes differential gene expression heatmaps, marker gene heatmaps, and other complex visualizations.
ralated.ipynb	Scripts for performing correlation analysis and generating corresponding correlation plots (e.g., gene-gene correlation).
LC3_motif.ipynb	Analysis and visualization of the LC3 motif prediction results. This notebook generates the figures related to motif patterns.
SPR_analysis.ipynb	Scripts to load, parse, and statistically analyze the results from Surface Plasmon Resonance (SPR) experiments. Also generates the resulting publication-quality plots.
ggplot.ipynb	A dedicated notebook for creating customized scatter plots and dot plots (e.g., for correlation data) using the ggplot2 framework in Python (plotnine) or R.
