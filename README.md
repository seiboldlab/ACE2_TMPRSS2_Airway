---
Title: "Type 2 and interferon inflammation strongly regulate SARS-CoV-2 related gene expression in the airway epithelium"
Author: Satria P. Sajuthi, Peter DeFord
---

This repository contains scripts used to produce the figures within the paper "Type 2 and interferon inflammation strongly regulate SARS-CoV-2 related gene expression in the airway epithelium". They are organized as follow:
1. preliminary analyses.
2. main population analyses.
3. population subset analyses. 
4. scRNA.
The necessary packages to load can be found at the top of each script. We used R version 3.5.1 for all analyses. 

### Bulk RNA-seq analysis on GALA II
1. Preprocessing of RNA-seq data
Raw sequencing reads were trimmed using skewer37 (v0.2.2) with the following parameter settings: end-quality=15, mean-quality=25, min=30. Trimmed reads were then aligned to the human reference genome GRCh38 using GSNAP38 (v20160501) with the following parameter settings: max-mismatches=0.05, indel-penalty=2, batch=3, expand-offsets=0, use-sarray=0, merge-distant-same-chr. Gene quantification was performed with htseq-count39 (v0.9.1) using iGenomes GRCh38 gene transcript model. Variance stabilization transformation (VST) implemented in DESeq240 (v1.22.2) was then carried out on the raw gene count matrix to create a variance stabilized gene expression matrix suitable for downstream analyses. 
