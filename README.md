---
Title: "Type 2 and interferon inflammation strongly regulate SARS-CoV-2 related gene expression in the airway epithelium"
Author: Satria P. Sajuthi, Peter DeFord
---

This repository contains scripts used to produce the figures within the paper "Type 2 and interferon inflammation strongly regulate SARS-CoV-2 related gene expression in the airway epithelium". They are organized as follow:
1. [Nasal airway epithelium brushing bulk RNA-seq analysis](#Nasal-airway-epithelium-brushing-bulk-RNA-seq-analysis).
2. [Analysis of scRNA-seq data from the nasal epithelial brushing](#Analysis-of-scRNA-seq-data-from-the-nasal-epithelial-brushing).
3. [Analysis of bulk RNA-seq data from IL-13 and HRV infected ALI nasal airway epithelial cultures](#Analysis-of-bulk-RNA-seq-data-from-IL-13-and-HRV-infected-ALI-nasal-airway-epithelial-cultures). 
4. [Analysis of scRNA-seq data from 10 day IL-13-stimulated and control tracheal cell ALI cultures](#Analysis-of-scRNA-seq-data-from-10-day-IL-13-stimulated-and-control-tracheal-cell-ALI-cultures).

The necessary packages to load can be found at the top of each script. We used R version 3.5.1 for all analyses. 

## __Nasal airway epithelium brushing bulk RNA-seq analysis__ 
* [WGCNA.R](1_Analysis%20on%20GALAII%20cohort/WGCNA.R) - performs weighted gene coexpression network analysis on 695 RNA-seq samples
* [T2_and_interferon_analysis.R](1_Analysis%20on%20GALAII%20cohort/2_T2_and_interferon_analysis.R) 
  * Association analyses between T2 related module eigengenes (MEpink and MEsaddlebrown) and expression of ACE and TMPRSS2 (Fig. 2b-g)
  * Association analyses between interferon related module eigengenes (MEtan and MEpurple) and expression of TMPRSS2 (Fig. 4a,b,d)
  * Hierarchical clustering based on the tan and saddlebrown modules genes (SFig. 1a, 2a)
* [Multivariate_Regression_Analysis.R](1_Analysis%20on%20GALAII%20cohort/3_Multivariate_Regression_Analysis.R) - performs multivariate regression analysis on the following models (Table 1):
  * ACE2 expr. ~ age + interferon_status + type2_status + sex + Asthma_Status + rs181603331
  * TMPRSS2 expr. ~ age + interferon_status + type2_status + sex + Asthma_Status + rs1475908 + rs74659079 + rs2838057
* [inVivo_virus_Diff_Gene_Expression.R](1_Analysis%20on%20GALAII%20cohort/4_inVivo_virus_Diff_Gene_Expression.R) - Performs differential expression analysis between virus infected individuals and uninfected individuals for two sets of viruses:
  * Human coronaviruses
  * Pooled human rhinovirus C, Influenza A, Influenza B, Orthopneumovirus, Metapneumovirus, Enterovirus, and Parainfluenza.
* [IPA_virus_expression.R](1_Analysis%20on%20GALAII%20cohort/5_IPA_virus_expression.R) - Compiles the output of IPA to produce a heatmap comparing the expression of canonical cytotoxic pathways in the two groups of virus infected individuals.
* [viral_GSEA.R](1_Analysis%20on%20GALAII%20cohort/6_viral_GSEA.R) - Performs Gene Set Enrichment Analysis of the differntially expressed genes that are specific to each group of viruses, or the shared DEGs.
  * Rank orders genes based on the relative expression in immune cell types vs other immune cells.
  * GSEA of Shared DEGs, CoV+ specific DEGs, or DEGs specific to the other viruses on these pre-ranked gene sets.

## __Analysis of scRNA-seq data from the nasal epithelial brushing__
* [Nasal_brush_10X_analysis.R](2_Nasal_brush_scRNA-seq/Nasal_brush_10X_analysis.R) - carries out analyses on the nasal brushing 10X scRNA-seq data including clustering, UMAP visualization of cell populations (Fig 1a), find marker genes among cell populations, violin plots of ACE2 and TMPRSS2 (Fig. 1b-c) 

## __Analysis of bulk RNA-seq data from IL-13 and HRV infected ALI nasal airway epithelial cultures__
* [HRV_IL13_ALI_analysis.R](3_HRV_IL13_in_vitro_ALI/HRV_IL13_ALI_analysis.R) - carries out differential expression analyses between paired IL-13-stimulated and control samples (N = 5 donors, Fig. 3b,c) and between paired HRV-infected and control samples (N = 5 donors, Fig. 4g,h)

## __Analysis of scRNA-seq data from 10 day IL-13-stimulated and control tracheal cell ALI cultures__
* [IL-13_singleCell.R](4_IL13_scRNA-seq/10X_IL-13_singleCell.R) - carries out analyses on the tracheal 10X scRNA-seq data including clustering, UMAP visualization of cell populations (Fig 3d), find marker genes among cell populations, differential expression analysis between unstim. and IL-13 stimulated sample (Fig. 3e,f) 
