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

### Nasal airway epithelium brushing bulk RNA-seq analysis 
* [WGCNA.R](1_Analysis%20on%20GALAII%20cohort/WGCNA.R) - performs weighted gene coexpression network analysis on 695 RNA-seq samples
* [T2_and_interferon_analysis.R](1_Analysis%20on%20GALAII%20cohort/2_T2_and_interferon_analysis.R) 
  * Association analyses between T2 related module eigengenes (MEpink and MEsaddlebrown) and expression of ACE and TMPRSS2 (Fig. 2b-g)
  * Association analyses between interferon related module eigengenes (MEtan and MEpurple) and expression of TMPRSS2 (Fig. 4a,b,d)
  * Hierarhical clustering based on the tan and saddlebrown modules (SFig. 1a, 2a)
* [Multivariate_Regression_Analysis.R](1_Analysis%20on%20GALAII%20cohort/3_Multivariate_Regression_Analysis.R) - performs multivariate regression analysis on the following models (Table 1):
  * ACE2 expr. ~ age + interferon_status + type2_status + sex + Asthma_Status + rs181603331
  * TMPRSS2 expr. ~ age + interferon_status + type2_status + sex + Asthma_Status + rs1475908 + rs74659079 + rs2838057

### Analysis of scRNA-seq data from the nasal epithelial brushing
* 

### Analysis of bulk RNA-seq data from IL-13 and HRV infected ALI nasal airway epithelial cultures


### Analysis of scRNA-seq data from 10 day IL-13-stimulated and control tracheal cell ALI cultures
