#!/usr/bin/env Rscript

"
This script performs the initial differential expression analysis 
between CoV+ and Severe+ each vs. uninfected individuals.

INPUTS:
- Phenotype data to use as covariates
- Raw count matrix
- VST normalized count matrix
- Virus Identification Pipeline results by Virus type
- WGCNA module eigengenes
- WGCNA module assignments by gene

OUTPUTS:
- Limma results for CoV+ vs Uninfected (Corona_v_uninfected.DEG.txt)
- Limma results for Severe+ vs Uninfected (Severe_v_uninfected.DEG.txt)
- Combined differential expression results (COV_Severe_Combined_DEGs.txt)
- DEG list with annotation of which infection group genes were significant
  (Combined_sigDEGs.txt)
- Panels for Figure 6:
  - Tan module eigengenes by infection status:
    Fig6A_tan_eigen_byInfection.pd
  - Purple module eigengenes by infection status:
    Fig6B_purple_eigen_byInfection.pdf
  - Scatterplot of differentially expressed genes in the two viruses:
    Fig6c_Severe_logFC_comp.pdf
  - Enrichment terms for shared genes, from Fig 6C:
    bothup_top_enrichr_Severe.txt
    bothdn_top_enrichr_Severe.txt
- Enrichment scores for upregulated DEGs in the purple and tan modules:
    (CoV_Severe_DEG_WGCNA_enrichment.txt)
"

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(enrichR))


#=================#
# LOAD INPUT DATA #
#=================#

##### Clinical data #####
phen = read.table("/Seibold/proj/GALA_720donors/ANALYSIS_FIXED_PHENO/TABLES/695.GALA.final.phen.with.add.txt", header=T, sep='\t', strings=F)
phen$asthma_status = factor(plyr::mapvalues(phen$control, from=c(0,1), to=c("asthmatic","HC")), levels=c("HC", "asthmatic"))
rownames(phen) = phen$SubjectID

##### Raw expression count matrix #####
counts = read.table("/Seibold/proj/GALA_720donors/ANALYSIS_FIXED_PHENO/TABLES/695_raw_counts.txt", sep='\t', header=T, strings=F, row.names=1)

##### Normalized expression data #####
vst = read.table("/Seibold/proj/GALA_720donors/ANALYSIS_FIXED_PHENO/TABLES/695_expr_vst.txt", sep='\t', header=T, strings=F, row.names=1)

##### Virus Finder Results #####
virus_counts = read.xlsx("DATA/vf2_new_merged_mas_edit_5.15.20.xlsx", rowNames=TRUE, sheet=1)
exclude = c("HR1872","HR1751","HR1801")
virus_counts = virus_counts[!rownames(virus_counts) %in% exclude,]


# Convert old Donor IDs to new ones
renames = list(
  HR5519 = "HR5579",
  HR1241 = "SJ1241",
  HR5125 = "HR1212",
  HR1251 = "HR1250",
  HR5398 = "SJ5398"
)

oldnames = rownames(virus_counts)
for (name in names(renames)) {
  if (name %in% oldnames){
    rownames(virus_counts)[which(oldnames == name)] = renames[[name]]
  }
}


#========================================#
# REFORMAT INPUT DATA FOR USE WITH LIMMA #
#========================================#

##### Identify infection groups #####
# Extract CoV+ individuals, and enforce high enough read counts
corona = rownames(virus_counts[virus_counts$For.Use == 1 & !is.na(virus_counts$For.Use), ])
corona = corona[corona %in% colnames(counts)]

# Extract individuals with a more severe virus, and enforce high enough read counts
virus_severe = rownames(virus_counts[virus_counts$For.Use == 2 & !is.na(virus_counts$For.Use), ])
virus_severe = virus_severe[virus_severe %in% colnames(counts)]

# Extract uninfected individuals
all_individuals = colnames(counts)
uninfected = all_individuals[!all_individuals %in% rownames(virus_counts)]


##### Combine CoV+ w/ uninfected #####
samples = c(corona, uninfected)
exp_design = data.frame(
    SubjectID = samples,
    Infection = factor(ifelse(samples %in% uninfected, "Uninfected", "Coronavirus"), levels=c("Uninfected", "Coronavirus")),
    age=phen[samples, "age"],
    gender=factor(phen[samples, "Male"]),
    asthma_status = phen[samples, "asthma_status"],
    stringsAsFactors=FALSE
    )
rownames(exp_design) = exp_design$SubjectID

##### Combine Severe+ w/ uninfected #####
samples2 = c(virus_severe, uninfected)
exp_design2 = data.frame(
    SubjectID = samples2,
    Infection = factor(ifelse(samples2 %in% uninfected, "Uninfected", "SevereVirus"), levels=c("Uninfected", "SevereVirus")),
    age=phen[samples2, "age"],
    gender=factor(phen[samples2, "Male"]),
    asthma_status = phen[samples2, "asthma_status"],
    stringsAsFactors=FALSE
    )
rownames(exp_design2) = exp_design2$SubjectID



##### Filter for genes with high expression #####
# Identify genes with enough representation in the dataset
good_genes = rowSums(counts >= 10) >= (ncol(counts) * 0.15)

# Subset the count matrix
counts_sub = counts[good_genes,samples]
counts_sub2 = counts[good_genes,samples2]

#======================================================#
# Differential expression analysis between CoV vs Ctrl #
#======================================================#

dge = DGEList(counts=counts_sub)

design = model.matrix(~ Infection + age + gender + asthma_status, data=exp_design)
keep = filterByExpr(dge, design)
dge = dge[keep, keep.lib.sizes=F]

dge = calcNormFactors(dge)

# Compute weights to handle class imbalance
logCPM = cpm(dge, log=TRUE)
weights = arrayWeights(logCPM, design, var.group=exp_design$Infection)

v = voom(dge, design, plot=F, weights=weights)

fit_CoV = lmFit(v, design, weights=weights)
fit_CoV = eBayes(fit_CoV)

res_fit_CoV = topTable(fit_CoV, coef="InfectionCoronavirus",  n=Inf, sort="p")
head(res_fit_CoV)

write.table(res_fit_CoV, file="./OUTPUT/Corona_v_uninfected.DEG.txt", sep='\t', quote=F, col.names=NA)

#======================================================#
# Differential expression analysis between Severe vs Ctrl #
#======================================================#

dge = DGEList(counts=counts_sub2)
design = model.matrix(~ Infection + age + gender + asthma_status, data=exp_design2)
keep = filterByExpr(dge, design)
dge = dge[keep, keep.lib.sizes=F]

dge = calcNormFactors(dge)

logCPM = cpm(dge, log=TRUE)
weights = arrayWeights(logCPM, design, var.group=exp_design2$Infection)

v = voom(dge, design, plot=F, weights=weights)

fit_Severe = lmFit(v, design, weights=weights)
fit_Severe = eBayes(fit_Severe)

res_fit_Severe = topTable(fit_Severe, coef="InfectionSevereVirus",  n=Inf, sort="p")
head(res_fit_Severe)

write.table(res_fit_Severe, file="./OUTPUT/Severe_v_uninfected.DEG.txt", sep='\t', quote=F, col.names=NA)


#============================================#
# Compare the output of the two sets of DEGs #
#============================================#

##### Combine the output of the DEG analyses: CoV vs Severe #####
common_genes_Severe = intersect(rownames(res_fit_Severe), rownames(res_fit_CoV))
common_res_Severe = cbind(res_fit_CoV[common_genes_Severe, c("logFC", "AveExpr", "adj.P.Val")], res_fit_Severe[common_genes_Severe,c("logFC", "AveExpr", "adj.P.Val")])
colnames(common_res_Severe) = c("logFC_CoV", "AveExpr_CoV", "adj.P.Val_CoV", "logFC_Severe", "AveExpr_Severe", "adj.P.Val_Severe")

write.table(common_res_Severe, file="./OUTPUT/COV_Severe_Combined_DEGs.txt", sep='\t', quote=F, col.names=NA)


#### Create DEG table for the supplement
sig_cov = common_res_Severe$adj.P.Val_CoV < 0.05 & abs(common_res_Severe$logFC_CoV) > 0.5
sig_other = common_res_Severe$adj.P.Val_Severe < 0.05 & abs(common_res_Severe$logFC_Severe) > 0.5
sig_status = ifelse(sig_cov & sig_other, "Both",
             ifelse(sig_cov, "CoV",
             ifelse(sig_other, "ORV", "Neither")))
common_res_sig = common_res_Severe[sig_cov | sig_other,]
common_res_sig$Significance = sig_status[sig_cov | sig_other]
common_res_sig = common_res_sig[order(common_res_sig$Significance, common_res_sig$adj.P.Val_CoV),]
colnames(common_res_sig) = c("logFC_CoV", "AveExpr_CoV", "adj.P.Val_CoV", "logFC_Other", "AveExpr_Other", "adj.P.Val_Other", "Significance")
write.table(common_res_sig, file="./OUTPUT/Combined_sigDEGs.txt", sep="\t", quote=F, col.names=NA)

##### Identify significant DEGs ####
# FDR of 0.05 and absolute logFC of 0.5
CoV_degs = rownames(res_fit_CoV[(res_fit_CoV$adj.P.Val < 0.05) & (abs(res_fit_CoV$logFC) > 0.5),])
Severe_degs = rownames(res_fit_Severe[(res_fit_Severe$adj.P.Val < 0.05) & (abs(res_fit_Severe$logFC) > 0.5),])
# Find the common and union sets of DEGs
common_degs_Severe = intersect(CoV_degs, Severe_degs)
all_degs_Severe = unique(c(CoV_degs, Severe_degs))

##### Determine direction and membership #####
gene_sets_Severe = data.frame(
    CoV_up = ifelse((all_degs_Severe %in% CoV_degs) & (res_fit_CoV[all_degs_Severe, "logFC"] > 0), 1, 0),
    CoV_down = ifelse((all_degs_Severe %in% CoV_degs) & (res_fit_CoV[all_degs_Severe, "logFC"] < 0), 1, 0),
    Severe_up = ifelse((all_degs_Severe %in% Severe_degs) & (res_fit_Severe[all_degs_Severe, "logFC"] > 0), 1, 0),
    Severe_down = ifelse((all_degs_Severe %in% Severe_degs) & (res_fit_Severe[all_degs_Severe, "logFC"] < 0), 1, 0)
    )
rownames(gene_sets_Severe) = all_degs_Severe
gene_sets_Severe$CoV = (gene_sets_Severe$CoV_up + gene_sets_Severe$CoV_down)
gene_sets_Severe$Severe = (gene_sets_Severe$Severe_up + gene_sets_Severe$Severe_down)

#==============================#
# Generate panels for Figure 6 #
#==============================#

uninfected_color = "#009E73"
CoV_color = "#FFC20A"
OTHER_color = "#D41159"
both_color = "#0C7BDC"
neither_color="gray"

##### Load additional data #####
# WGCNA Results
eigengenes = read.table("/Seibold/proj/GALA_720donors/ENDOTYPING/WGCNA_PAM/TABLES/695.WGCNA.ME3.txt")
wgcna_modules = read.table("/Seibold/proj/GALA_720donors/ENDOTYPING/WGCNA_PAM/TABLES/695.WGCNA.gene2module.with.cor3.txt",
                           header=T, stringsAsFactors = F)
rownames(wgcna_modules) = wgcna_modules$gene

# Reformat data to be easier for plotting
qc_design_Severe_samples_Severe = c(uninfected, corona, virus_severe)
qc_design_Severe = data.frame(
    SubjectID = qc_design_Severe_samples_Severe,
    Infection = factor(ifelse(qc_design_Severe_samples_Severe %in% corona, "CoV", ifelse(qc_design_Severe_samples_Severe %in% virus_severe, "Severe", "Uninfected")), levels=c("Uninfected", "Severe", "CoV")),
    stringsAsFactors=FALSE
    )

##### Panel A: Tan module eigengene by infection status #####
qc_design_Severe$tan  = unlist(eigengenes[qc_design_Severe$SubjectID, "MEtan"])

# Compute t-test for infected vs. uninfected
tan_Severe_t = t.test(y = qc_design_Severe[qc_design_Severe$Infection == "Uninfected", "tan"], x=qc_design_Severe[qc_design_Severe$Infection == "Severe", "tan"])
tan_CoV_t = t.test(y = qc_design_Severe[qc_design_Severe$Infection == "Uninfected", "tan"], x=qc_design_Severe[qc_design_Severe$Infection == "CoV", "tan"])

# Plot boxplot, with t-test values
pdf("./IMAGES/Fig6A_tan_eigen_byInfection_Severe.pdf", width=180/25.4/2, height=180/25.4/2)
ggboxplot(qc_design_Severe, x="Infection", y="tan", add="dotplot", 
          add.params=list(binwidth=diff(range(qc_design_Severe[,"tan"]))/100, alpha=0.5), 
          outlier.shape=NA, main="Interferon response", ylab="Eigengene") +
    annotate("text", x=1, y=max(qc_design_Severe[,"tan"])*1.1, size = 3,
             label="dEg:\np-value:") +
    annotate("text", x=3, y=max(qc_design_Severe[,"tan"])*1.1, size = 3,
             label=sprintf("%0.3f\n%0.1e", diff(range(tan_CoV_t$estimate)), tan_CoV_t$p.value)) +
    annotate("text", x=2, y=max(qc_design_Severe[,"tan"])*1.1, size = 3,
             label=sprintf("%0.3f\n%0.1e", diff(range(tan_Severe_t$estimate)), tan_Severe_t$p.value)) +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

##### Panel B: Purple module eigengene by infection status #####
qc_design_Severe$purple  = unlist(eigengenes[qc_design_Severe$SubjectID, "MEpurple"])

# Compute t-test for infected vs. uninfected
purp_Severe_t = t.test(y = qc_design_Severe[qc_design_Severe$Infection == "Uninfected", "purple"], x=qc_design_Severe[qc_design_Severe$Infection == "Severe", "purple"])
purp_CoV_t = t.test(y = qc_design_Severe[qc_design_Severe$Infection == "Uninfected", "purple"], x=qc_design_Severe[qc_design_Severe$Infection == "CoV", "purple"])

# Plot boxplot, with t-test values
pdf("./IMAGES/Fig6B_purple_eigen_byInfection_Severe.pdf", width=180/25.4/2, height=180/25.4/2)
ggboxplot(qc_design_Severe, x="Infection", y="purple", add="dotplot", 
          add.params=list(binwidth=diff(range(qc_design_Severe[,"purple"]))/100, alpha=0.5), 
          outlier.shape=NA, main="Cytotoxic immune response", ylab="Eigengene") +
    annotate("text", x=1, y=max(qc_design_Severe[,"purple"])*1.1, size = 3,
             label="dEg:\np-value:") +
    annotate("text", x=3, y=max(qc_design_Severe[,"purple"])*1.1, size = 3,
             label=sprintf("%0.3f\n%0.1e", diff(range(purp_CoV_t$estimate)), purp_CoV_t$p.value)) +
    annotate("text", x=2, y=max(qc_design_Severe[,"purple"])*1.1, size = 3,
             label=sprintf("%0.3f\n%0.1e", diff(range(purp_Severe_t$estimate)), purp_Severe_t$p.value)) +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

##### Panels E and F
qc_design_Severe$ACE2    = unlist(vst["ACE2",qc_design_Severe$SubjectID])
qc_design_Severe$IL6     = unlist(vst["IL6",qc_design_Severe$SubjectID])
panel_genes = list(E="ACE2", F="IL6")
for (panel in names(panel_genes)) {
    gene = panel_genes[[panel]]
    # Plot boxplot, with t-test values
    pdf(sprintf("./IMAGES/Fig6%s_%s_byInfection_Severe.pdf", panel, gene), width=180/25.4/2, height=180/25.4/2)
    g = ggboxplot(qc_design_Severe, x="Infection", y=gene, add="dotplot", 
              add.params=list(binwidth=diff(range(qc_design_Severe[,gene]))/100, alpha=0.5), 
              outlier.shape=NA, main=gene, ylab="VST Normalized Expression") +
        annotate("text", x=1, y=max(qc_design_Severe[,gene])*1.1, size = 3,
                 label="log2FC:\np-value:") +
        annotate("text", x=3, y=max(qc_design_Severe[,gene])*1.1, size = 3,
                 label=sprintf("%0.3f\n%0.1e", common_res_Severe[gene,'logFC_CoV'], common_res_Severe[gene,'adj.P.Val_CoV'])) +
        annotate("text", x=2, y=max(qc_design_Severe[,gene])*1.1, size = 3,
                 label=sprintf("%0.3f\n%0.1e", common_res_Severe[gene,'logFC_Severe'], common_res_Severe[gene,'adj.P.Val_Severe'])) +
        theme(plot.title = element_text(hjust = 0.5))
    print(g)
    dev.off()
}

write.table(common_res_Severe[rownames(gene_sets_Severe[gene_sets_Severe$CoV == 1 & gene_sets_Severe$Severe == 1,]),], file="./OUTPUT/COV_Severe_Combined_DEGs.Common.noSample.txt", sep='\t', quote=F, col.names=NA)
write.table(common_res_Severe[rownames(gene_sets_Severe[gene_sets_Severe$CoV == 1 & gene_sets_Severe$Severe == 0,]),], file="./OUTPUT/COV_Severe_Combined_DEGs.CoV_only.noSample.txt", sep='\t', quote=F, col.names=NA)
write.table(common_res_Severe[rownames(gene_sets_Severe[gene_sets_Severe$CoV == 0 & gene_sets_Severe$Severe == 1,]),], file="./OUTPUT/COV_Severe_Combined_DEGs.Severe_only.noSample.txt", sep='\t', quote=F, col.names=NA)

write.table(common_res_Severe[rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1,]),], file="./OUTPUT/COV_Severe_Combined_DEGs.CoV_onlyUp.noSample.txt", sep='\t', quote=F, col.names=NA)
write.table(common_res_Severe[rownames(gene_sets_Severe[gene_sets_Severe$Severe_up == 1,]),], file="./OUTPUT/COV_Severe_Combined_DEGs.Severe_onlyUp.noSample.txt", sep='\t', quote=F, col.names=NA)

######################################################################
######################################################################

# Find enrichment terms for the shared genes
## Function to limit the output to 8 example genes
condense_genes = function(x){
    y = strsplit(x, ";")[[1]]
    n = length(y)
    out = c()
    for (i in 1:n){
        if (i==8) {
            out = c(out, y[i])
            break
        }
        else if (i==n) sep = ""
        else sep = ", "
        out = c(out, y[i], sep)
    }
    paste(out, collapse="")
}


## Use enrichR to get gene enrichments
### From 2018 GO terms
databases = c("GO_Biological_Process_2018")

### Separate shared genes by up- or down-regulated
both_up_enrich_Severe = enrichr(genes=rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1 & gene_sets_Severe$Severe_up == 1,]), databases=databases)[["GO_Biological_Process_2018"]]
both_dn_enrich_Severe = enrichr(genes=rownames(gene_sets_Severe[gene_sets_Severe$CoV_down == 1 & gene_sets_Severe$Severe_down == 1,]), databases=databases)[["GO_Biological_Process_2018"]]

### Parse for diverse terms, instead of three of the same thing
bothup_top_enrichr_Severe = both_up_enrich_Severe[order(both_up_enrich_Severe$Adjusted.P.value), c("Term", "Adjusted.P.value", "Genes")][1:20,]
bothdn_top_enrichr_Severe = both_dn_enrich_Severe[order(both_dn_enrich_Severe$Adjusted.P.value), c("Term", "Adjusted.P.value", "Genes")][1:20,]
### Reformat for the table
bothup_top_enrichr_Severe[,"Example Genes"] = unlist(sapply(bothup_top_enrichr_Severe$Genes, condense_genes))
bothdn_top_enrichr_Severe[,"Example Genes"] = unlist(sapply(bothdn_top_enrichr_Severe$Genes, condense_genes))


write.table(bothup_top_enrichr_Severe[, c("Term", "Adjusted.P.value", "Example Genes", "Genes")], 
            file="./OUTPUT/bothup_top_enrichr_Severe.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(bothdn_top_enrichr_Severe[, c("Term", "Adjusted.P.value", "Example Genes", "Genes")], 
            file="./OUTPUT/bothdn_top_enrichr_Severe.txt", sep="\t", quote=F, col.names=T, row.names=F)

#=================================================#
# Figure 6 panel C, logFC of CoV vs. logFC of ORV #
#=================================================#

logFC_thresh = 0.5
FDR_thresh = 0.05

common_res_Severe = common_res_Severe[order(apply(common_res_Severe[,c('adj.P.Val_Severe','adj.P.Val_CoV')],1,mean), decreasing=T),]

cov_sig_bool = common_res_Severe$adj.P.Val_CoV < FDR_thresh & abs(common_res_Severe$logFC_CoV) > logFC_thresh
severe_sig_bool = common_res_Severe$adj.P.Val_Severe < FDR_thresh & abs(common_res_Severe$logFC_Severe) > logFC_thresh
sig_degs = cov_sig_bool | severe_sig_bool
significance = ifelse(cov_sig_bool & severe_sig_bool, 'Shared',
                 ifelse(severe_sig_bool, 'ORV',
                 ifelse(cov_sig_bool, 'CoV', 'Neither')))
common_res_Severe$Sig_status = significance

pdf(sprintf("./IMAGES/Fig6c_Severe_logFC_comp.pdf", FDR_thresh, logFC_thresh))
r = cor.test(common_res_Severe$logFC_CoV[sig_degs], common_res_Severe$logFC_Severe[sig_degs])
g = ggscatter(common_res_Severe, x='logFC_Severe', y='logFC_CoV', 
              color="Sig_status", pal=c(CoV_color, neither_color, OTHER_color, both_color),
              main=sprintf("R = %0.2f (p = %0.2e)", r$estimate, r$p.value)
)
print(g)
dev.off()

##### Use METASOFT to compare the logFC values #####

common_genes = rownames(common_res_Severe[sig_degs,])
metasoft_input = data.frame(
name = common_genes,
CoV_logFC = res_fit_CoV[common_genes, "logFC"],
CoV_SE = (sqrt(fit_CoV$s2.post) * fit_CoV$stdev.unscaled)[common_genes, "InfectionCoronavirus"],
Severe_logFC = res_fit_Severe[common_genes, "logFC"],
Severe_SE = (sqrt(fit_Severe$s2.post) * fit_Severe$stdev.unscaled)[common_genes, "InfectionSevereVirus"],
stringsAsFactors=F
)

write.table(metasoft_input, file=sprintf("./OUTPUT/Severe_v_CoV.MetasoftInput.%f.%f.txt", FDR_thresh, logFC_thresh), sep='\t', quote=F, col.names=F, row.names=F)

system(sprintf("
  module load metasoft
  java -jar $METASOFT -input ./OUTPUT/Severe_v_CoV.MetasoftInput.%1$f.%2$f.txt -pvalue_table /software/cgeh/metasoft/default/install/HanEskinPvalueTable.txt -output ./OUTPUT/Severe_v_CoV.MetasoftOutput.%1$f.%2$f.txt -log ./OUTPUT/meta.log
  ", FDR_thresh, logFC_thresh))

meta_res = read.table(sprintf("./OUTPUT/Severe_v_CoV.MetasoftOutput.%f.%f.txt", FDR_thresh, logFC_thresh), skip=1, strings=F)
colnames(meta_res) = c("GENE","NUM_STUDY","PVALUE_FE","BETA_FE","STD_FE","PVALUE_RE","BETA_RE","STD_RE","PVALUE_RE2","STAT1_RE2","STAT2_RE2","PVALUE_BE","I_SQUARE","Q","PVALUE_Q","TAU_SQUARE","STUDY1_P","STUDY2_P","STUDY1_M","STUDY2_M")
meta_res_rmNA = meta_res[!is.na(meta_res$PVALUE_Q),]
meta_res_rmNA$adjPVALUE_Q = p.adjust(meta_res_rmNA$PVALUE_Q, method="BH")
meta_res_rmNA = meta_res_rmNA[order(meta_res_rmNA$PVALUE_Q),]
write.table(meta_res_rmNA, file=sprintf("./OUTPUT/Severe_v_CoV.MetasoftOutput.%f.%f.txt", FDR_thresh, logFC_thresh), sep='\t', quote=F, col.names=T, row.names=F)

#================#
# Data Summaries #
#================#

##### Hypergeometric test for enrichment of DEGs in tan and purple modules #####
# Genes in the relevant WGCNA networks
interferon_genes = rownames(wgcna_modules[wgcna_modules$module == "tan",])
cytotoxic_genes = rownames(wgcna_modules[wgcna_modules$module == "purple",])

# Upregulated genes in each group
CoV_only_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1 & gene_sets_Severe$Severe_up == 0,])
Severe_only_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 0 & gene_sets_Severe$Severe_up == 1,])
both_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1 & gene_sets_Severe$Severe_up == 1,])
genesets_Severe = list(CoV=CoV_only_up, Severe=Severe_only_up, Shared=both_up)

# Hypergeometric test
hypergeom_df_Severe = data.frame(
    Module=character(0),
    Infection=character(0),
    p.value=numeric(0),
    stringsAsFactors=FALSE)
N = nrow(wgcna_modules[wgcna_modules$module != "grey" & !is.na(wgcna_modules$cor),])

m = length(cytotoxic_genes)
n = N-m
for (setName in names(genesets_Severe)) {
    geneset = genesets_Severe[[setName]]
    k = length(geneset)
    x = sum(geneset %in% cytotoxic_genes)
    p.value = phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    hypergeom_df_Severe[nrow(hypergeom_df_Severe)+1,] = c("Purple", setName, p.value)
}

m = length(interferon_genes)
n = N-m
for (setName in names(genesets_Severe)) {
    geneset = genesets_Severe[[setName]]
    k = length(geneset)
    x = sum(geneset %in% interferon_genes)
    p.value = phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    hypergeom_df_Severe[nrow(hypergeom_df_Severe)+1,] = c("Tan", setName, p.value)
}

write.table(hypergeom_df_Severe, "./OUTPUT/CoV_Severe_DEG_WGCNA_enrichment.txt", sep="\t", quote=FALSE, row.names=FALSE)

#================================#
# Save session data to use later #
#================================#

save.image("./OUTPUT/CoV_DEG.RData")
