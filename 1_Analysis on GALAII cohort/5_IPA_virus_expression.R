#!/usr/bin/env Rscript

desc = "
This script takes the output from the IPA canonical pathway
analysis, and the differential expression, and generates a
heatmap for Figure 6e.

INPUTS:
- IPA Canonical Pathway gene overlaps with upregulated genes in 
  CoV infected individuals for the following pathways:
  - CTLA4 Signaling in Cytotoxic T Lymphocytes
  - Cytotoxic T Lymphocyte-meditaed Apoptosis of Target Cells
  - Natural Killer Cell Signaling
  - OX40 Signaling Pathway
  - iCOS-iCOSL Signaling in T Helper Cells

OUTPUTS:
- Panels for Figure 6:
  - Fig6E: Heatmap of logFC values for top upregulated genes, 
           and genes from the pathways above.
           virus_logFC_heatmap.pdf
"

suppressPackageStartupMessages(library(RSkittleBrewer))
suppressPackageStartupMessages(library(pheatmap))
args = commandArgs(trailingOnly=TRUE)

#=======================================#
# Load saved objects from previous step #
#=======================================#

load("./OUTPUT/CoV_DEG.RData")

#===========================#
# Prepare data for plotting #
#===========================#

# Estabish categories of genes from the venn diagram
cov_specific_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1 & gene_sets_Severe$Severe == 0,])
Severe_specific_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV == 0 & gene_sets_Severe$Severe_up == 1,])
shared_up = rownames(gene_sets_Severe[gene_sets_Severe$CoV_up == 1 & gene_sets_Severe$Severe_up == 1,])

# Filter the genes down to some examples with highest fold changes
top_shared = shared_up[order(common_res_Severe[shared_up, "logFC_CoV"], decreasing=T)][1:10]
top_cov = cov_specific_up[order(common_res_Severe[cov_specific_up, "logFC_CoV"], decreasing=T)][1:10]
top_Severe = Severe_specific_up[order(common_res_Severe[Severe_specific_up, "logFC_Severe"], decreasing=T)][1:10]

# Pull out the genes from the IPA Canonical Pathway analysis
ipa_path = args[[1]]#"./FIXED3_IPA_Pathway_Enrichment/"
pathway_files = c(
    "CTLA4 Signaling in Cytotoxic T Lymphocytes",
    "Cytotoxic T Lymphocyte-mediated Apoptosis of Target Cells",
    "Natural Killer Cell Signaling",
    "OX40 Signaling Pathway",
    "iCOS-iCOSL Signaling in T Helper Cells"
  )

pathway_genes = list()
for (path in pathway_files){
    df = read.table(sprintf("%s/%s.txt", ipa_path, path), sep="\t", skip=2, header=TRUE, stringsAsFactors=FALSE)
    # These are genes that were upregulated in CoV+ that overlapped this pathway
    genes = df[,"Symbol"]
    pathway_genes[[path]] = unlist(genes)
}
## Combine genes from all pathways, and sort, since there is overlap between the pathways.
enriched_genes = unique(unlist(pathway_genes))
enriched_genes = enriched_genes[order(
  apply(common_res_Severe[enriched_genes, c("logFC_Severe", "logFC_CoV")], 1, function(x){diff(x)*x[2]}),
  decreasing=T)]

# Put all of the genes together
genes_of_interest1 = c("ACE2", "IL6")
genes_of_interest = c(top_shared,    top_cov,    top_Severe,
                      genes_of_interest1, enriched_genes)
genes_of_interest = genes_of_interest[genes_of_interest %in% rownames(common_res_Severe)]

#==================#
# Plotting helpers #
#==================#

reds = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(1000)
breakscale = seq(0.0,2.5, length.out=length(reds)+1)
tropical = RSkittleBrewer('tropical')

#===================#
# Generate the plot #
#===================#

# Establish color bars
## By WGCNA module
gene_cols = ifelse(genes_of_interest %in% cytotoxic_genes, "purple",
            ifelse(genes_of_interest %in% interferon_genes, "tan",
            'white'))
## By IPA Pathway
i = length(tropical)
color_guide = list(Network=c(tan='tan', purple='purple', white=NA))
enriched_colors = data.frame(genes=genes_of_interest, stringsAsFactors=FALSE)
for (path in pathway_files){
    enriched_colors[,path] = factor(ifelse(genes_of_interest %in% pathway_genes[[path]], tropical[i], "white"), levels=c('white', tropical[i]))
    color_guide[[path]] = c()
    color_guide[[path]][tropical[i]] = tropical[i]
    color_guide[[path]]['white'] = NA
    i = i - 1
}
dups = duplicated(enriched_colors$genes)
enriched_colors$genes[dups] = paste0(enriched_colors$genes[dups], '.1')
rownames(enriched_colors) = enriched_colors$genes
enriched_colors$genes = gene_cols
colnames(enriched_colors)[1] = "Network"

# as.matrix(enriched_colors)[1:4, 1:4]

plotMat = t(common_res_Severe[genes_of_interest, c("logFC_Severe","logFC_CoV")])
plotMat[plotMat > 2.5] = 2.5
plotMat[plotMat < 0] = 0

dups = duplicated(colnames(plotMat))
colnames(plotMat)[dups] = paste0(colnames(plotMat)[dups], '.1')

# The heatmap
pdf("./IMAGES/virus_logFC_heatmap_Severe.pdf", width=15, height=4)
pheatmap(plotMat,
         color=reds, kmeans_k=NA, breaks=breakscale,
         border_color=NA, cellwidth=NA, cellheight=NA,
         scale='none', cluster_rows=F, cluster_cols=F, labels_row=c("Severe", "CoV"),
         legend=T, legend_breaks=seq(0, 2.5, length.out=6),
         annotation_col=enriched_colors, annotation_names_col = T, annotation_legend=F, 
         annotation_colors = color_guide, fontsize_col=8
    )
dev.off()
