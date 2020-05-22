#!/usr/bin/env Rscript

"
This script imports data from GEO (Accessions GSE3982 & GSE22886)
for flow sorted immune cell types, batch corrects, extracts the
cell types of interest, and then performs 1-vs-all differential
expression for each of the cell types. The the script
takes genes ranked by 1-vs-all comparisons for each of the
flow sorted immune cell types and performs Gene Set Enrichment 
Analysis (GSEA) for the genes that were specifically upregulated in 
CoV+ individuals, Severe+ individuals, or both.

Based in part on code by Lando Ringel and Nathan Dyjack

INPUT:
- Ranked genes for 1-vs-all comparisons (output from 3-rankCells.R)
- Upregulated differentially expressed genes
  - Shared between CoV+ and Severe+
  - CoV+ enhanced DEGs
  - Severe+ enhanced DEGs

OUTPUT:
- GSEA summary table
- Components required for Figure 6, panels F & G
"

#immune cell differential expression'
options(width=250, stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(hgu133a.db))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(fgsea))

set.seed(8994)

#=========================#
# Define helper functions #
#=========================#

limma_geo_diffExpr<-function(group1, allgroups, gset, all=T, group2=NA){
        # performs differential expression between group1 and other groups
        # args:
        #       group1 - Character string or vector containing positive LFC group ID's (i.e. "CD8plus_T_cells_No_Stim" or c('NK_cells_IL15', 'NK_cells_IL2'))
        #       allgroups - Vector containing each sample's group label (i.e. c('Th2', 'Th2', 'Th2', 'Neutrophil', 'Neutrophil', 'NK_cells_IL15'))
        #                   Must be in same order as samples in gset.
        #       gset - Gene expression data-set downloaded from GEO using getGEO
        #       all - TRUE/FALSE compare group1 vs all other groups in allgroups?
        #       group2 - (all=F) other group(s) to compare group1 samples with
        # returns:
        #       table -  containing differentially expressed genes
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
                  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
                  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) {
                ex[which(ex <= 0)] <- NaN
                exprs(gset) <- log2(ex)
        }
        rm(qx)
        if(length(group1) > 1){
                group1<-paste(group1, collapse="|")
        }
        if(!is.na(group2) & all == F){
                if(length(group2) > 1){
                        group2<-paste(group2, collapse="|")
                }
                gset<-gset[, grepl(paste(c(group1, group2), collapse="|"), allgroups)]
                allgroups<-allgroups[grepl(paste(c(group1, group2), collapse="|"), allgroups)]
        }
        desiigner<-ifelse(grepl(group1, allgroups), 1, 0)
        gset$description<-as.factor(paste0("G", desiigner))
        design<-model.matrix(~description+0, gset)
        colnames(design)<-levels(gset$description)
        fit<-lmFit(gset, design)
        cont.matrix<-makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=dim(fit2)[1])
        return(tT)
}

rm_dup_annot<-function(GeneSet, xx, rank){
        # updates annotation from GEO2R differential expression
        # args:
        #       GeneSet - table containing GEO2R differential expression (probeID, Gene Symbol, LFC, p-value)
        #       xx - list containing mappedkeys
        #       rank - select duplicates by 'pval' or 'lfc'
        GeneSet2<-GeneSet
        cat("total probes: ", nrow(GeneSet2), "\n")
        GeneSet2[,'Gene.symbol2']<-as.character(xx[as.character(GeneSet2$ID)])
        # remove multi-matched probes that do not have any selected hgu133aSYMBOL symbol 
        GeneSet_rm<-GeneSet2[!(grepl("///",GeneSet2$Gene.symbol) & GeneSet2$Gene.symbol2 == 'NULL'),]
        # remove all non-annotated probes that do not have any matched hgu133aSYMBOL symbol 
        GeneSet_rm<-GeneSet_rm[!(GeneSet_rm$Gene.symbol == "" & GeneSet_rm$Gene.symbol2 == "NULL"),]
        # GeneSet_rm[] <- lapply(GeneSet_rm, as.character)
        GeneSet_rm$Gene.symbol<-as.character(GeneSet_rm$Gene.symbol)
        # replace all non-annotated probes with matched hgu133aSYMBOL symbol
        GeneSet_rm$Gene.symbol[GeneSet_rm$Gene.symbol == ""] <- GeneSet_rm[GeneSet_rm$Gene.symbol =="",'Gene.symbol2']
        # replace all multi-matched probes with matched hgu133aSYMBOL symbol
        GeneSet_rm$Gene.symbol[grepl("///",GeneSet_rm$Gene.symbol)]<-as.character(GeneSet_rm[grepl("///",GeneSet_rm$Gene.symbol),'Gene.symbol2'])
        # replace all NULL hgu133aSYMBOL matches with original Gene.Symbol annotation
        GeneSet_rm$Gene.symbol[GeneSet_rm$Gene.symbol != GeneSet_rm$Gene.symbol2 & GeneSet_rm$Gene.symbol2 != 'NULL']<-GeneSet_rm$Gene.symbol2[GeneSet_rm$Gene.symbol != GeneSet_rm$Gene.symbol2 & GeneSet_rm$Gene.symbol2 != 'NULL']
        # get duplicated genes 
        duplicated_genes<-GeneSet_rm[duplicated(GeneSet_rm$Gene.symbol), 'Gene.symbol']
        temp_dup<-GeneSet_rm[GeneSet_rm$Gene.symbol %in% duplicated_genes, ]
        temp_dup<-temp_dup[order(temp_dup$Gene.symbol),]
        # merge(temp_dup, aggregate(logFC ~ Gene.symbol, FUN=function(logFC){return(max(abs(logFC)))}, data=temp_dup), by=c('logFC', 'Gene.symbol'), all.y=T)
        # get gene with largest (or smallest) LogFC 
        if(rank == 'lfc'){
                temp_dup2<-do.call(rbind,lapply(split(temp_dup, temp_dup$Gene.symbol),function(genes) genes[which.max(abs(as.numeric(genes$logFC))),]))
                # get genes without duplicates
                GeneSet_rm2<-GeneSet_rm[!GeneSet_rm$Gene.symbol %in% duplicated_genes,]
                GeneSet_noDup<-rbind(temp_dup2,GeneSet_rm2)
                GeneSet_noDup<-GeneSet_noDup[order(as.numeric(GeneSet_noDup[, 'logFC']), decreasing=T), ]
                cat("total probes: ", nrow(GeneSet_noDup), "\n")
                return(GeneSet_noDup)
        }
        else if(rank == 'pval'){
                temp_dup2<-do.call(rbind,lapply(split(temp_dup, temp_dup$Gene.symbol),function(genes) genes[which.min(abs(as.numeric(genes$adj.P.Val))),]))
                # get genes without duplicates
                GeneSet_rm2<-GeneSet_rm[!GeneSet_rm$Gene.symbol %in% duplicated_genes,]
                GeneSet_noDup<-rbind(temp_dup2,GeneSet_rm2)
                GeneSet_noDup<-GeneSet_noDup[order(as.numeric(GeneSet_noDup[, 'adj.P.Val']), decreasing=F), ]
                cat("total probes: ", nrow(GeneSet_noDup), "\n")
                return(GeneSet_noDup)
        }
}

#====================================#
# IMPORT AND PREPARE EXPRESSION DATA #
#====================================#

# Download data from GEO
## Garvan immune cells
gset_GSE3982<-getGEO("GSE3982", GSEMatrix =TRUE, AnnotGPL=TRUE)
gset_GSE3982<-gset_GSE3982[[1]]
gset_GSE3982.2<-gset_GSE3982

## iris immune cells
gset_GSE22886<-getGEO("GSE22886", GSEMatrix =TRUE, AnnotGPL=TRUE)
gset_GSE22886 <- gset_GSE22886[[1]]
# remove corrupted (stimulated) cell
gset_GSE22886.2<-gset_GSE22886[, !grepl("PlasmaCell|NKcell-IL|MemoryTcell|Bcell-Memory", pData(gset_GSE22886)[,'title'])]

## Combine the two datasets
combine_immune<-BiocGenerics::combine(gset_GSE3982.2, gset_GSE22886.2)

## get gene annotation data
geneData<-fData(combine_immune)
phenoData<-pData(combine_immune)
phenoData[,'project']<-c(rep("garvan", ncol(gset_GSE3982.2)), rep("iris", ncol(gset_GSE22886.2)))

# log2 transform the expression data
ex <- exprs(combine_immune)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    exprs(combine_immune) <- log2(ex) 
}
rm(qx, ex)

# Perform batch correction
not_garvan<-data.frame(sample_set=ifelse(phenoData[,'project'] == "garvan", "garvan", "other"))
modcombat = model.matrix(~1, data=not_garvan)
exprs(combine_immune)<-ComBat(dat=exprs(combine_immune), batch=not_garvan[,'sample_set'], mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# Get symbol map for chip HGU133A
x<-hgu133aSYMBOL
mapped_probes<-mappedkeys(x)
xx<-as.list(x[mapped_probes])

#==================================#
# 1-vs-all differential expression #
#==================================#

# Identifiers for each cell type
grep_codes = list(
        NK.Cell="NK",
        CD8T.Cell="CD8Tcell",
        CD4T.Cell = "CD4Tcell",
        Dendritic.Cell="endritic",
        B.Cell="Bcell|B cell",
        Neutrophil="Neutrophil",
        Basophil.Mast.Cell = "Basophil|mast",
        Eosinophils="Eosinophils",
        Macrophage="Macrophage",
        Monocyte.Day0="Monocyte-Day0",
        Monocyte.Day1="Monocyte-Day1"
    )

# Run limma on each comparison
for (gcode in names(grep_codes)){
    celltypes_replace = ifelse(grepl(grep_codes[[gcode]], phenoData$title), gcode, "other")
    cell_results = limma_geo_diffExpr(group1=gcode, allgroups=celltypes_replace, gset=combine_immune, all=T)
    cell_results = rm_dup_annot(GeneSet=cell_results, xx=xx, rank='pval')
    cell_rnk = cell_results[order(cell_results$logFC,decreasing=T),]
    write.table(cell_rnk[,c('Gene.symbol','logFC')],sprintf("./OUTPUT/%s_lfc_full.rnk",gcode),row.names=F,quote=F,sep="\t",col.names=F)
}

#===========#
# LOAD DATA #
#===========#

shared_degs = read.table("./OUTPUT/COV_Severe_Combined_DEGs.Common.noSample.txt")
cov_degs = read.table("./OUTPUT/COV_Severe_Combined_DEGs.CoV_only.noSample.txt")
hrv_degs = read.table("./OUTPUT/COV_Severe_Combined_DEGs.Severe_only.noSample.txt")

ranked_gene_dir = "./OUTPUT/"

#===========#
# PREP DATA #
#===========#

shared_up   = rownames(shared_degs[shared_degs$logFC_CoV > 0,])
cov_up   = rownames(cov_degs[cov_degs$logFC_CoV > 0,])
hrv_up   = rownames(hrv_degs[hrv_degs$logFC_Severe > 0,])

pathways = list(
                Shared=as.list(shared_up), 
                CoV_specific=as.list(cov_up), 
                Severe_specific=as.list(hrv_up))

celltypes = c(
        "NK.Cell",
        "CD8T.Cell",
        "CD4T.Cell",
        "Dendritic.Cell",
        "B.Cell",
        "Neutrophil",
        "Basophil.Mast.Cell",
        "Eosinophils",
        "Macrophage",
        "Monocyte.Day0",
        "Monocyte.Day1"
    )

#=====================================#
# RUN GSEA AND GENERATE SUMMARY TABLE #
#=====================================#

# Initialize empty data frame
fgsea_summary = data.frame(
    celltype=character(0),
    CoV_ES=numeric(0),
    CoV_p=numeric(0),
    Severe_ES=numeric(0),
    Severe_p=numeric(0),
    Shared_ES=numeric(0),
    Shared_p=numeric(0),
    stringsAsFactors=FALSE)

for (celltype in celltypes){
    # Read in ranked genes and convert to named numeric vector
    ranked_genes = read.table(sprintf("%s%s_lfc_full.rnk", ranked_gene_dir, celltype),sep="\t",header=F,row.names=1)
    stats = ranked_genes[,1]
    names(stats) = rownames(ranked_genes)

    # Run FGSEA on these genes, vs. the DEGs
    fgseaRes = fgseaMultilevel(pathways=pathways, stats=stats, scoreType='pos', eps=0)
    
    # Use the plot function to extract the enrichment scores
    g_shared = plotEnrichment(pathways[[1]], stats)
    g_cov = plotEnrichment(pathways[[2]], stats)
    g_hrv = plotEnrichment(pathways[[3]], stats)

    # Pull out the enrichment score and p-value for each
    cov_es = fgseaRes[[1,"ES"]]
    hrv_es = fgseaRes[[2,"ES"]]
    shared_es = fgseaRes[[3,"ES"]]
    cov_p = fgseaRes[[1,"padj"]]
    hrv_p = fgseaRes[[2,"padj"]]
    shared_p = fgseaRes[[3,"padj"]]

    result = c(celltype, cov_es, cov_p, hrv_es, hrv_p, shared_es, shared_p)
    fgsea_summary[nrow(fgsea_summary)+1,] = result
}

# Write out summary table
write.table(fgsea_summary, "./OUTPUT/celltype_enrichment_summary_Severe.txt", sep="\t", quote=F, row.names=F)

#======================================#
# Generate plots for CD8T and NK cells #
#======================================#

##### Plotting helpers #####
# Color palette
logFC_colors = colorRampPalette(c("#377eb8", "white", "#e41a1c"))(1000)

#######################
##### CD8+ T-CELL #####
#######################
celltype = "CD8T.Cell"

# Read in logFC genes
ranked_genes = read.table(sprintf("%s%s_lfc_full.rnk", ranked_gene_dir, celltype),sep="\t",header=F,row.names=1)
stats = ranked_genes[,1]
names(stats) = rownames(ranked_genes)

# Identify where the the logFC goes from (+) to (-)
zeropoint = NA
for (i in 1:(length(stats)-1)){if (stats[i] > 0 & stats[i+1]<=0){zeropoint = i}}

# Repeat the FGSEA
fgseaRes = fgseaMultilevel(pathways=pathways, stats=stats, scoreType='pos', eps=0)

# Extract the leading edge
goi = unlist(fgseaRes[[1,"leadingEdge"]])[1:30][c(T,F,F)]
fgseaRes$leadingEdge = sapply(fgseaRes$leadingEdge, function(x){paste(unlist(x), collapse=",")})
write.table(fgseaRes, sprintf("%s/%s_GSEA_results.final_Severe.txt", ranked_gene_dir, celltype), sep="\t", quote=F, row.names=F)

# Use the plotting function to get access to the curve for custom plotting
g_shared = plotEnrichment(pathways[[1]], stats)
g_cov = plotEnrichment(pathways[[2]], stats)
g_hrv = plotEnrichment(pathways[[3]], stats)

# Plot the logFC values as a gradient colorbar
pdf(sprintf("./IMAGES/%s_logFCs.final_Severe.pdf", celltype), width=180/25.4, height=20/25.4)
    maxlogFC = max(abs(stats))
    breakscale = c(-maxlogFC-1, seq(-1.5,1.5, length.out=length(logFC_colors)-1),maxlogFC+1)
    barcolors = logFC_colors[as.numeric(cut(stats, breaks=breakscale))]
    par(mar = c(0, 4.1, 0, 2.1))
    barplot(rep(1, length(stats)), width=1, space=0, xaxt='n', ann=FALSE, col=barcolors, border=NA, axes=FALSE)
    abline(v=zeropoint, lty=2, col='gray')
dev.off()

# Plot the enrichment score curve, with the following attributes:
# - Representative leading edge genes for CoV+
# - The genes in each set notated as vertical bars below the curve
pdf(sprintf("./IMAGES/%s_GSEA_DEGs.final_Severe.pdf", celltype), width=180/25.4, height=140/25.4)
    par(mar = c(0, 4.1, 2, 2.1))
    # Min and Max enrichment score values, for setting limits
    min_val = min(c(g_cov$data[,2], g_hrv$data[,2], g_shared$data[,2]))
    max_val = max(c(g_cov$data[,2], g_hrv$data[,2], g_shared$data[,2]))

    # Plot the enrichment scores
    plot(g_shared$data[,1], g_shared$data[,2], lty=1, col='black', type="l", ylim=c(min_val-0.15, max_val),
        ylab="Enrichment Score", xlab=sprintf("%s gene rank",celltype), main=NA,
        axes=FALSE)
    lines(g_cov$data[,1], g_cov$data[,2], lty=1, col='darkgoldenrod1')
    lines(g_hrv$data[,1], g_hrv$data[,2], lty=1, col='red')

    # Identify positions for annotating text for leading edge
    gene_indices = which(names(stats) %in% goi)
    gene_indicesb = which(g_cov$data[,1] %in% gene_indices)

    points(g_cov$data[,1][gene_indicesb], g_cov$data[,2][gene_indicesb], lty=1, col='darkorchid4', pch=20)
    text(g_cov$data[,1][gene_indicesb], 
         g_cov$data[,2][gene_indicesb],
         labels=names(stats)[gene_indices],
         pos = 4, 
         adj=c(0,0.5))

    # Add vertical lines at position of each gene in each set below the plot
    for (i in g_shared$data[,1]){lines(c(i,i), c(min_val-0.05,min_val), col='black')}
    for (i in g_cov$data[,1]){lines(c(i,i), c(min_val-0.1,min_val-0.05), col='darkgoldenrod1')}
    for (i in g_hrv$data[,1]){lines(c(i,i), c(min_val-0.15,min_val-0.1), col='red')}

    # Set up the axes
    abline(h=0)
    axis(side = 2, at = seq(-0.4,0.4,0.2), labels = formatC(seq(-0.4,0.4,0.2),digits=1,format="f"),las=2)

    # Add legend and zero point
    legend("topright", c("Shared", "CoV Specific", "Severe Specific"), lty=1, col=c('black', 'darkgoldenrod1', 'red'))
    abline(v=zeropoint, lty=2, col='gray')
dev.off()

###############################
##### NATURAL KILLER CELL #####
###############################
celltype = "NK.Cell"

# Read in logFC genes
ranked_genes = read.table(sprintf("%s%s_lfc_full.rnk", ranked_gene_dir, celltype),sep="\t",header=F,row.names=1)
stats = ranked_genes[,1]
names(stats) = rownames(ranked_genes)

# Identify where the the logFC goes from (+) to (-)
zeropoint = NA
for (i in 1:(length(stats)-1)){if (stats[i] > 0 & stats[i+1]<=0){zeropoint = i}}

# Repeat the FGSEA
fgseaRes = fgseaMultilevel(pathways=pathways, stats=stats, scoreType='pos', eps=0)

# Extract the leading edge
goi = unlist(fgseaRes[[1,"leadingEdge"]])[1:30][c(T,F,F)]
goi2 = unlist(fgseaRes[[3,"leadingEdge"]])[1:30][c(T,F,F)]
fgseaRes$leadingEdge = sapply(fgseaRes$leadingEdge, function(x){paste(unlist(x), collapse=",")})
write.table(fgseaRes, sprintf("%s/%s_GSEA_results.final_Severe.txt", ranked_gene_dir, celltype), sep="\t", quote=F, row.names=F)

# Use the plotting function to get access to the curve for custom plotting
g_shared = plotEnrichment(pathways[[1]], stats)
g_cov = plotEnrichment(pathways[[2]], stats)
g_hrv = plotEnrichment(pathways[[3]], stats)

# Plot the logFC values as a gradient colorbar
pdf(sprintf("./IMAGES/%s_logFCs.final_Severe.pdf", celltype), width=180/25.4, height=20/25.4)
    maxlogFC = max(abs(stats))
    breakscale = c(-maxlogFC-1, seq(-1.5,1.5, length.out=length(logFC_colors)-1),maxlogFC+1)
    barcolors = logFC_colors[as.numeric(cut(stats, breaks=breakscale))]
    par(mar = c(0, 4.1, 0, 2.1))
    barplot(rep(1, length(stats)), width=1, space=0, xaxt='n', ann=FALSE, col=barcolors, border=NA, axes=FALSE)
    abline(v=zeropoint, lty=2, col='gray')
dev.off()

# Plot the enrichment score curve, with the following attributes:
# - Representative leading edge genes for CoV+ and Shared genes
# - The genes in each set notated as vertical bars below the curve
pdf(sprintf("./IMAGES/%s_GSEA_DEGs.final_Severe.pdf", celltype), width=180/25.4, height=140/25.4)
    par(mar = c(0, 4.1, 2, 2.1))
    # Min and Max enrichment score values, for setting limits
    min_val = min(c(g_cov$data[,2], g_hrv$data[,2], g_shared$data[,2]))
    max_val = max(c(g_cov$data[,2], g_hrv$data[,2], g_shared$data[,2]))

    # Plot the enrichment scores
    plot(g_shared$data[,1], g_shared$data[,2], lty=1, col='black', type="l", ylim=c(min_val-0.15, max_val),
        ylab="Enrichment Score", xlab=sprintf("%s gene rank",celltype), main=NA,
        axes=FALSE)
    lines(g_cov$data[,1], g_cov$data[,2], lty=1, col='darkgoldenrod1')
    lines(g_hrv$data[,1], g_hrv$data[,2], lty=1, col='red')

    # Identify positions for annotating text for leading edge
    gene_indices = which(names(stats) %in% goi)
    gene_indicesb = which(g_cov$data[,1] %in% gene_indices)
    gene_indices2 = which(names(stats) %in% goi2)
    gene_indices2b = which(g_shared$data[,1] %in% gene_indices2)

    points(g_cov$data[,1][gene_indicesb], g_cov$data[,2][gene_indicesb], lty=1, col='darkorchid4', pch=20)
    points(g_shared$data[,1][gene_indices2b], g_shared$data[,2][gene_indices2b], lty=1, col='gray30', pch=20)
    text(g_cov$data[,1][gene_indicesb], 
         g_cov$data[,2][gene_indicesb],
         labels=names(stats)[gene_indices],
         pos = 4, 
         adj=c(0,0.5))
    text(g_shared$data[,1][gene_indices2b], 
         g_shared$data[,2][gene_indices2b],
         labels=names(stats)[gene_indices2],
         pos = 2, 
         adj=c(1,0.5))

    # Add vertical lines at position of each gene in each set below the plot
    for (i in g_shared$data[,1]){lines(c(i,i), c(min_val-0.05,min_val), col='black')}
    for (i in g_cov$data[,1]){lines(c(i,i), c(min_val-0.1,min_val-0.05), col='darkgoldenrod1')}
    for (i in g_hrv$data[,1]){lines(c(i,i), c(min_val-0.15,min_val-0.1), col='red')}

    # Set up the axes
    abline(h=0)
    axis(side = 2, at = seq(-0.4,0.4,0.2), labels = formatC(seq(-0.4,0.4,0.2),digits=1,format="f"),las=2)

    # Add legend and zero point line
    legend("topright", c("Shared", "CoV Specific", "Severe Specific"), lty=1, col=c('black', 'darkgoldenrod1', 'red'))
    abline(v=zeropoint, lty=2, col='gray')
dev.off()
