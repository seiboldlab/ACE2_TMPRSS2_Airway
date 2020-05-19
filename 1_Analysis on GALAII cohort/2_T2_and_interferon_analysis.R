library(ggplot2)
library(cowplot)
library(ggpubr)
library(heatmap3)
library(dendextend)
library(ggbeeswarm)

### load normalized expression data
vstMat <- read.table("695_expr_vst.txt", header=T, sep='\t', row.names=1)
vstMat <- as.matrix(vstMat)

### load Module Eigengene values
MEs <- read.table("695.WGCNA.ME3.txt", sep='\t', header=T)

### remove the grey module
MEs_no_grey <- MEs[,-19]

### load gene2module table
gene2mod <- read.table("695.WGCNA.gene2module.with.cor3.txt", header=T)

phen <- data.frame(MEs, t(vstMat[c("ACE2", "TMPRSS2"), rownames(MEs)]))

### compute correlation between module eigengenes
MEs_cor <- cor(MEs_no_grey, MEs_no_grey)

MEs_cor_test <- apply(MEs_no_grey, 2, function(x) apply(MEs_no_grey, 2, function(y) cor.test(x, y)$p.value))
MEs_cor.df <- MEs_cor

for(i in 1:(dim(MEs_cor.df)[2]-1)) {
    for(j in (i+1):(dim(MEs_cor.df)[2])) {
        MEs_cor.df[j,i] <- MEs_cor_test[j,i]
    }
}

# Supplementary table: full correlation matrix of network eigengenes
write.table(MEs_cor.df, file="MEs_cor.txt", sep='\t', quote=F)

### tan: virus
### purple: immune (ACE2)
### pink: mucin (TMPRSS2)
### saddlebrown: core type2

### Figure 2b-f and Figure 4a,b
pdf("ACE2_TMPRSS2_MEs_correlation.pdf", width=3.5, height=3.5)
g <- ggscatter(phen, x = "ACE2", y = "MEtan", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "ACE2", y = "MEpurple", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "ACE2", y = "MEpink", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "ACE2", y = "MEsaddlebrown", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "TMPRSS2", y = "MEtan", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "TMPRSS2", y = "MEpurple", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "TMPRSS2", y = "MEpink", add = "reg.line", conf.int=T) + stat_cor()
print(g)
g <- ggscatter(phen, x = "TMPRSS2", y = "MEsaddlebrown", add = "reg.line", conf.int=T) + stat_cor()
print(g)
dev.off()

#==========================================#
# Perform hierarchical clustering to make: #
# 1. tan : Interferon high and low groups  #
# 2. saddlebrown: T2 high and low groups   #
#==========================================#

# hclust based on tan 
# tan: 296 genes

tan_genes <- as.character(gene2mod[gene2mod$module=="tan", "gene"])

tan_hm <- heatmap3(vstMat[tan_genes,], scale="row", col=mycols, breaks=breakscale, method="ward.D2", keep.dendro=T)

ct_tan <- cutree(tan_hm$Colv, k=4, order_clusters_as_data=F)
tan_WGCNA_dat_mean <- apply(vstMat[tan_genes,], 2, mean)

mean(tan_WGCNA_dat_mean[names(ct_tan[ct_tan==1])])
mean(tan_WGCNA_dat_mean[names(ct_tan[ct_tan==2])])
mean(tan_WGCNA_dat_mean[names(ct_tan[ct_tan==3])])
mean(tan_WGCNA_dat_mean[names(ct_tan[ct_tan==4])])

interferon_low_donors <- names(ct_tan[ct_tan==3 | ct_tan==4])
interferon_high_donors <- names(ct_tan[ct_tan==1 | ct_tan==2])

# hclust based on saddlebrown 
# saddlebrown 156 genes

sb_genes <- as.character(gene2mod[gene2mod$module=="saddlebrown", "gene"])
sb_hm <- heatmap3(vstMat[sb_genes,], scale="row", col=mycols, breaks=breakscale, method="ward.D2", keep.dendro=T)

ct_sb <- cutree(sb_hm$Colv, k=2, order_clusters_as_data=F)
sb_WGCNA_dat_mean <- apply(vstMat[sb_genes,], 2, mean)

mean(sb_WGCNA_dat_mean[names(ct_sb[ct_sb==1])])
mean(sb_WGCNA_dat_mean[names(ct_sb[ct_sb==2])])

type2_low_donors <- names(ct_sb[ct_sb==1])
type2_high_donors <- names(ct_sb[ct_sb==2])

interferon_type2.df <- data.frame(donorID=colnames(vstMat))
rownames(interferon_type2.df) <- interferon_type2.df$donorID
interferon_type2.df$interferon_status <- NA
interferon_type2.df[interferon_low_donors, "interferon_status"] <- "interferon_low"
interferon_type2.df[interferon_high_donors, "interferon_status"] <- "interferon_high"
interferon_type2.df$type2_status <- NA
interferon_type2.df[type2_low_donors, "type2_status"] <- "type2_low"
interferon_type2.df[type2_high_donors, "type2_status"] <- "type2_high"

write.table(interferon_type2.df, file="interferon_type2.status.txt", sep='\t', row.names=F, quote=F)

#### breakscale for heatmap3
mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale <- c(-8,seq(-1.8,1.8, length.out=length(mycols)-1),8)

# Supplementary Figure 1a, 2a
pdf("./IMAGES/interferon_type2_heatmap.pdf", width=12, height=18)
sb_hm <- heatmap3(vstMat[sb_genes,], scale="row", col=mycols, breaks=breakscale, method="ward.D2", keep.dendro=T, ColSideColors=plyr::mapvalues(virus_type2.df$type2_status, from=c("type2_low", "type2_high"), to=c("blue", "red")))
tan_hm <- heatmap3(vstMat[tan_genes,], scale="row", col=mycols, breaks=breakscale, method="ward.D2", keep.dendro=T, ColSideColors=plyr::mapvalues(virus_type2.df$virus_status, from=c("virus_low", "virus_high"), to=c("blue", "red")))
dev.off()

#===================================================#
# ACE2 and TMPRSS2 expression based on type2 status #
#===================================================#

dat <- data.frame(interferon_type2.df, t(vstMat[c("ACE2", "TMPRSS2"), rownames(interferon_type2.df)]))
dat$type2_status <- factor(dat$type2_status, levels=c("type2_low", "type2_high"))
dat$ACE2_norm <- t(norm_dat["ACE2", rownames(interferon_type2.df)])
dat$TMPRSS2_norm <- t(norm_dat["TMPRSS2", rownames(interferon_type2.df)])

# Figure 2d and 2g
pdf("./IMAGES/ACE2_TMPRSS2_by_type2_status.pdf", width=3.5, height=3.5)
g <- ggboxplot(dat, x = "type2_status", y = "ACE2") + stat_compare_means()
print(g)
g <- ggplot(dat, aes(x=type2_status, y=ACE2)) + geom_boxplot(aes(fill=type2_status), outlier.shape=NA)
# + geom_beeswarm(alpha=0.4, size=0.1)
g <- g + scale_fill_manual(values=c("#bebebe","#ff6347")) + theme(legend.position="none")
print(g)
g <- ggboxplot(dat, x = "type2_status", y = "TMPRSS2") + stat_compare_means()
print(g)
g <- ggplot(dat, aes(x=type2_status, y=TMPRSS2)) + geom_boxplot(aes(fill=type2_status), outlier.size=0.1, outlier.alpha=0.4) + geom_beeswarm(alpha=0.4, size=0.1)
g <- g + scale_fill_manual(values=c("#bebebe","#ff6347")) + theme(legend.position="none")
print(g)
dev.off()

ACE2_L2FC <- log2(mean(dat[dat$type2_status=="type2_high", "ACE2_norm"]) / mean(dat[dat$type2_status=="type2_low", "ACE2_norm"]))
TMPRSS2_L2FC <- log2(mean(dat[dat$type2_status=="type2_high", "TMPRSS2_norm"]) / mean(dat[dat$type2_status=="type2_low", "TMPRSS2_norm"]))

#============================================#
# ACE2 expression based on interferon status #
#============================================#

# Peter could you add this section
# Figure 4d