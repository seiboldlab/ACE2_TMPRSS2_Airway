library("WGCNA")
library(heatmap3)

# load expression data
raw_counts <- read.table("695_raw_counts.txt", header=T, sep="\t", row.names=1)
raw_counts <- as.matrix(raw_counts)
vstMat <- read.table("695_expr_vst.txt", header=T, sep='\t', row.names=1)
vstMat <- as.matrix(vstMat)

# remove lowly expressed genes
# 17473 genes remained
good_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) * .15)

vstMat_rm <- vstMat[good_genes,]

#############
### WGCNA ###
#############

####check genes
expMat <- t(vstMat_rm)
gsg <- goodSamplesGenes(expMat, verbose=3)

if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(colnames(expMat)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(expMat)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    expMat = expMat[gsg$goodSamples, gsg$goodGenes]
}

#Now cluster donors (in contrast to clustering genes later on...)
sampleTree = hclust(dist(expMat), method = "average")

# Plot the sample tree; tips will be ordered to some extent by library size,
# Look for outliers along this continuum
pdf("695.SampleTree_outliers.wgcna.pdf", height=48, width=64)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

enableWGCNAThreads()

#First, choose a set of candidate powers to look at
powers = c(1:30)

# Call the network topology analysis function
sft = pickSoftThreshold(expMat, powerVector = powers, verbose = 5, networkType="signed")

# Plot the results:
pdf("695.pickSoftThresholdPlots.wgcna.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# pick soft threshold
softPower <- sft$powerEstimate + 1

####step by step WGCNA

#1. Create Similarity Matrix
pearson <- WGCNA::cor(as.matrix(expMat),method="pearson")

#2. Convert Similarity Matrix to Adjacency Matrix using Adjacency Function
adjacency.p <- adjacency.fromSimilarity(pearson,type = "signed",power=softPower)

#3. Convert Adjacency to TOM dissimilarity
TOMdissim.p <- 1 - TOMsimilarity(adjacency.p,TOMType = "signed",TOMDenom = "min")

#4. Perform hierarchical clustering on the dissimilarity matrix
geneTree = hclust(as.dist(TOMdissim.p), method = "average")

#### ds=0 -> x=0.64
#### ds=1 -> x=0.73
#### ds=2 -> x=0.82
#### ds=3 -> x=0.91
deepsplits <- c(0.64, 0.73, 0.82, 0.91)

modules <- vector("list")
modcolors <- vector("list")
MEs <- vector("list")
datKME <- vector("list")
gene2module_with_cor <- vector("list")
hub_genes <- vector("list")

#5. cut tree based on hierarchical clustering
for(i in 1:length(deepsplits)) {
    modules[[i]] = cutreeDynamic(dendro = geneTree,method='hybrid', distM = TOMdissim.p,pamStage=T, pamRespectsDendro = T,
        maxCoreScatter=min(geneTree$height)+deepsplits[i]*(max(geneTree$height)-min(geneTree$height)),
        minGap=(1-deepsplits[i])*(max(geneTree$height)-min(geneTree$height)),
        cutHeight=quantile(geneTree$height,.99), minClusterSize=20)
    modcolors[[i]] = labels2colors(modules[[i]])
}

#table(modcolors)

all.modcolors <- do.call( "cbind", modcolors)
colnames(all.modcolors) <- paste0("ds",0:(length(deepsplits)-1))

# dendogram
pdf("695.dendrogram_signed_wgcna.pdf",width=14,height=7)
#plotDendroAndColors(dendro=geneTree,colors=cbind(modcolors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,main=' ')
plotDendroAndColors(dendro=geneTree, colors=all.modcolors, dendroLabels = FALSE,main=' ')
dev.off()

##########################################
####### compute MEs, kME, hubgenes #######
##########################################

for(i in 1:length(deepsplits)) {
    MEs[[i]] <- moduleEigengenes(expMat, modcolors[[i]])$eigengenes
    rownames(MEs[[i]]) <- rownames(expMat)

    write.table(MEs[[i]],  file=paste0("695.WGCNA.ME",i,".txt"), sep='\t', quote=F)

    # create table for genes and their corresponding module
    #gene2module[[i]] = data.frame(gene=colnames(expMat), module=modcolors)
    #gene2module <- na.omit(gene2module)
    #write.table(gene2module, file="./TABLES/695.WGCNA.gene2module.txt",sep="\t",quote=F,row.names=F)

    datKME[[i]] <- signedKME(expMat, MEs[[i]])
    write.table(datKME[[i]], file=paste0("695.WGCNA.KME", i, ".txt"),sep="\t",quote=F)

    gene2module_with_cor[[i]] <- data.frame(gene=colnames(expMat), module=modcolors[[i]])
    gene2module_with_cor[[i]] <- na.omit(gene2module_with_cor[[i]])
    gene2module_with_cor[[i]]$cor <- NA

    for(j in unique(gene2module_with_cor[[i]]$module)) {
        kME_name <- paste0("kME",j)
        idx <- which(gene2module_with_cor[[i]]$module==j)
        gene.idx <- as.character(gene2module_with_cor[[i]][idx,"gene"])
        #gene.idx <- gsub("-",".",gene.idx)
        gene2module_with_cor[[i]]$cor[idx] <- datKME[[i]][gene.idx,kME_name]
        #print(kME_name)
    }
    write.table(gene2module_with_cor[[i]], file=paste0("695.WGCNA.gene2module.with.cor",i,".txt"),sep="\t",quote=F,row.names=F)

    ###hub genes
    ADJ1 <- abs(cor(expMat, use="p"))^softPower
    Alldegrees1 <- intramodularConnectivity(ADJ1, gene2module_with_cor[[i]]$module)
    hub_genes.df <- data.frame()
    for(j in 1:length(unique(gene2module_with_cor[[i]]$module))){
        Alldegrees1$Module <- gene2module_with_cor[[i]]$module
        tmp <- Alldegrees1[Alldegrees1$Module == unique(gene2module_with_cor[[i]]$module)[j], ]
        hub_genes.df <- rbind(hub_genes.df, head(tmp[order(tmp$kWithin, decreasing=T),], n=nrow(tmp)))
    }
    hub_genes[[i]] <- hub_genes.df
    write.table(hub_genes[[i]], file=paste0("695.WGCNA.hub_genes",i,".txt"), sep="\t",quote=F)
}

###########################
##### module dendogram ####
###########################

moduleTree <- vector("list")
for(i in 1:length(deepsplits)) {
    moduleTree[[i]] <- hclust(dist(t(MEs[[i]])), method = "average")
}

pdf("695.ModuleTree.wgcna.pdf")
for(i in 1:length(deepsplits)) {
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(moduleTree[[i]], main = paste0("Module clustering, ds=", i-1), sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
}
dev.off()

