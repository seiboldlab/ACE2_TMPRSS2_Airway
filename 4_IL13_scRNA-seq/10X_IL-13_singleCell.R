library(Seurat)
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("CommonFunctions.r")

# 1. First bring in the four count tables from cell ranger
# Datasets
dat_bsa<-Read10X("/Seibold/proj/Single_Cell_Seq/10X_SPDEF_KO/190228_A00405_0082_BHHW35DSXX/MK_COUNT_TABLE/scrb_BSA/outs/filtered_feature_bc_matrix")
dat_il13<-Read10X("/Seibold/proj/Single_Cell_Seq/10X_SPDEF_KO/190228_A00405_0082_BHHW35DSXX/MK_COUNT_TABLE/scrb_IL13/outs/filtered_feature_bc_matrix")

# 2. Make summary stats table for each of the samples and decide what to keep
# Make list of summary tables
sumTabList<-list()
namesVec<-c("dat_bsa","dat_il13")
for(i in 1:length(namesVec)){
	currData<-as.matrix(get(namesVec[i]))
	sumTabList[[length(sumTabList) + 1]]<-as.data.frame(matrix(NA,nrow=ncol(currData),ncol=4))
	names(sumTabList)[[length(sumTabList)]]<-namesVec[i]
	colnames(sumTabList[[length(sumTabList)]])<-c("ID","nUMI","nGene","propMT")
	sumTabList[[length(sumTabList)]]$ID<-colnames(currData)
	sumTabList[[length(sumTabList)]]$nUMI<-colSums(currData)
	sumTabList[[length(sumTabList)]]$nGene<-colSums(currData != 0)
	mtGenes<-grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^MRPS",rownames(currData),value=T)
	sumTabList[[length(sumTabList)]]$propMT<-colSums(currData[mtGenes,]) / sumTabList[[length(sumTabList)]]$nUMI
	riboGenes<-grep("^RPL|^RPS",rownames(currData),value=T)
	sumTabList[[length(sumTabList)]]$propRibo<-colSums(currData[riboGenes,]) / sumTabList[[length(sumTabList)]]$nUMI
}

### Export these tables to xcel
write.xlsx(sumTabList, file = "10X_SCBL_summaryStats.xlsx")

########################## Histograms of summary stats

#nGenes
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")
#Equal size bins
bins <- function(x){
	c(seq(min(x),max(x),by = 200),max(x))
}
pdf("10X_SCBL_nGene_hist_separate.pdf",width=8,height=4)
#dev.new(width=8,height=4)
par(mfrow=c(1,2))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$nGene,breaks=bins,las=1,main = names(sumTabList)[[i]],xlab="Number of genes",col=colVec[i],ylab="")
	abline(v=quantile(currData$nGene,0.99),col="red")
	abline(v=1500,col="red")
}	
dev.off()

#nUMI
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")
#Equal size bins
bins <- function(x){
	c(seq(min(x),max(x),by = 2000),max(x))
}
pdf("10X_SCBL_nUMI_hist_separate.pdf",width=8,height=4)
#dev.new(width=8,height=4)
par(mfrow=c(1,2))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$nUMI,breaks=bins,las=1,main = names(sumTabList)[[i]],xlab="Number of UMIs",col=colVec[i],ylab="")
	abline(v=quantile(currData$nUMI,0.99),col="red")
#	abline(v=quantile(currData$nUMI,0.01),col="red")
}	
dev.off()

#propMT
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")

pdf("10X_SCBL_propMT_hist_separate.pdf",width=8,height=4)
#dev.new(width=8,height=4)
par(mfrow=c(1,2))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$propMT,breaks=50,las=1,main = names(sumTabList)[[i]],xlab="Proportion MT genes",col=colVec[i],ylab="")
	abline(v=0.3,col="red")
}	
dev.off()

# 3. Filter out cells, merge datasets, and filter out genes

#Identify the samples that are outliers to be removed
removeList<-list()
for(i in 1:length(namesVec)){
removeList[[length(removeList) + 1]]<-sumTabList[[i]]$ID[which(sumTabList[[i]]$nGene > quantile(sumTabList[[i]]$nGene,0.99) | 
	sumTabList[[i]]$nGene < 1500 | sumTabList[[i]]$nUMI > quantile(sumTabList[[i]]$nUMI,0.99) | sumTabList[[i]]$propMT > 0.3)]
names(removeList)[[length(removeList)]]<-namesVec[i]
} 

#Original number of cells 
NumCellsVec<-sapply(sumTabList,function(x)nrow(x))
# dat_bsa dat_il13 
#    5971     5401 

#Number gone (29,661 total tossed)
NumGoneVec<-sapply(removeList,function(x)length(x)) 
# dat_bsa dat_il13 
#    1717     2686 

#Proportion gone
PropGoneVec<-NumGoneVec / NumCellsVec
#  dat_bsa  dat_il13 
#0.2875565 0.4973153 

#Number remaining (38,678 total left)
NumRemainVec<-NumCellsVec - NumGoneVec
# dat_bsa dat_il13 
#    4254     2715 

#Proportion remaining
PropRemainVec<-NumRemainVec / NumCellsVec
#  dat_bsa  dat_il13 
#0.7124435 0.5026847 

##### Filter out samples
dat_bsaf<-as.matrix(dat_bsa)[,-which(colnames(dat_bsa) %in% removeList$dat_bsa)]
dat_il13f<-as.matrix(dat_il13)[,-which(colnames(dat_il13) %in% removeList$dat_il13)]

##### Make summary stat tables based on filtered cells
sumTabList_f<-list()
for(i in 1:length(namesVec)){
	currData<-get(paste(namesVec[i],"f",sep=""))
	sumTabList_f[[length(sumTabList_f) + 1]]<-sumTabList[[i]][which(sumTabList[[i]]$ID %in% colnames(currData)),]
}

############### Combine datasets and remove genes
#Merge samples into 1
SCBLf<-cbind(dat_bsaf,dat_il13f)

#Add an identifier to barcodes, so there will be no duplicates for each donor
colnames(SCBLf)<-c(paste(colnames(dat_bsaf),"1",sep="_"),paste(colnames(dat_il13f),"2",sep="_"))

#Filter out genes (toss 223 genes; 33315 genes remaining
SCBLf<-SCBLf[-c(grep("MTAT|MT-|MTCO|MTCY|MTERF|MTND|MTRF|MTRN|MRPL|MRPS|RPL|RPS",rownames(SCBLf))),]

#Initially, just toss genes with zero expression (toss 9,885; 23,430 remaining)
SCBLf<-SCBLf[-which(rowSums(SCBLf) == 0),]

#Make meta.data table
SCBL_metadata<-data.frame("treatment"=c(rep("BSA",ncol(dat_bsaf)),rep("IL13",ncol(dat_il13f))),
	"nUMI"=c(sumTabList_f[[1]]$nUMI,sumTabList_f[[2]]$nUMI),
	"nGene"=c(sumTabList_f[[1]]$nGene,sumTabList_f[[2]]$nGene),
	"propMT"=c(sumTabList_f[[1]]$propMT,sumTabList_f[[2]]$propMT),
	"propRibo"=c(sumTabList_f[[1]]$propRibo,sumTabList_f[[2]]$propRibo),
	row.names=colnames(SCBLf))

# 4. SEURAT analysis using the full filtered dataset
#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.6.1 && env R_MAX_VSIZE=700Gb R CMD BATCH integrateFull_scTransform_3000Features.r"
#### min cells equal 0.1% of data
min.cells <- floor(0.001 * dim(SCBLf)[2]) 
SCBL <- CreateSeuratObject(counts = SCBLf, min.cells = min.cells, min.features = 0, meta.data = SCBL_metadata, project = "SCBL")

#Integration

#To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.
SCBL.list<-SplitObject(SCBL,split.by="treatment")

#Do Seurat's scTransform plug in in liu of the normal normalization and feature selection (~20 minutes on cluster)
for (i in 1:length(SCBL.list)) {
    SCBL.list[[i]] <- SCTransform(SCBL.list[[i]], verbose = FALSE)
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that 
#all necessary Pearson residuals have been calculated
SCBL.features <- SelectIntegrationFeatures(object.list = SCBL.list, nfeatures = 3000)
SCBL.list <- PrepSCTIntegration(object.list = SCBL.list, anchor.features = SCBL.features, verbose = TRUE)

#Now find integration "anchors" (uses CCA) based on 30 dimensions (~3 hours on cluster)
SCBL.anchors <- FindIntegrationAnchors(SCBL.list, anchor.features = SCBL.features, dims = 1:30, normalization.method = "SCT")

#Now integrate the datasets (~2 hours on cluster)
SCBL_int <- IntegrateData(anchorset = SCBL.anchors, dims = 1:30, normalization.method = "SCT")

#Switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(SCBL_int) <- "integrated"

#Run pca
SCBL_int<- RunPCA(SCBL_int, npcs = 30, verbose = TRUE)

#Elbow
ElbowPlot(SCBL_int,ndims=30) #How many PCs

#Look at loadings for each dim
VizDimLoadings(object = SCBL_int, dims = 1:2)

#Look at PCs
DimPlot(SCBL_int, reduction="pca")

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(SCBL_int, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(SCBL_int, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(SCBL_int, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(SCBL_int, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(SCBL_int, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()

#Don't forget to normalize the RNA data in case I need it
SCBL_int<-NormalizeData(SCBL_int, verbose = FALSE,assay="RNA")

#Clustering
#Input
nPCs<-20
colorVec<-c("blue","royalblue4","mediumpurple3","paleturquoise3","paleturquoise2","lightsteelblue1","lightskyblue","purple",
	"slategrey","darkseagreen","turquoise4","peachpuff4","brown","brown2","yellow","deeppink","deeppink4","lightcoral","pink","yellowgreen",
	"rosybrown","black","darkgoldenrod","chocolate","gold","greenyellow","yellow4","yellow3","wheat","azure3","sienna4",
	"darkolivegreen","darkolivegreen1","forestgreen","darkslategray","springgreen")
	
#Run umap
SCBL_int <- RunUMAP(SCBL_int, reduction.use = "pca", dims=1:nPCs, n.neighbors = 18, min.dist = 0.6)

#Now do clustering
SCBL_int <- FindNeighbors(object = SCBL_int, dims = 1:nPCs, k.param = 10)
SCBL_int <- FindClusters(SCBL_int, reduction.type="pca", resolution=0.15, algorithm=1) #0.26 for 11 clusters, 0.15 for 8 clusters

#pdf("SCBL_UMAP_clusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(SCBL_int, reduction='umap', label=T, pt.size = 0.3, cols=colorVec)
dev.off()

#Plot by donor, nGene, etc
p1<-DimPlot(SCBL_int, reduction = "umap", group.by = "treatment",pt.size=2, cols=c("midnightblue","red")) + NoAxes()
p2<-FeaturePlot(SCBL_int, "nGene",cols=c("gray", "blue"),pt.size=1,min.cutoff="q5",max.cutoff="q95") + NoAxes()
p3<-FeaturePlot(SCBL_int, "propMT",cols=c("gray", "blue"),pt.size=1,min.cutoff="q5",max.cutoff="q95") + NoAxes()
p4<-FeaturePlot(SCBL_int, "nUMI",cols=c("gray", "blue"),pt.size=1,min.cutoff="q5",max.cutoff="q95") + NoAxes()
pdf("SCBL_UMAP_treatment.pdf",height=10,width=14)
plot_grid(p1)
dev.off()
pdf("SCBL_UMAP_nGene.pdf",height=5,width=8)
plot_grid(p2)
dev.off()
pdf("SCBL_UMAP_propMT.pdf",height=5,width=8)
plot_grid(p3)
dev.off()
pdf("SCBL_UMAP_nUMI.pdf",height=5,width=8)
plot_grid(p4)
dev.off()

#Rename clusters to be previous names (8 clusters)
treat.col<-data.frame("clusters_8"=SCBL_int@meta.data$integrated_snn_res.0.15,row.names=rownames(SCBL_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "5")],]<-"c1"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "6")],]<-"c1c"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "1")],]<-"c2"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "2")],]<-"c4a"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "0")],]<-"c4b"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "3")],]<-"c5"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "7")],]<-"c6"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.15 == "4")],]<-"c8"
SCBL_int<-AddMetaData(SCBL_int,metadata=treat.col,col.name="clusters_8")

#Colors
colorVec<-c("#1619BA","mediumpurple3","#3281FF","tomato","pink","orange","tan","green3")

#pdf("T15_UMAP_clusters_8.pdf",height=5,width=7)
dev.new(height=5,width=6)
DimPlot(SCBL_int, reduction='umap', label=T, pt.size = 0.8, cols=colorVec,group.by="clusters_8") + NoLegend()

#Rename clusters to be previous names (11 clusters)
treat.col<-data.frame("clusters_11"=SCBL_int@meta.data$integrated_snn_res.0.23,row.names=rownames(SCBL_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "4")],]<-"c1"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "9")],]<-"c1c"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "1")],]<-"c2"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "2")],]<-"c4a"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "0")],]<-"c4b"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "7")],]<-"c4c"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "3")],]<-"c5a"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "8")],]<-"c5b"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "10")],]<-"c6"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "6")],]<-"c8a"
treat.col[rownames(SCBL_int@meta.data)[which(SCBL_int@meta.data$integrated_snn_res.0.23 == "5")],]<-"c8b"
SCBL_int<-AddMetaData(SCBL_int,metadata=treat.col,col.name="clusters_11")

#Colors
colorVec<-c("#1619BA","mediumpurple3","#3281FF","tomato","pink","violet","orange","saddlebrown","tan","greenyellow","green3")

#pdf("T15_UMAP_clusters_11.pdf",height=5,width=7)
dev.new(height=5,width=6)
DimPlot(SCBL_int, reduction='umap', label=T, pt.size = 0.8, cols=colorVec,group.by="clusters_11") + NoLegend()

################### DEGs and Enrichments to define major clusters

#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.6.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp_8.r"
#Logistic regression, with treatment as a latent variable
Idents(SCBL_int) <- SCBL_int@meta.data$clusters_8
clusterList<-sort(as.character(unique(Idents(SCBL_int))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(SCBL_int, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, 
    	only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
}

#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.6.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp_11.r"
Idents(SCBL_int) <- SCBL_int@meta.data$clusters_11
clusterList<-sort(as.character(unique(Idents(SCBL_int))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(SCBL_int, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, 
    	only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
}
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c8a", ident.2="c8b", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "wilcox")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c8b", ident.2="c8a", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "wilcox")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c5a", ident.2="c5b", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c5b", ident.2="c5a", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c4b", ident.2="c4c", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c4c", ident.2="c4b", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c4c", ident.2="c5a", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
markers.list[[length(markers.list) + 1]]<-FindMarkers(SCBL_int, ident.1="c5a", ident.2="c4c", min.pct=0.1, 
	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "treatment")
	
######### Consolidate and do Enrichr (8 clusters)
Idents(SCBL_int) <- SCBL_int@meta.data$clusters_8
comparisonVec<-sort(as.character(unique(Idents(SCBL_int))))
#Run OrganizeDF_list
SCBL_1againstAll_8_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec)
#Write all to single txt file
write.table(SCBL_1againstAll_8_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/SCBL_1againstAll_8_DEGs.txt")
#Write to excel
writeDEGsToExcel(DEG_table = SCBL_1againstAll_8_DEGs, savedFile = "DEG_Tables/SCBL_1againstAll_8_DEGs.xlsx")
#Export table for Enrichr
SCBL_1againstAll_8_DEGs_forEnrichr<-SCBL_1againstAll_8_DEGs[which(SCBL_1againstAll_8_DEGs$p_adj_FDR < 0.05),] #Only good DEGs
#Write all to single txt file
write.table(SCBL_1againstAll_8_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/SCBL_1againstAll_8_DEGs_forEnrichr.txt")
#Do enrichments
DEG_table<-read.table("Enrichr/SCBL_1againstAll_8_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"SCBL_1againstAll_8_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

######### Consolidate and do Enrichr (11 clusters)
Idents(SCBL_int) <- SCBL_int@meta.data$clusters_11
comparisonVec<-sort(as.character(unique(Idents(SCBL_int))))
comparisonVec<-append(comparisonVec,c("c8a_v_c8b","c8b_v_c8a","c5a_v_c5b","c5b_v_c5a","c4b_v_c4c","c4c_v_c4b",
	"c4c_v_c5a","c5a_v_c4c"))
#Run OrganizeDF_list
SCBL_1againstAll_11_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec)
#Write all to single txt file
write.table(SCBL_1againstAll_11_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/SCBL_1againstAll_11_DEGs.txt")
#Write to excel
writeDEGsToExcel(DEG_table = SCBL_1againstAll_11_DEGs, savedFile = "DEG_Tables/SCBL_1againstAll_11_DEGs.xlsx")
#Export table for Enrichr
SCBL_1againstAll_11_DEGs_forEnrichr<-SCBL_1againstAll_11_DEGs[which(SCBL_1againstAll_11_DEGs$p_adj_FDR < 0.05),] #Only good DEGs
#Write all to single txt file
write.table(SCBL_1againstAll_11_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/SCBL_1againstAll_11_DEGs_forEnrichr.txt")
#Do enrichments
DEG_table<-read.table("Enrichr/SCBL_1againstAll_11_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"SCBL_1againstAll_11_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

######################### Visualize covid gene expression across cell types in response to IL-13
	
#Get vector to give groups of cells
clust.col<-data.frame("clusters_11_treat"=paste(SCBL_int@meta.data$clusters_11,SCBL_int@meta.data$treatment,sep="_"),
	row.names=rownames(SCBL_int@meta.data))
clust.col[,1]<-factor(clust.col[,1],levels=c(levels(clust.col[,1])[1:2],levels(clust.col[,1])[5:6],
	levels(clust.col[,1])[3:4],levels(clust.col[,1])[7:21]))
SCBL_int<-AddMetaData(SCBL_int,metadata=clust.col,col.name="clusters_11_treat")

#Get dataset
sct<-as.matrix(GetAssayData(SCBL_int,assay="RNA"))

#Now plot
#pdf("Coronovirus_genes_SPDEFKO_scramble_singleCell.pdf",width=4,height=5)
#BSA plots
dev.new(width=5.3,height=6.3)
par(mfrow=c(2,1),bty="l")
beeswarm(sct["ACE2",]~SCBL_int@meta.data$clusters_11_treat, 
	col= c("black","black"), pch = 16, cex = 0.5, method="swarm", corral="random",las=1,
	ylab="Normalized expression",xlab="",labels=c("",""),add=F)
boxplot(sct["ACE2",]~SCBL_int@meta.data$clusters_11_treat,
	las=1,main="ACE2",col=c("grey","tomato"),names=F,outline=F,
	ylim=c(0,max(sct["ACE2",])),add=T,boxwex=0.35,medcol="white")	
vioplot(sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[1])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[2])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[3])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[4])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[5])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[7])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[8])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[9])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[10])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[13])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[14])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[15])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[16])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[20])],
	sct["ACE2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[21])],
	col="grey",drawRect=F,add=T)
text(x=c(1:21),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c(rep(c("CTRL","IL-13"),9),"CTRL","CTRL","IL-13"),
	srt=45,adj=0.8,xpd=T,,cex=0.7)
#0.0048 - down

beeswarm(sct["TMPRSS2",]~SCBL_int@meta.data$clusters_11_treat, 
	col= c("black","black"), pch = 16, cex = 0.5, method="swarm", corral="random",las=1,
	ylab="Normalized expression",xlab="",labels=c("",""),add=F)
boxplot(sct["TMPRSS2",]~SCBL_int@meta.data$clusters_11_treat,
	las=1,main="TMPRSS2",col=c("grey","tomato"),names=F,outline=F,
	ylim=c(0,max(sct["TMPRSS2",])),add=T,boxwex=0.35,medcol="white")
vioplot(sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[1])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[2])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[3])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[4])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[5])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[6])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[7])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[8])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[9])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[10])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[11])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[12])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[13])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[14])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[15])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[16])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[17])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[18])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[19])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[20])],
	sct["TMPRSS2",which(SCBL_int@meta.data$clusters_11_treat==levels(SCBL_int@meta.data$clusters_11_treat)[21])],
	col="grey",drawRect=F,add=T)	
text(x=c(1:21),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c(rep(c("CTRL","IL-13"),9),"CTRL","CTRL","IL-13"),
	srt=45,adj=0.8,xpd=T,cex=0.7)
#0.0048 - down

#Get p values
##### ACE2
pvals_ACE2<-list()
clusterList<-levels(SCBL_int@meta.data$clusters_11_treat)
for(i in c(1,3,5,7,9,13,15,20)){
	pvals_ACE2[[length(pvals_ACE2) + 1]]<-FindMarkers(SCBL_int, features = "ACE2", ident.1=clusterList[i+1], ident.2=clusterList[i],min.pct=-1, 
	logfc.threshold=0, only.pos=F, assay="RNA", test.use = "wilcox")
}

##### TMPRSS2
pvals_TMP<-list()
for(i in c(1,3,5,7,9,11,13,15,17,20)){
	pvals_TMP[[length(pvals_TMP) + 1]]<-FindMarkers(SCBL_int, features = "TMPRSS2", ident.1=clusterList[i+1], ident.2=clusterList[i],min.pct=-1, 
	logfc.threshold=0, only.pos=F, assay="RNA", test.use = "wilcox")
}