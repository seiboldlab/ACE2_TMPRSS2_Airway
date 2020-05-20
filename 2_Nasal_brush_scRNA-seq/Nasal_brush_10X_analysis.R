## R/3.5.1

library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
library(future)
library(cowplot)

### set parallelization settings
plan("multiprocess", workers = 8)

#### sample name
sample_name <- "NasalBrush"

dat <- Read10X(data.dir = "/Seibold/proj/10X_Mapping/200211_A00405_0207_AH2FMFDSXY/NasalBrush/outs/filtered_feature_bc_matrix")

sDat <- CreateSeuratObject(counts = dat, min.cells = 3, min.features = 100, project = sample_name)

mito.genes <- grep(pattern = "^MT-|^MRPS|^MRPL", x = rownames(x = sDat), value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(sDat, slot='counts')[mito.genes,])/Matrix::colSums(x=GetAssayData(sDat, slot='counts'))

ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(x = sDat), value = TRUE)
percent.ribo <- Matrix::colSums(x=GetAssayData(sDat, slot='counts')[ribo.genes,])/Matrix::colSums(x=GetAssayData(sDat, slot='counts'))

sDat[['percent.mito']] <- percent.mito
sDat[['percent.ribo']] <- percent.ribo

pdf(paste0("./IMAGES/", sample_name, "_init.QC.plots.pdf"))
VlnPlot(sDat, c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol=2)
hist(sDat@meta.data$nFeature_RNA, breaks=50)
hist(sDat@meta.data$nCount_RNA, breaks=50)
hist(sDat@meta.data$percent.mito, breaks=50)
hist(sDat@meta.data$percent.ribo, breaks=50)
g <- FeatureScatter(sDat, "nCount_RNA", "percent.mito")
g <- g + geom_hline(yintercept=0.25, color="blue")
print(g)
FeatureScatter(sDat, "nCount_RNA", "percent.ribo")
FeatureScatter(sDat, "nCount_RNA", "nFeature_RNA")
dev.off()

sDat <- subset(x = sDat, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mito < 0.25)
sDat <- NormalizeData(sDat)
sDat <- SCTransform(sDat, vars.to.regress = "percent.mito", variable.features.n=5000, verbose = FALSE)

#====================#
# Initial clustering #
#====================#

sDat <- RunPCA(object = sDat, npcs = 30, verbose = FALSE)

pdf(paste0("./IMAGES/", sample_name, "_seurat.PCA.pdf"))
VizDimLoadings(object = sDat, dims = 1:2)
DimPlot(sDat, reduction="pca")
DimHeatmap(sDat, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(sDat, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(sDat, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(sDat, dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(sDat, dims = 25:30, cells = 500, balanced = TRUE)
ElbowPlot(sDat, ndims=30)
dev.off()

nPCAs <- 18

sDat <- FindNeighbors(object = sDat, dims = 1:nPCAs)

sDat <- RunUMAP(object = sDat, reduction = "pca", dims = 1:nPCAs)
#sDat <- RunTSNE(object = sDat, reduction = "pca", dims = 1:nPCAs)

sDat0 <- FindClusters(sDat, reduction.type="pca", resolution=0.2, algorithm=3)
sDat1 <- FindClusters(sDat, reduction.type="pca", resolution=0.4, algorithm=3)
sDat2 <- FindClusters(sDat, reduction.type="pca", resolution=0.6, algorithm=3)
sDat3 <- FindClusters(sDat, reduction.type="pca", resolution=0.8, algorithm=3)
sDat4 <- FindClusters(sDat, reduction.type="pca", resolution=1, algorithm=3)

pdf(paste0("./IMAGES/", sample_name, "_sDat.UMAP.pdf"))
DimPlot(sDat0, reduction='umap')
DimPlot(sDat1, reduction='umap')
DimPlot(sDat2, reduction='umap')
DimPlot(sDat3, reduction='umap')
DimPlot(sDat4, reduction='umap')
FeaturePlot(object=sDat, features=c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"))
dev.off()

pdf(paste0("./IMAGES/", sample_name, "_QC_VlnPlot.pdf"), width=18)
VlnPlot(object=sDat0, features="nCount_RNA")
VlnPlot(object=sDat0, features="nFeature_RNA")
VlnPlot(object=sDat0, features="percent.mito")
VlnPlot(object=sDat0, features="percent.ribo")
dev.off()

#==============#
# Find markers #
#==============#

# Choose resolution=0.8 

markers.list <- list()
nClusters <- 18

for(i in 1:nClusters) {
    markers.list[[i]] <- FindMarkers(sDat3, ident.1=(i-1), min.pct=0.1, logfc.threshold=0.25,  only.pos=T, assay="RNA")
    markers.list[[i]] <- markers.list[[i]][markers.list[[i]]$p_val_adj < 0.05,]
}

wb <- createWorkbook()
for(i in 1:nClusters) {
    if(dim(markers.list[[i]])[1]!=0) {
        addWorksheet(wb, sheetName=paste0("cl_",i-1))
        writeData(wb, x=markers.list[[i]], sheet=paste0("cl_",i-1), rowNames=T)
    }
}
saveWorkbook(wb, file=paste0("./TABLES/", sample_name, "_markers.xlsx"), overwrite=T)

#====================================#
# define epithelial and immune cells #
#====================================#
immune_cells <- rownames(sDat3@meta.data[sDat3@meta.data$seurat_clusters %in% c(14,17),])
ion_cells <- rownames(sDat3@meta.data[sDat3@meta.data$seurat_clusters==16,])
epi_cells <- setdiff(rownames(sDat3@meta.data), c(ion_cells, immune_cells))

# second filtering for epithelial (without ionocytes) cells 
sDat_epi <- subset(x = sDat3, cells=epi_cells)
sDat_epi <- subset(x = sDat_epi, subset = nFeature_RNA > 1000)

#===========================#
# Second step of clustering #
#===========================#

final_cell_names <- c(rownames(sDat_epi@meta.data), ion_cells, immune_cells)
sDat_final <- subset(sDat3, cells=final_cell_names)
sDat_final <- NormalizeData(sDat_final)
sDat_final <- SCTransform(sDat_final, vars.to.regress = "percent.mito", variable.features.n=5000, verbose = FALSE)

sDat_final <- RunPCA(object = sDat_final, npcs = 30, verbose = FALSE)

pdf(paste0("./IMAGES/", sample_name, "_seurat_sDat_final.PCA.pdf"))
VizDimLoadings(object = sDat_final2, dims = 1:2)
DimPlot(sDat_final, reduction="pca")
DimHeatmap(sDat_final, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(sDat_final, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(sDat_final, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(sDat_final, dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(sDat_final, dims = 25:30, cells = 500, balanced = TRUE)
ElbowPlot(sDat_final, ndims=30)
dev.off()

nPCAs <- 30

sDat_final <- FindNeighbors(object = sDat_final, dims = 1:nPCAs, k.param=10)
sDat_final <- RunUMAP(object = sDat_final, reduction = "pca", dims = 1:nPCAs)

sDat_final_0 <- FindClusters(sDat_final, reduction.type="pca", resolution=0.2, algorithm=1)
sDat_final_1 <- FindClusters(sDat_final, reduction.type="pca", resolution=0.3, algorithm=1)
sDat_final_2 <- FindClusters(sDat_final, reduction.type="pca", resolution=0.4, algorithm=1)
sDat_final_3 <- FindClusters(sDat_final, reduction.type="pca", resolution=0.5, algorithm=1)
sDat_final_4 <- FindClusters(sDat_final, reduction.type="pca", resolution=0.6, algorithm=1)

cell_type_cols=RColorBrewer::brewer.pal(12, "Paired")[c(1:10,12)]

# Find markers 

markers.list <- list()
nClusters <- 15

for(i in 1:nClusters) {
    markers.list[[i]] <- FindMarkers(sDat_final_2, ident.1=(i-1), min.pct=0.1, logfc.threshold=0.25,  only.pos=T, assay="RNA")
    markers.list[[i]] <- markers.list[[i]][markers.list[[i]]$p_val_adj < 0.05,]
}

wb <- createWorkbook()
for(i in 1:nClusters) {
    if(dim(markers.list[[i]])[1]!=0) {
        addWorksheet(wb, sheetName=paste0("cl_",i-1))
        writeData(wb, x=markers.list[[i]], sheet=paste0("cl_",i-1), rowNames=T)
    }
}
saveWorkbook(wb, file=paste0("./TABLES/", sample_name, "_sDat_final_markers.xlsx"), overwrite=T)

### immune
T_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==12,])
mast_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==10,])
DCs <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==14,])

### epi
Ionocytes_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==11,])
Basal_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==0,])
Secretory_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==2 | sDat_final_2@meta.data$seurat_clusters==3,])
Ciliated_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==1 | sDat_final_2@meta.data$seurat_clusters==9,])
LQ_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==5 | sDat_final_2@meta.data$seurat_clusters==7,])
MUC5B_pos_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==4,])
SMG_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==6,])
MUC5AC_pos_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==8,])
Prolif_Basal_cells <- rownames(sDat_final_2@meta.data[sDat_final_2@meta.data$seurat_clusters==13,])

sDat_final_2@meta.data$cell_type <- "Basal/early sec."
sDat_final_2@meta.data[T_cells, "cell_type"] <- "T cells"
sDat_final_2@meta.data[mast_cells, "cell_type"] <- "Mast cells"
sDat_final_2@meta.data[DCs, "cell_type"] <- "DCs"
sDat_final_2@meta.data[ionocytes_cells, "cell_type"] <- "Ionocytes"
sDat_final_2@meta.data[Secretory_cells, "cell_type"] <- "Secretory"
sDat_final_2@meta.data[Ciliated_cells, "cell_type"] <- "Ciliated"
sDat_final_2@meta.data[MUC5B_pos_cells, "cell_type"] <- "Mucus sec."
sDat_final_2@meta.data[SMG_cells, "cell_type"] <- "Glandular mucus sec."
sDat_final_2@meta.data[MUC5AC_pos_cells, "cell_type"] <- "MUC5AC mucus sec."
sDat_final_2@meta.data[Prolif_Basal_cells, "cell_type"] <- "Prolif. Basal"
sDat_final_2@meta.data[LQ_cells, "cell_type"] <- "Indeterminate"

Idents(sDat_final_2) <- sDat_final_2@meta.data$cell_type

cell_type_cols <- c("#bebebe", RColorBrewer::brewer.pal(12, "Paired")[c(1:10,12)])

# Fig 1a
pdf(paste0("./IMAGES/", sample_name, "_sDat_final.UMAP.pdf"))
DimPlot(sDat_final_2, reduction='umap', cols=cell_type_cols)
dev.off()

# Fig 1b
tiff(paste0("./IMAGES/", sample_name, "_sDat_final.ACE2_VlnPlot.tiff"), res=300, width=6*300, height=4*300, compression="lzw")
VlnPlot(sDat_final_2, c("ACE2"), assay="RNA")
dev.off()
# Fig 1c
tiff(paste0("./IMAGES/", sample_name, "_sDat_final.TMPRSS2_VlnPlot.tiff"), res=300, width=6*300, height=4*300, compression="lzw")
VlnPlot(sDat_final_2, c("TMPRSS2"), assay="RNA")
dev.off()

