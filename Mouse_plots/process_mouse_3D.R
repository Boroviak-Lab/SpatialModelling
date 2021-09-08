library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")

set.seed(1)

#Script for loading and normalisnig data for mouse analyses
saveext = "./Mouse3D/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

BS<-read.table("../Data/Mouse3D_ALLKEY.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("../Data/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)

mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
isinvit <- BS$NotRef
labs <- BS$ID2
labs <- labs[which(isinvit>-1)]

raw_counts3<- raw_counts3[,which(isinvit>-1)]
mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data$Emb <- BS$Embryo
mouse_data$species <- "mouse"
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)
Idents(mouse_data) <- BS$Embryo #labs
mouse_data$ID <- labs
mouse_data$Emb <- Idents(mouse_data)
#mouse_data <- subset(mouse_data, idents = c("Myo_CS7","ReGland_CS5","ReGland_CS7","Gland_CS5","Gland_CS6","Gland_CS7","4-cell_CS2","8-cell_CS2","Zy_CS1","cMor_CS3"), invert = TRUE)
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 20000)
mouse_data <- ScaleData(mouse_data, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mouse_data<- RunPCA(mouse_data, npcs = 30, verbose = FALSE)
mouse_data <- RunUMAP(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- RunTSNE(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- FindNeighbors(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- FindClusters(mouse_data, resolution = 0.5)
mouse_data$Cl <- Idents(mouse_data)


p <- DimPlot(mouse_data, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_ClLab","_4.pdf",sep=""),width = 10, height = 8,p)

Idents(mouse_data) <- mouse_data$Emb

p <- DimPlot(mouse_data, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) #+NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs","_4.pdf",sep=""),width = 10, height = 8,p)

FeaturePlot(mouse_data,  reduction = "umap", features = "Sox2", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/DimRed/SOX2","_4.pdf",sep=""),width = 10, height = 8)

FeaturePlot(mouse_data,  reduction = "umap", features = "Pou5f1", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/DimRed/OCT4","_4.pdf",sep=""),width = 10, height = 8)

FeaturePlot(mouse_data,  reduction = "umap", features = "Sox17", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/DimRed/SOX17","_4.pdf",sep=""),width = 10, height = 8)


raw_counts3<-read.table("../Data/featurecountsAllMouse3D_all.csv",sep=",",header = T, row.names=1)
mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
isinvit <- BS$IsEmb
labs <- BS$ID2
labs2 <- BS$PosID
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
raw_counts3<- raw_counts3[,which(isinvit>0)]
mouse_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data$species <- "mouse"
mouse_data <- subset(mouse_data, subset = nFeature_RNA > 0)
Idents(mouse_data) <- labs
mouse_data$ID <- labs
mouse_data$ID2 <- labs2
#mouse_data <- subset(mouse_data, idents = c("Myo_CS7","ReGland_CS5","ReGland_CS7","Gland_CS5","Gland_CS6","Gland_CS7","4-cell_CS2","8-cell_CS2","Zy_CS1","cMor_CS3"), invert = TRUE)
mouse_data <- NormalizeData(mouse_data, verbose = FALSE)
mouse_data <- FindVariableFeatures(mouse_data, selection.method = "vst", nfeatures = 20000)
mouse_data <- ScaleData(mouse_data, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
mouse_data<- RunPCA(mouse_data, npcs = 30, verbose = FALSE)
mouse_data <- RunUMAP(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- RunTSNE(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- FindNeighbors(mouse_data, reduction = "pca", dims = 1:20)
mouse_data <- FindClusters(mouse_data, resolution = 0.5)
mouse_data$Cl <- Idents(mouse_data)


p <- DimPlot(mouse_data, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_ClLab2","_4.pdf",sep=""),width = 10, height = 8,p)

Idents(mouse_data) <- mouse_data$ID

p <- DimPlot(mouse_data, pt.size = 4, reduction = "umap", label = TRUE, repel = TRUE) +NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs2","_4.pdf",sep=""),width = 45, height = 40,p)

FeaturePlot(mouse_data,  reduction = "umap", features = "Sox2", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/Markers/SOX2","2_4.pdf",sep=""),width = 10, height = 8)

FeaturePlot(mouse_data,  reduction = "umap", features = "Pou5f1", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/Markers/OCT4","2_4.pdf",sep=""),width = 10, height = 8)

FeaturePlot(mouse_data,  reduction = "umap", features = "Sox17", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 8)
ggsave(filename=paste(saveext,"/Markers/SOX17","2_4.pdf",sep=""),width = 10, height = 8)

saveRDS(mouse_data, file = paste(saveext,"mouse3D.rds",sep=""))

write.csv(as.data.frame(GetAssayData(object = mouse_data, slot = 'data')), file=paste(saveext,"/RNA_assay.csv",sep=""))
write.csv(as.data.frame(Idents(mouse_data)), file=paste(saveext,"/Header.csv",sep=""))

Idents(mouse_data) <- mouse_data$ID2
AvExp <- AverageExpression(object = mouse_data)
write.csv(as.data.frame(AvExp$RNA*100), file=paste(saveext,"/AverageExp.csv",sep=""))

#Write data out, this will be used for 3D models
Dsubs <- AvExp$RNA[c("Sox2","Pou5f1","Nanog","Prdm14","Fbxo2","Tdgf1","Sox17","Mixl1","Hnf4a","T","Eomes","Mixl1","Foxa2","Otx2"),] #GetAssayData(object = AvExp , slot = 'data')[c("Sox2","Pou5f1","Nanog","Prdm14","Fbxo2","Tdgf1","Sox17","Mixl1","Hnf4a","T"),]
write.csv(as.data.frame(Dsubs), file=paste(saveext,"/AverageExpSubs.csv",sep=""))

