library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")

library(data.table)

#Script for analysing all good shots from the marmoset data

set.seed(1)

#              _______ Human
#     ________|
#     |       |_______ Macaque
#_____|
#|    |
#|    |_______________ Marmoset
#|
#|____________________ Mouse


#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "DimRed/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))


BS<-read.table("../Data/marmKey.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("../Data/marmData.csv",sep=",",header = T, row.names=1)
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
isinvit <- BS$All * BS$QC
labs1 <- BS$Carnegie.Stage
labs2 <- BS$Annotation2
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]

#Load raw counts
marmoset_data <- CreateSeuratObject(counts = raw_counts3[,which(isinvit>0)], assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "marmoset"
marmoset_data$divergence1 <- "newworld"
marmoset_data$divergence2 <- "primate"
marmoset_data$CS <- labs
marmoset_data$CT <- labs2
#Don't really neet to subset already QCd
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
#Normalise the data
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
#Get variable genes (reduce number)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)

#Rename to keep variable convention for the join species modelling
mammal.combined <- marmoset_data
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

#Dimesionality reduction and clustering
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20, perplexity = 50)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 10) #k.param reduce for later applications? 
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)

write.csv(as.data.frame( GetAssayData(object = mammal.combined, slot = 'scale.data') ), file=paste(saveext,"/NormData.csv",sep=""))
write.csv(as.data.frame( GetAssayData(object = mammal.combined, slot = 'data') ), file=paste(saveext,"/RawData.csv",sep=""))

#write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsCl.csv",sep=""))
#write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNE.csv",sep=""))
#write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddings.csv",sep=""))
#write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey.csv",sep=""))
#write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCA.csv",sep=""))

ElbowPlot(mammal.combined,ndims = 30)
ggsave(filename=paste(saveext,"/DimRed/Components",".pdf",sep=""),width = 4, height = 4)

VizDimLoadings(mammal.combined, dims = 1:6, nfeatures = 30, col = "blue",reduction = "pca", projected = FALSE, balanced = FALSE,ncol = NULL, combine = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_loadings",".pdf",sep=""),width = 8, height = 20)

#Identify markers based on cluster
#mammal.combined.markerscl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 100,limitsize = FALSE)
#write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))

Cls <- Idents(mammal.combined)
Idents(mammal.combined) <- labs2

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2")
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

stagelab <- LETTERS[mammal.combined$orig.ident]
stagelab[1:length(stagelab)] <- 0
stagelab[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
stagelab[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
stagelab[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
stagelab[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
stagelab[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
stagelab[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
stagelab[which(Idents(mammal.combined)=="VE_CS4")] <- 1
stagelab <- as.factor(stagelab)

names(stagelab)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = stagelab,col.name = 'cell.orig')

DimPlot(mammal.combined, reduction = "pca",  cols = coluse, label = TRUE,repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_base_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)

DimPlot(mammal.combined, reduction = "tsne", cols = coluse, label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)
