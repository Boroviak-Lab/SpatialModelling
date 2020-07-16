library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "CrossSpeciesAlignment/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#List of cell cycle markers loaded in with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Load marmoset data key
cyBS<-read.table("../Data/cyKey.csv",sep=",", header = T, row.names=1)
BS<-read.table("../Data/marmKey.csv",sep=",",header = T, row.names=1)

#Load monkey
raw_counts1<-read.table("../Data/cyData.csv",sep=",",header = T, row.names=1)
cynomolgous_data <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
cyisinvit <- cyBS$AllNoES
cylabs <- cyBS$LABEL_4 
cylabs <- cylabs[which(cyisinvit>0)]
cylabs2 <- cyBS$LABEL_2
cylabs2 <- cylabs2[which(cyisinvit>0)]
raw_counts1<- raw_counts1[,which(cyisinvit>0)]
cynomolgous_data <- CreateSeuratObject(counts = raw_counts1, assay = "RNA",min.cells = 0, min.features = 0)
cynomolgous_data$species <- "2) Cynomolgous"
cynomolgous_data$divergence1 <- "Cyno"
cynomolgous_data <- subset(cynomolgous_data, subset = nFeature_RNA > 0)
cynomolgous_data$divergence2 <- "primate"
cynomolgous_data <- NormalizeData(cynomolgous_data, verbose = FALSE)
cynomolgous_data <- FindVariableFeatures(cynomolgous_data, selection.method = "vst", nfeatures = 20000)
#cynomolgous_data <- ScaleData(cynomolgous_data, verbose = FALSE)
Idents(cynomolgous_data) <- cylabs
cynomolgous_data$ID <- cylabs
cynomolgous_data$ID2 <- cylabs2

#Load the marmoset data ... to get ordering
raw_counts3<-read.table("../Data/marmData.csv",sep=",",header = T, row.names=1)
isinvit <- BS$vsCynoMin * BS$QC
labs <- BS$Annotation2
labs <- labs[which(isinvit>0)]
raw_counts3<- raw_counts3[,which(isinvit>0)]
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "1) Marmoset"
marmoset_data$divergence1 <- "Marm"
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)
marmoset_data <- ScaleData(marmoset_data, verbose = FALSE)
Idents(marmoset_data) <- labs
marmoset_data$ID <- labs
marmoset_data$ID2 <- labs

#Load human dataset
hBS<-read.table("../Data/humanKey.csv",sep=",",header = T, row.names=1)
#Cell labels
hisinvit <- which(hBS$Type5 != "MISSING" )
hlabs <- hBS$Type5[hisinvit]
hlabs2 <- hBS$Type2[hisinvit]
raw_countsA1<-read.table("../Data/humanData.csv",sep=",",header = T, row.names=1)
raw_countsA1 <- raw_countsA1[,2:dim(raw_countsA1)[2]]
human_dataA1 <- CreateSeuratObject(counts = na.omit(raw_countsA1[,hisinvit]), assay = "RNA",min.cells = 0, min.features = 0)
human_dataA1$species <- "3) Human (in vitro)"
human_dataA1$divergence1 <- "Human_invit"
human_dataA1 <- subset(human_dataA1, subset = nFeature_RNA > 0)
human_dataA1 <- NormalizeData(human_dataA1, verbose = FALSE)
human_dataA1 <- FindVariableFeatures(human_dataA1, selection.method = "vst", nfeatures = 20000)
human_dataA1$ID <- hlabs
human_dataA1$ID2 <- hlabs2
Idents(human_dataA1) <- hlabs

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(marmoset_data, cynomolgous_data, human_dataA1), dims = 1:20, anchor.features = 4000)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#mammal.combined <- CellCycleScoring(mammal.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

types <- Idents(mammal.combined)

#Print our loadings on PCA
Loadis <- as.data.frame(Loadings(mammal.combined, reduction = "pca")[, 1:2])
p <- ggplot(Loadis, aes(x=PC_1, y=PC_2)) + geom_point() + geom_text(label=rownames(Loadis))+ theme_classic()
ggsave(filename=paste(saveext,"/DimRed/PCA_Loadings_2D",".pdf",sep=""),width = 40, height = 40,p)

cType <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS3","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo_CS7","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7")
BaseCol <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#967700","#E6E600","#E6E600","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2")

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(Idents(mammal.combined)=="ICM_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Hyp_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Epi_CS3")] <- 1
cluster_letters[which(Idents(mammal.combined)=="Tb_CS4")] <- 1
cluster_letters[which(Idents(mammal.combined)=="VE_CS4")] <- 1
cluster_letters <- as.factor(cluster_letters)

names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')

p <- ElbowPlot(mammal.combined,ndims = 20, reduction = "pca")
ggsave(filename=paste(saveext,"/DimRed/Components",".pdf",sep=""),width = 4, height = 8, p)

p <- DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_4.pdf",sep=""),width = 42, height = 8,p)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_12.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 10, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = coluse, pt.size = 12, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_12.pdf",sep=""),width = 42, height = 8)


#Align to other annotations
Idents(mammal.combined) <- mammal.combined$ID2

colKey1 <- levels(factor(mammal.combined$ID))
ID1 <- factor(mammal.combined$ID) #factor(c(as.character(labs), as.character(cylabs),  as.character(hlabs)   ) )
Idents(mammal.combined) <- factor(as.character(mammal.combined$ID2)) #factor(c(as.character(labs), as.character(cylabs2),  as.character(hlabs2)   ) )
colKey2 <- levels(Idents(mammal.combined))
ID2 <- mammal.combined$ID2 #factor(c(as.character(labs), as.character(cylabs2),  as.character(hlabs2)   ) )

ind1 <- integer(length(colKey2))
for (i in 1:length(colKey2)) {
  ind1[i] <- which(colKey1==ID1[which(ID2==colKey2[i])[1]])
}

necols <- coluse[ind1]



DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 4, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 6, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 8, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 10, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 12, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab","_supp_12.pdf",sep=""),width = 42, height = 8)

DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 4, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_4.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 6, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_6.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 8, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_8.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 10, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_10.pdf",sep=""),width = 42, height = 8)
DimPlot(mammal.combined,  shape.by = 'cell.orig', cols = necols, pt.size = 12, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab","_supp_12.pdf",sep=""),width = 42, height = 8)




#Pull out marm dataset
split <- SplitObject(mammal.combined, split.by = "species")
AllMarkers <- FindAllMarkers(split$`1) Marmoset`, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
IntG <- intersect( intersect(rownames(marmoset_data), rownames(cynomolgous_data)), rownames(human_dataA) )

dir.create(paste(saveext,"/EpiCS3/",sep=""))
uGenes1 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Epi_CS3")] )
for (i in 1:length(uGenes1)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes1[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/EpiCS3/Markerscatter_",as.character(uGenes1[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/HypCS3/",sep=""))
uGenes2 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Hyp_CS3")] )
for (i in 1:length(uGenes2)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes2[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/HypCS3/Markerscatter_",as.character(uGenes2[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/TbCS3/",sep=""))
uGenes3 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS3")] )
for (i in 1:length(uGenes3)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes3[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/TbCS3/Markerscatter_",as.character(uGenes3[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/EmDisc_CS5/",sep=""))
uGenes4 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="EmDisc_CS5")] )
for (i in 1:length(uGenes4)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes4[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/EmDisc_CS5/Markerscatter_",as.character(uGenes4[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/SYS_CS5/",sep=""))
uGenes5 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="SYS_CS5")] )
for (i in 1:length(uGenes5)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes5[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/SYS_CS5/Markerscatter_",as.character(uGenes5[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/Tb_CS5/",sep=""))
uGenes6 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS5")] )
for (i in 1:length(uGenes6)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes6[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/Tb_CS5/Markerscatter_",as.character(uGenes6[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/ExMes_CS5/",sep=""))
uGenes7 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="ExMes_CS5")] )
for (i in 1:length(uGenes7)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes7[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
  ggsave(filename=paste(saveext,"/ExMes_CS5/Markerscatter_",as.character(uGenes7[i]),".pdf",sep=""),width = 35, height = 8)
}


#Example plot
FeaturePlot(mammal.combined, reduction = "pca", features = "POU5F1", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste(saveext,"/Markerscatter_POU5F1.pdf",sep=""),width = 35, height = 8)

FeaturePlot(mammal.combined, reduction = "pca", features = "T", combine=TRUE, cols =  c("lightgrey", "black"), split.by = "species", pt.size = 8)
ggsave(filename=paste(saveext,"/Markerscatter_TBXT.pdf",sep=""),width = 35, height = 8)

