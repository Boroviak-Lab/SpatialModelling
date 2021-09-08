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
library(pheatmap)

set.seed(1)

saveext = "./Maternal/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Load in the data
marmoset_dataInVivo2 <- readRDS('../Data/InVivo.rds')

annotationL <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3",
                 "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7",
                 "PGC_CS5","PGC_CS6","PGC_CS7",
                 "Am_CS5","Am_CS6","Am_CS7",
                 "Hyp_CS3",
                 "VE_CS5","VE_CS6","VE_CS7",
                 "SYS_CS5","SYS_CS6","SYS_CS7",
                 "ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7",
                 "Tb_CS3",
                 "Tb_CS5","Tb_CS6","Tb_CS7","Stroma","Gland","Remodelled","Myo")

#Sort the levels
Idents(marmoset_dataInVivo2) <- factor(Idents(marmoset_dataInVivo2), levels = annotationL)
#Look at maternal
marmoset_dataInVivo3 <- subset(marmoset_dataInVivo2,idents=c("Stroma","Remodelled","Gland","Myo"),invert=TRUE)

Maternal <- subset(marmoset_dataInVivo2,idents=c("Stroma","Remodelled","Gland","Myo"),invert=FALSE)
AllMarkerGenes <- FindAllMarkers(Maternal,test.use = "MAST", only.pos = TRUE)
write.csv(as.data.frame(AllMarkerGenes), file=paste(saveext,"/MaternalDE.csv",sep=""))

Maternal <- FindVariableFeatures(Maternal, selection.method = "vst", nfeatures = 20000)
Maternal <- ScaleData(Maternal, verbose = FALSE)
Maternal <- RunPCA(Maternal, npcs = 20, verbose = FALSE)
Maternal <- RunUMAP(Maternal, reduction = "pca", dims = 1:20, n.neighbors = 20)
Maternal <- RunTSNE(Maternal, reduction = "pca", dims = 1:20, perplexity = 30)
Maternal <- FindNeighbors(Maternal, reduction = "pca", dims = 1:20)

cType <- c("Gland_CS5","Gland_CS6","Gland","Stroma_CS5","Stroma","Remodelled_CS5","Remodelled_CS6","Remodelled","Stalk_CS6","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7")
BaseCol <- c("#B3B2B2","#A5A4A3","#969593","#DFDFDF","#D0D0D0","#878684","#797775","#6A6866","#754C24","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#603813","#E6E600","#BFBF04","#999903","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#754C24")

colind <- integer( length( levels(Idents(Maternal)) )  )
for (i in 1:length( levels(Idents(Maternal)) ) ) {
  colind[i] <- which(cType==levels(Idents(Maternal))[i])
}
coluse <- BaseCol[colind]

DimPlot(Maternal, cols = coluse, pt.size = 4, reduction = "umap",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Mat.pdf",sep=""),width = 10, height = 8)
DimPlot(Maternal, cols = coluse, pt.size = 4, reduction = "pca",label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Mat.pdf",sep=""),width = 10, height = 8)
DimPlot(Maternal, cols = coluse, pt.size = 4, reduction = "tsne", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_Mat.pdf",sep=""),width = 10, height = 8)


FeaturePlot(Maternal,  reduction = "umap", features = "SOX17", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_SOX17.pdf",sep=""),width = 16, height = 16)

