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

saveext = "./LineageInference/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

marmoset_dataInVivo2 <- readRDS('../Data/InVivo.rds')

#Colorcode for plotting
cType <- c("Other","Gland_CS5","Gland_CS6","Gland","Stroma_CS5","Stroma","Remodelled_CS5","Remodelled_CS6","Remodelled","Stalk_CS6","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Other")
BaseCol <- c("#f5f2d0","#B3B2B2","#A5A4A3","#969593","#DFDFDF","#D0D0D0","#878684","#797775","#6A6866","#754C24","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#603813","#E6E600","#BFBF04","#999903","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#754C24","#f5f2d0")

colind <- integer( length( levels(Idents(marmoset_dataInVivo2)) )  )
for (i in 1:length( levels(Idents(marmoset_dataInVivo2)) ) ) {
  colind[i] <- which(cType==levels(Idents(marmoset_dataInVivo2))[i])
}
coluse <- BaseCol[colind]

#Split the data by time and realign
marmoset_dataInVivo2$Cells <- Idents(marmoset_dataInVivo2)

devStage <- as.character(marmoset_dataInVivo2$Stage)
Loc <- as.character(marmoset_dataInVivo2$LOC)
devStage[grepl("E15C1_",Loc)] <- "CS6_1"
devStage[grepl("E15C2_",Loc)] <- "CS6_2"


Pre_d <- subset(marmoset_dataInVivo2,idents=c("Tb_CS3","Epi_CS3","Hyp_CS3"))
Pre_d$Dat <- "CS3"
CS5_d <- subset(marmoset_dataInVivo2,idents=c("Tb_CS5","ExMes_CS5","EmDisc_CS5","SYS_CS5","Am_CS5","VE_CS5","Am_CS5","PGC_CS5"))
CS5_d$Dat <- "CS5"
CS7_d <- subset(marmoset_dataInVivo2,idents=c("Tb_CS7","ExMes_CS7","EmDisc_CS7","SYS_CS7","Am_CS7","Stalk_CS7","PGC_CS7","VE_CS7"))
CS7_d$Dat <- "CS7"

Idents(marmoset_dataInVivo2) <- devStage
CS6_1 <- subset(marmoset_dataInVivo2,idents=c("CS6_1"))
Idents(CS6_1) <- CS6_1$Cells
CS6_1 <- subset(CS6_1,idents=c("Tb_CS6","ExMes_CS6","EmDisc_CS6","SYS_CS6","Am_CS6","VE_CS6","PGC_CS6","Stalk_CS6"))
CS6_1$Dat <- "CS6_1"

Idents(marmoset_dataInVivo2) <- devStage
CS6_2 <- subset(marmoset_dataInVivo2,idents=c("CS6_2"))
Idents(CS6_2) <- CS6_2$Cells
CS6_2 <- subset(CS6_2,idents=c("Tb_CS6","ExMes_CS6","EmDisc_CS6","SYS_CS6","Am_CS6","VE_CS6","PGC_CS6","Stalk_CS6"))
CS6_2$Dat <- "CS6_2"

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(Pre_d,CS5_d,CS6_1,CS6_2,CS7_d), dims = 1:20, anchor.features = 6000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 20, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
AllGraph1 <- mammal.combined@graphs$integrated_snn
#saveRDS(mammal.combined,file=paste(saveext,"/DimRed/Aligned6k.rds",sep=""))


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dat", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_alignstages6k.pdf",sep=""),width = 32, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dat", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_alignstages6k.pdf",sep=""),width = 32, height = 8)
DimPlot(mammal.combined, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dat", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_alignstages6k.pdf",sep=""),width = 32, height = 8)

#pca <- mammal.combined[["pca"]]
#varExp = (pca@stdev)^2 / pca@misc$total.variance
#0.027226789 0.024295445 0.019411336 0.011525806

#(pca@stdev)^2 / sum((pca@stdev)^2)
#0.15624749 0.13942527 0.11139663
