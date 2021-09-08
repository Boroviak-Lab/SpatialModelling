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

saveext = "./EmbryoHeatmaps/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

#Load in the data
marmoset_dataInVivo2 <- readRDS('../Data/InVivo.rds')

#X <- GetAssayData(marmoset_dataInVivo2, assay = "RNA")
#write.csv(as.data.frame(X), file=paste(saveext,"/X.csv",sep=""))
#write.csv(as.data.frame(marmoset_dataInVivo2$LOC), file=paste(saveext,"/XLOC.csv",sep=""))

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

#Colorcode for plotting
cType <- c("Other","Gland_CS5","Gland_CS6","Gland","Stroma_CS5","Stroma","Remodelled_CS5","Remodelled_CS6","Remodelled","Stalk_CS6","2307","2308","Other","Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3","Tb_CS3","Hyp_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","VE_CS5","VE_CS6","VE_CS7","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS5","Tb_CS6","Tb_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS7","PGC_CS5","PGC_CS6","PGC_CS7","Gland_CS5","Gland_CS6","Gland_CS7","ReGland_CS5","ReGland_CS7","Myo","VE_CS4","Epi_CS4","Tb_CS4","Hyp_CS4","EmDisc_CS6/7","ExMes_CS6/7","PGC_CS6/7","SYS_CS6/7","Tb_CS6/7","Am","Amnoid_bead","BMP_MEF","BMP_noMEF","EmD","EmDisc","ActA_MEF","ActA_noMEF","SB43_MEF","CHIR_MEF","FGF_noMEF","Am_CS5_PGC","Am_CS5_ExMesPGC","EmDisc_CS5_Am","PGC_CS6","Stalk_CS7_PGC","EmDisc_CS7_PGC","EmDisc_CS7_Am","Am_CS6_EmDisc","EmDisc_CS5_Gast","EmDisc_CS6_Am","EmDisc_CS6_Gast","EmDisc_CS7_Gast","Stalk_CS7_PGC","Gland_CS6_","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","VE_CS6","VE_CS7","Am_CS7_EmDisc","EmDisc_Gast_CS6","EmDisc_Gast_CS7","Other")
BaseCol <- c("#f5f2d0","#B3B2B2","#A5A4A3","#969593","#DFDFDF","#D0D0D0","#878684","#797775","#6A6866","#754C24","lightgrey","lightgrey","lightgrey","#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF","#BF0489","#E6B500","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#F04C04","#D74404","#BF3C04","#E68600","#d17600","#BF7104","#921FE6","#8017c2","#7108a6","#e6c800","#c49a00","#967700","#603813","#E6E600","#BFBF04","#999903","#C0C0C0","#C0C0C0","#C0C0C0","#A9A9A9","#A9A9A9","#E6E6FA","#F04C04","#00BFBF","#BF0489","#E6B500","#0767DA","#c49a00","#E6E600","#d17600","#8017c2","#7b3294","#c2a5cf","#a6dba0","#008837","#ca0020","#f4a582","#92c5de","#0571b0","#d8b365","#5ab4ac","#4d4d4d","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#E6E600","#0233BF","#5F54C7","#0c9cf5","#0767DA","#0767DA","#0233BF","#E6E600","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0",'#D74404',"black","#5F54C7","#0767DA","#0233BF","#754C24","#f5f2d0")

#Idents(marmoset_dataInVivo2) <- marmoset_dataInVivo2marmoset_dataInVivo2$Cells

colind <- integer( length( levels(Idents(marmoset_dataInVivo2)) )  )
for (i in 1:length( levels(Idents(marmoset_dataInVivo2)) ) ) {
  colind[i] <- which(cType==levels(Idents(marmoset_dataInVivo2))[i])
}
coluse <- BaseCol[colind]

#PCA plots
DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "pca", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_All.pdf",sep=""),width = 20, height = 16)

#Nonlinear
marmoset_dataInVivo2 <- RunUMAP(marmoset_dataInVivo2, reduction = "pca", dims = 1:20, n.neighbors = 20)
DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "umap", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_All.pdf",sep=""),width = 20, height = 16)

marmoset_dataInVivo2 <- RunTSNE(marmoset_dataInVivo2, reduction = "pca", dims = 1:20, perplexity = 40)
DimPlot(marmoset_dataInVivo2, cols = coluse, pt.size = 4, reduction = "tsne", split.by = "Dataset", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_All.pdf",sep=""),width = 20, height = 16)

#Plot some example markers
FeaturePlot(marmoset_dataInVivo2,  reduction = "tsne", features = "TFAP2A", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_TFAP2A_alignstages.pdf",sep=""),width = 20, height = 16)
FeaturePlot(marmoset_dataInVivo2,  reduction = "tsne", features = "PDGFRA", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_PDGFRA_alignstages.pdf",sep=""),width = 20, height = 16)
FeaturePlot(marmoset_dataInVivo2,  reduction = "tsne", features = "POU5F1", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_POU5F1_alignstages.pdf",sep=""),width = 20, height = 16)
FeaturePlot(marmoset_dataInVivo2,  reduction = "tsne", features = "ISL1", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_ISL1_alignstages.pdf",sep=""),width = 20, height = 16)
FeaturePlot(marmoset_dataInVivo2,  reduction = "tsne", features = "TFAP2C", split.by = "Dataset", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_TFAP2C_alignstages.pdf",sep=""),width = 20, height = 16)

#Generate some volcano plots
VlnPlot(marmoset_dataInVivo2,features = "TFAP2A",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_TFAP2C_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "TFAP2C",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_TFAP2C_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "POU5F1",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_POU5F1_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "SOX2",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_SOX2_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "ISL1",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_ISL1_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "PDGFRA",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_PDGFRA_alignstages.pdf",sep=""),width = 20, height = 16)
VlnPlot(marmoset_dataInVivo2,features = "VTCN1",pt.size=4,idents=c("Epi_CS3","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","Am_CS5","Am_CS6","Am_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7"),cols = c("#00BFBF","#0c9cf5","#0767DA","#0233BF","#877bd6","#5F54C7","#1A0873","#e6c800","#c49a00","#967700"))
ggsave(filename=paste(saveext,"/Markers/Violin_VTCN1_alignstages.pdf",sep=""),width = 20, height = 16)


#Write out average expression values for cell types (CPM)
#AvExp <- AverageExpression(marmoset_dataInVivo2)
#write.csv(as.data.frame(AvExp$RNA *100), file=paste(saveext,"/AverageExp_CPM.csv",sep=""))

#Remove maternal cells and calculate average expression
marmoset_dataInVivo3 <- subset(marmoset_dataInVivo2,idents=c("Stroma","Remodelled","Gland","Myo"),invert=TRUE)
AvExp <- AverageExpression(marmoset_dataInVivo3)

#Marker genes (gene panel)
FortyEight <- c("OOEP","TCL1A","WEE2",
                "IL1RN","NOV","ZNF80",
                "SPIC","ESRRB","STAT3",
                "KLF17","SOX15","POU5F1","NANOG","SOX2","SFRP2","DNMT3B","T",
                "TFAP2C","TFAP2A","VTCN1","PRDM1","PRDM14","NANOS3",
                "GATA6","GATA4","PDGFRA","SOX17","OTX2","APOA1","TTR","APOB","HAND2","TBX4","HGF",
                "JAM2","GATA3","GATA2","CGB3","CGA") 


#Colours for heatmap plots
mycolors <- c("#56E600","#48BF00","#00BF30","#009926","#00E6E6","#00BFBF",
              "#0c9cf5","#0767DA","#0233BF",
              "#E6E600","#E6E600","#E6E600",
              "#877bd6","#5F54C7","#1A0873",
              "#E6B500",
              "#F04C04","#D74404","#BF3C04",
              "#E68600","#d17600","#BF7104",
              "#e6c800","#c49a00","#967700","#c49a00","#967700",
              "#BF0489",
              "#921FE6", "#8017c2","#7108a6")

annotationL <- c("Zy_CS1","4-cell_CS2","8-cell_CS2","cMor_CS2","ICM_CS3","Epi_CS3",
                 "EmDisc_CS5","EmDisc_CS6","EmDisc_CS7",
                 "PGC_CS5","PGC_CS6","PGC_CS7",
                 "Am_CS5","Am_CS6","Am_CS7",
                 "Hyp_CS3",
                 "VE_CS5","VE_CS6","VE_CS7",
                 "SYS_CS5","SYS_CS6","SYS_CS7",
                 "ExMes_CS5","ExMes_CS6","ExMes_CS7","Stalk_CS6","Stalk_CS7",
                 "Tb_CS3",
                 "Tb_CS5","Tb_CS6","Tb_CS7")

exps <- annotationL
mat_breaks <- seq(0, 2, length.out = 20)

#Get genes and samples to plot
X <- (AvExp$RNA[FortyEight,exps])

annotation_col = data.frame(Stage = factor(colnames(X)))
rownames(annotation_col) <- colnames(X)
names(mycolors) <- colnames(X)
anno_colors <- list(Stage = mycolors)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

#Normalised/scaled expression values 
mat_breaks <- seq(0, 100, length.out = 20)
pheatmap(X*100,color =  redblue1(20),gaps_col=c(5,9,12,15,19,22,27), gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/Fullheamtap",".pdf",sep="") ,width=10,height=10)
#Plot the row normalised gene expression values
mat_breaks <- seq(-2, 2, length.out = 20)
pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(5,9,12,15,19,22,27),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/Fullheamtapscale",".pdf",sep=""),width=10,height=10 )
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(5,9,12,15,19,22,27),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/Fullheamtapscale_2",".pdf",sep=""),width=10,height=10 )

#Addtional panel of markers
ReviewerMarkers <- c(
  "FGF1", "FGF2", "FGF3","FGF4", "FGF5", "FGF6", "FGF7","FGF8","FGF9","FGF10","FGF11","FGF12","FGF13", "FGF14", "FGF16", "FGF17", "FGF18", "FGF19", "FGF20", "FGF21", "FGF22", "FGF23",
  "IGF1", "IGF1R",
  "WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10A","WNT10B","WNT11","WNT16",
  "RSPO1","RSPO2","RSPO3","RSPO4",
  "CER1", "SFRP1","SFRP2","SFRP4","SFRP5", "GSK3B", "CTNNB1", "AXIN2", "CDX1", "CDX2", "LEF1")

X <- (AvExp$RNA[ReviewerMarkers,exps])
annotation_col = data.frame(Stage = factor(colnames(X)))
rownames(annotation_col) <- colnames(X)
names(mycolors) <- colnames(X)
anno_colors <- list(Stage = mycolors)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

mat_breaks <- seq(0, 50, length.out = 20)
pheatmap(X*100,color =  redblue1(20),gaps_col=c(5,9,12,15,19,22,27),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE, cluster_cols=FALSE,  filename = paste(saveext,"/Fullheamtap_reviewer",".pdf",sep="") ,width=10,height=10)
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(X+1),color =redblue1(20),gaps_col=c(5,9,12,15,19,22,27),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=TRUE, cluster_cols=FALSE, scale = "row",  filename = paste(saveext,"/Fullheamtapscale_reviewer_cl",".pdf",sep=""),width=10,height=10 )
