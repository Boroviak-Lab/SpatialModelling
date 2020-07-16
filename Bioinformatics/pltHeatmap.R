library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")

library(data.table)

set.seed(1)

saveext = "Heatmap/"
dir.create(saveext)

BS<-read.table("../Data/marmKey.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("../Data/marmData.csv",sep=",",header = T, row.names=1)
#marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)

isinvit <- BS$All * BS$QC
labs <- BS$Type
labs2 <- BS$Annotation2
#labs3 <- BS$LABEL_3 
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
#labs3 <- labs3[which(isinvit>0)]

#Load raw counts
#raw_counts3<-read.table("../Data/marmData.csv",sep=",",header = T, row.names=1)
raw_counts3<- raw_counts3[,which(isinvit>0)]
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "marmoset"
marmoset_data$divergence1 <- "newworld"
marmoset_data$divergence2 <- "primate"
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
Idents(marmoset_data) <- labs2
#marmoset_data <- subset(marmoset_data, idents = c("Myo_CS7","ReGland_CS5","ReGland_CS7","Gland_CS5","Gland_CS6","Gland_CS7"), invert = TRUE)
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)

mammal.combined <- marmoset_data
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

#Dimesionality reduction and clustering
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 10) #k.param reduce for later applications? 
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)

#Marker genes
FortyEight <- c("OOEP","TCL1A","WEE2",
                "IL1RN","NOV","ZNF80",
                "SPIC","ESRRB","STAT3",
                "KLF17","SOX15","POU5F1","NANOG","SOX2","SFRP2","DNMT3B","T",
                "TFAP2C","TFAP2A","HOXD3","PRDM1","PRDM14","NANOS3",
                "GATA6","GATA4","PDGFRA","SOX17","OTX2","APOA1","TTR","APOB","HAND2","TBX4","HGF",
                "JAM2","GATA3","GATA2","CGB3","CGA") 

#Force order in heatmap
order <- as.character(Idents(mammal.combined))
order[which(order=="Zy_CS1")] <- "00_Zy_CS1"
order[which(order=="4-cell_CS2")] <- "01_4-cell_CS2"
order[which(order=="8-cell_CS2")] <- "02_8-cell_CS2"
order[which(order=="cMor_CS3")] <- "03_cMor_CS2"
order[which(order=="ICM_CS3")] <- "04_ICM_CS3"
order[which(order=="Epi_CS3")] <- "05_Epi_CS3"
order[which(order=="EmDisc_CS5")] <- "06_EmDisc_CS5"
order[which(order=="EmDisc_CS6")] <- "06_EmDisc_CS6"
order[which(order=="EmDisc_CS7")] <- "06_EmDisc_CS7"
order[which(order=="PGC_CS5")] <- "07_PGC_CS5"
order[which(order=="PGC_CS6")] <- "07_PGC_CS6"

order[which(order=="Am_CS5")] <- "08_Am_CS5"
order[which(order=="Am_CS6")] <- "08_Am_CS6"
order[which(order=="Am_CS7")] <- "08_Am_CS7"

order[which(order=="Hyp_CS3")] <- "09_Hyp_CS3"
order[which(order=="VE_CS5")] <- "10_VE_CS5"
order[which(order=="VE_CS6")] <- "10_VE_CS6"
order[which(order=="SYS_CS5")] <- "11_SYS_CS5"
order[which(order=="SYS_CS6")] <- "11_SYS_CS6"
order[which(order=="SYS_CS7")] <- "11_SYS_CS7"

order[which(order=="ExMes_CS5")] <- "12_ExMes_CS5"
order[which(order=="ExMes_CS6")] <- "12_ExMes_CS6"
order[which(order=="ExMes_CS7")] <- "12_ExMes_CS7"
order[which(order=="ExMes_stalk_CS7")] <- "13_ExMes_stalk_CS7"

order[which(order=="Tb_CS3")] <- "14_Tb_CS3"
order[which(order=="Tb_CS5")] <- "14_Tb_CS5"
order[which(order=="Tb_CS6")] <- "14_Tb_CS6"
order[which(order=="Tb_CS7")] <- "14_Tb_CS7"

order[which(order=="Gland_CS5")] <- "16_Gland_CS5"
order[which(order=="Gland_CS6")] <- "16_Gland_CS6"
order[which(order=="Gland_CS7")] <- "16_Gland_CS7"

order[which(order=="ReGland_CS5")] <- "15_ReGland_CS5"
order[which(order=="ReGland_CS7")] <- "15_ReGland_CS7"
order[which(order=="Myo_CS7")] <- "17_Myo_CS7"

Idents(mammal.combined) <- as.factor(order)

mat_breaks <- seq(0, 2, length.out = 20)

mycolors <- c("#56E600",
"#48BF00",
"#00BF30",
"#009926",
"#00E6E6",
"#00BFBF",
"#0c9cf5",
"#0767DA",
"#0233BF",
"#E6E600",
"#E6E600",
"#877bd6",
"#5F54C7",
"#1A0873",
"#E6B500",
"#F04C04",
"#D74404",
"#E68600",
"#d17600",
"#BF7104",
"#e6c800",
"#c49a00",
"#967700",
"#967700",
"#BF0489",
"#921FE6",
"#8017c2",
"#7108a6",
"#A9A9A9","#A9A9A9",
"#C0C0C0","#C0C0C0","#C0C0C0","#E6E6FA")

annotationL <- c("Zy_CS1",
"4-cell_CS2",
"8-cell_CS2",
"cMor_CS2",
"ICM_CS3",
"Epi_CS3",
"EmDisc_CS5",
"EmDisc_CS6",
"EmDisc_CS7",
"PGC_CS5",
"PGC_CS6",
"Am_CS5",
"Am_CS6",
"Am_CS7",
"Hyp_CS3",
"VE_CS5",
"VE_CS6",
"SYS_CS5",
"SYS_CS6",
"SYS_CS7",
"ExMes_CS5",
"ExMes_CS6",
"ExMes_CS7",
"ExMes_stalk_CS7",
"Tb_CS3",
"Tb_CS5",
"Tb_CS6",
"Tb_CS7",
"ReGland_CS5","ReGland_CS7",
"Gland_CS5","Gland_CS6","Gland_CS7","Myo_CS7")

avexp  <- AverageExpression(object = mammal.combined, slot = "data")

a <- avexp$RNA

annotation_col = data.frame(Stage = factor(colnames(a)))
rownames(annotation_col) <- colnames(a)

names(mycolors) <- colnames(a)
anno_colors <- list(Stage = mycolors)

redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(a[FortyEight,],color =  redblue1(20),gaps_col=c(4,14,20,24,28),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "Heatmap/HM1.pdf")

mat_breaks <- seq(-.8, .8, length.out = 20)
pheatmap(a[FortyEight,],color =  redblue1(20),gaps_col=c(4,14,20,24,28),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = "Heatmap/HM2.pdf")

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a[FortyEight,]+1),color =  redblue1(20),gaps_col=c(4,14,20,24,28),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "Heatmap/M5.pdf")

mat_breaks <- seq(-.8, .8, length.out = 20)
pheatmap(log2(a[FortyEight,]+1),color =  redblue1(20),gaps_col=c(4,14,20,24,28),gaps_row=c(3,6,9,17,23,34,39),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = "Heatmap/HM6.pdf")
