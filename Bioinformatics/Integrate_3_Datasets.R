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

saveext = "CCA_3_dataset/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/Markers/OldList/",sep=""))

dir.create(paste(saveext,"/DimRed/",sep=""))

cType <-c("4-cell_CS2",
          "8-cell_CS2",
          "Am_CS5",
          "Am_CS5_EmDiscExMes",
          "Am_CS5_ExMes",
          "Am_CS5_ExMes_PGC",
          "Am_CS5_ExMesEmDisc",
          "Am_CS5_ExMesPGC",
          "Am_CS5_PGC",
          "Am_CS5_PGC_ExMes",
          "Am_CS6",
          "Am_CS6_EmDisc",
          "Am_CS6_ExMes",
          "Am_CS7",
          "Am_CS7_EmDisc",
          "Am_CS7_ExMes",
          "Am_CS7_ExMes_Tb",
          "Am_CS7_Stalk",
          "Am_CS7_SYS",
          "cMor_CS2",
          "EmDisc_CS5",
          "EmDisc_CS5_Am",
          "EmDisc_CS5_AmExMes",
          "EmDisc_CS5_Gast",
          "EmDisc_CS5_TbExMes",
          "EmDisc_CS5_VE",
          "EmDisc_CS6",
          "EmDisc_CS6_Am",
          "EmDisc_CS6_ExMes",
          "EmDisc_CS6_Gast",
          "EmDisc_CS6_mix",
          "EmDisc_CS6_SYS",
          "EmDisc_CS6_VE",
          "EmDisc_CS6_VE_ExMes",
          "EmDisc_CS7",
          "EmDisc_CS7_Am",
          "EmDisc_CS7_Gast",
          "EmDisc_CS7_PGC",
          "EmDisc_CS7_VE",
          "EmDisc_Gast_CS6",
          "EmDisc_Gast_CS7",
          "Epi_CS3",
          "ExMes_CS5",
          "ExMes_CS5_Am",
          "ExMes_CS5_Am_EmDisc",
          "ExMes_CS5_ReStroma",
          "ExMes_CS5_Tb",
          "ExMes_CS6",
          "ExMes_CS6_Am",
          "ExMes_CS6_Am_EmDisc",
          "ExMes_CS6_EmDisc",
          "ExMes_CS6_EmDisc_Am_VE",
          "ExMes_CS6_SYS",
          "ExMes_CS6_SYS_Tb",
          "ExMes_CS6_Tb",
          "ExMes_CS6_Tb_QC",
          "ExMes_CS7",
          "ExMes_CS7_Am",
          "ExMes_CS7_EmDisc",
          "ExMes_CS7_F",
          "ExMes_CS7_SYS",
          "ExMes_CS7_Tb",
          "Stalk_CS7",
          "F_CS7_",
          "Gland_CS5",
          "Gland_CS5_Stroma",
          "Gland_CS6",
          "Gland_CS7",
          "Gland_CS7_Mix",
          "Hyp_CS3",
          "ICM_CS3",
          "Myo_CS7",
          "NL_CS7_",
          "Oviduct_CS7",
          "PGC_CS5",
          "PGC_CS6",
          "PGC_CS6_ExMes",
          "PGC_CS7",
          "QC_CS6",
          "ReGland_CS5",
          "ReGland_CS7",
          "ReGland_CS7_Tb",
          "ReStroma_CS5",
          "ReStroma_CS6",
          "ReStroma_CS7",
          "Stalk_CS7",
          "Stalk_CS7_Am",
          "Stalk_CS7_PGC",
          "Stroma_CS5",
          "Stroma_CS5_Tb",
          "Stroma_CS6",
          "Stroma_CS7",
          "Stroma_CS7_Tb",
          "SYS_CS5",
          "SYS_CS6",
          "SYS_CS6_EmDisc",
          "SYS_CS6_ExMes",
          "SYS_CS6_VE",
          "SYS_CS7",
          "SYS_CS7_Am",
          "SYS_CS7_ExMes",
          "SYS_CS7_Tb",
          "Tb_abembryonal_CS5",
          "Tb_abembryonal_CS7",
          "Tb_CS3",
          "Tb_CS5",
          "Tb_CS5_Am",
          "Tb_CS5_ExMe",
          "Tb_CS5_ExMes",
          "Tb_CS5_Gland",
          "Tb_CS5_ReStroma",
          "Tb_CS5_Stroma",
          "Tb_CS6",
          "Tb_CS6_Am",
          "Tb_CS6_ExMes",
          "Tb_CS6_ReStroma",
          "Tb_CS6_Stroma",
          "Tb_CS7",
          "Tb_CS7_ExMes",
          "Tb_CS7_Maternal",
          "Tb_CS7_ReStroma",
          "Tb_CS7_Stroma",
          "VE_CS5",
          "VE_CS5_EmDisc",
          "VE_CS5_Tb",
          "VE_CS6",
          "VE_CS6_EmDisc",
          "VE_CS6_EmDisc_ExMes",
          "VE_CS6_ExMes",
          "VE_CS6_ExMes_SYS",
          "VE_CS7",
          "VE_CS7_EmDisc",
          "Zy_CS1",
          "VE_CS4",
          "Epi_CS4",
          "Tb_CS4",
          "Hyp_CS4",
          "EmDisc_CS6/7",
          "ExMes_CS6/7",
          "PGC_CS6/7",
          "SYS_CS6/7",
          "Tb_CS6/7",
          "EmDiscPS_CS6/7",
          "cylPGC",
          "PGC_CS6/7",
          "PGC_CS5",
          "PGC_E50",
          "Ectoderm",
          "Hemogenic Endothelial Progenitors",
          "Endoderm",
          "Advanced Mesoderm",
          "Primitive Streak",
          "YS Mesoderm",
          "Axial Mesoderm",
          "Erythrocytes",
          "Emergent Mesoderm",
          "Epiblast",
          "Nascent Mesoderm",
          "EmDiscPS_CS5",
          "EmDiscPS_CS6",
          "PGC",
          "Stalk_CS6",
          "Other",
          "pPGC",
          "EVTb_CS5",
          "EVTb_CS6",
          "STb_CS5",
          "STb_CS6",
          "STb_CS7")

BaseCol <- c("#48BF00",
             "#00BF30",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#877bd6",
             "#5F54C7",
             "#5F54C7",
             "#5F54C7",
             "#1A0873",
             "#1A0873",
             "#1A0873",
             "#1A0873",
             "#1A0873",
             "#1A0873",
             "#009926",
             "#0c9cf5",
             "#0c9cf5",
             "#0c9cf5",
             "#0c9cf5",
             "#0c9cf5",
             "#0c9cf5",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0767DA",
             "#0233BF",
             "#0233BF",
             "#0233BF",
             "#0233BF",
             "#0233BF",
             "#0767DA",
             "#0233BF",
             "#00BFBF",
             "#e6c800",
             "#e6c800",
             "#e6c800",
             "#e6c800",
             "#e6c800",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#c49a00",
             "#967700",
             "#967700",
             "#967700",
             "#967700",
             "#967700",
             "#967700",
             "#603813",
             "#000000",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#E6B500",
             "#00E6E6",
             "#E6E6FA",
             "#000000",
             "#000000",
             "#E6E600",
             "#BFBF04",
             "#BFBF04",
             "#999903",
             "#000000",
             "#A9A9A9",
             "#A9A9A9",
             "#A9A9A9",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#603813",
             "#603813",
             "#603813",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#C0C0C0",
             "#E68600",
             "#d17600",
             "#d17600",
             "#d17600",
             "#d17600",
             "#BF7104",
             "#BF7104",
             "#BF7104",
             "#BF7104",
             "#921FE6",
             "#7108a6",
             "#BF0489",
             "#921FE6",
             "#921FE6",
             "#921FE6",
             "#921FE6",
             "#921FE6",
             "#921FE6",
             "#921FE6",
             "#8017c2",
             "#8017c2",
             "#8017c2",
             "#8017c2",
             "#8017c2",
             "#7108a6",
             "#7108a6",
             "#7108a6",
             "#7108a6",
             "#7108a6",
             "#F04C04",
             "#F04C04",
             "#F04C04",
             "#D74404",
             "#D74404",
             "#D74404",
             "#D74404",
             "#D74404",
             "#BF3C04",
             "#BF3C04",
             "#56E600",
             "#F04C04",
             "#00BFBF",
             "#BF0489",
             "#E6B500",
             "#0767DA",
             "#c49a00",
             "#BFBF04",
             "#d17600",
             "#8017c2",
             "#0767DA",
             "#999903",
             "#BFBF04",
             "#E6E600",
             "#999903",
             "#1A0873",
             "#4d4d4d",
             "#BF3C04",
             "#662506",
             "#ec7014",
             "#BF7104",
             "#BF7104",
             "#4d4d4d",
             "#993404",
             "#0233BF",
             "#cc4c02",
             "#0c9cf5",
             "#0767DA",
             "#E6E600",
             "#754C24",
             "lightgrey",
             "#E6E600",
             "#dcd0ff",
             "#734f96",
             "#7A4988",
             "#630436",
             "#601A35")

mammal.combined <- readRDS("../Data/Three_dataset.rds")

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 26, height = 8)
mammal.combined$Cells2 <- Idents(mammal.combined)

#Subplots
IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="EmDisc_CS5")] <- "EmDisc_CS5"
IDs[which(mammal.combined$Cells2=="EmDisc_CS6")] <- "EmDisc_CS6"
IDs[which(mammal.combined$Cells2=="EmDisc_CS7")] <- "EmDisc_CS7"
IDs[which(mammal.combined$Cells2=="EmDisc_CS6/7")] <- "EmDisc_CS6/7"
IDs[which(mammal.combined$Cells2=="EmDiscPS_CS5")] <- "EmDiscPS_CS5"
IDs[which(mammal.combined$Cells2=="EmDiscPS_CS6")] <- "EmDiscPS_CS6"
IDs[which(mammal.combined$Cells2=="EmDiscPS_CS6/7")] <- "EmDiscPS_CS6/7"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_embryonalonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_embryonalonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_embryonalonly",".pdf",sep=""),width = 26, height = 8)



IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="Am_CS5")] <- "Am_CS5"
IDs[which(mammal.combined$Cells2=="Am_CS6")] <- "Am_CS6"
IDs[which(mammal.combined$Cells2=="Am_CS7")] <- "Am_CS7"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_amonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_amonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_amonly",".pdf",sep=""),width = 26, height = 8)


IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="ExMes_CS5")] <- "ExMes_CS5"
IDs[which(mammal.combined$Cells2=="ExMes_CS6")] <- "ExMes_CS6"
IDs[which(mammal.combined$Cells2=="ExMes_CS7")] <- "ExMes_CS7"
IDs[which(mammal.combined$Cells2=="ExMes_CS6/7")] <- "ExMes_CS6/7"

Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_exmesonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_exmesonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_exmesonly",".pdf",sep=""),width = 26, height = 8)



IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="Tb_CS5")] <- "Tb_CS5"
IDs[which(mammal.combined$Cells2=="Tb_CS6")] <- "Tb_CS6"
IDs[which(mammal.combined$Cells2=="Tb_CS7")] <- "Tb_CS7"
IDs[which(mammal.combined$Cells2=="Tb_CS6/7")] <- "Tb_CS6/7"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_tbonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_tbonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_tbonly",".pdf",sep=""),width = 26, height = 8)



IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="SYS_CS5")] <- "SYS_CS5"
IDs[which(mammal.combined$Cells2=="SYS_CS6")] <- "SYS_CS6"
IDs[which(mammal.combined$Cells2=="SYS_CS7")] <- "SYS_CS7"
IDs[which(mammal.combined$Cells2=="SYS_CS6/7")] <- "SYS_CS6/7"
IDs[which(mammal.combined$Cells2=="VE_CS5")] <- "VE_CS5"
IDs[which(mammal.combined$Cells2=="VE_CS6")] <- "VE_CS6"
IDs[which(mammal.combined$Cells2=="VE_CS7")] <- "VE_CS7"
Idents(mammal.combined) <- IDs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_sysveonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_sysveonly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_sysveonly",".pdf",sep=""),width = 26, height = 8)

#Back to full anotations
Idents(mammal.combined) <- mammal.combined$Cells2

###################################################
#Some plots 
fullmarkerlist <- c("ANPEP","KLF17","DAPP1","CDCA7L","FABP7","KHDC3L","EMP3","FGF4","INA","LDHC","PDE1B","RND1","SOX15",
                    "SOX2","PDZD4","ATP1B2","CA14","CARD10","CLDN11","CXCL12","DHDH","GPC4","JADE1","LRRN1","SPRY4","USP44","ZIC3",
                    "GSG1","GATM","HNF1B","STX7","ADD3","ANKRD1","APOA1","DENND2C","ELOVL7",
                    "SOX17","NODAL","ESAM","FN1","FST","HNF4A","HPN","ITGA2","LAMA4","LHX1","POSTN","TRPA1","GC",
                    "FXYD3","NOTO","SLC30A2","DYRK3","FABP3","HPGD","PLBD1","RMDN2","SLC6A15","SLC15A2","UNC93A",
                    "HAND2","TBX4","HGF","ADGRA2","C1QTNF3","CD226","FOXF1","OLFML2B","PROM1","SPON1","GDF10","WNT5A","SEMA3A",
                    "CGA","PRR9","CA12","FHDC1","MLLT1","NR2F2","PHOSPHO1","RAB20")

DefaultAssay(mammal.combined) <- "RNA"
for (i in 1:length( fullmarkerlist ) ) {
  FeaturePlot(mammal.combined,  reduction = "tsne", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/DimRed/Markers/TSNE_", fullmarkerlist[i],"_RNA.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined,  reduction = "umap", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/DimRed/Markers/umap_", fullmarkerlist[i],"_RNA.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined,  reduction = "pca", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/DimRed/Markers/PCA_", fullmarkerlist[i],"_RNA.pdf",sep=""),width = 26, height = 8)
  
}

DefaultAssay(mammal.combined) <- "integrated"
for (i in 1:length( fullmarkerlist ) ) {
  FeaturePlot(mammal.combined,  reduction = "tsne", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/TSNE_", fullmarkerlist[i],"_Int.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined,  reduction = "umap", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/umap_", fullmarkerlist[i],"_Int.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined,  reduction = "pca", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/PCA_", fullmarkerlist[i],"_Int.pdf",sep=""),width = 26, height = 8)
}
###################################################

#Now do putative filtering of PGCs (based on our criteria - NANOS3+)
#DefaultAssay(mammal.combined) <- "RNA"
#Idents(mammal.combined) <- mammal.combined$Cells2
#emb <- WhichCells(mammal.combined, idents=c("EmDisc_CS5","EmDisc_CS6","EmDiscPS_CS5","EmDiscPS_CS6","Am_CS5","Am_CS6"))
#pPGC <- intersect(intersect(WhichCells(mammal.combined, expression = NANOS3 > 0),colnames(human_dataA1)), emb)
#DefaultAssay(mammal.combined) <- "integrated"

IDs <- as.character(Idents(mammal.combined))
IDs[1:length(IDs)] <- "Other"
IDs[which(mammal.combined$Cells2=="PGC_CS5")] <- "PGC_CS5"
IDs[which(mammal.combined$Cells2=="PGC_CS6")] <- "PGC_CS6"
IDs[which(mammal.combined$Cells2=="PGC_CS7")] <- "PGC_CS7"
IDs[which(mammal.combined$Cells2=="PGC_CS6/7")] <- "PGC_CS6/7"
IDs[which(mammal.combined$Cells2=="PGC_E50")] <- "PGC_E50"
Idents(mammal.combined) <- IDs
#Idents(mammal.combined,cells=pPGC) <- "pPGC"
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_pgconly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_pgconly",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_pgconly",".pdf",sep=""),width = 26, height = 8)


#Now we can pull out Am and PGC 
Idents(mammal.combined) <- mammal.combined$Cells2
DefaultAssay(mammal.combined) <- "integrated"

mammal.combined_sub <- subset(mammal.combined,idents = c("Am_CS5","Am_CS6","Am_CS7","PGC_CS5","PGC_CS6","PGC_CS7","PGC_CS6/7","PGC_E50"),invert = TRUE)
colind <- integer( length( levels(Idents(mammal.combined_sub)) )  )
for (i in 1:length( levels(Idents(mammal.combined_sub)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined_sub))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined_sub, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 4,  label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_noAmPGC",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined_sub, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 4, label.size = 2,  split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_noAmPGC",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined_sub, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 4, label.size = 2, split.by = "DataOrder", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_noAmPGC",".pdf",sep=""),width = 26, height = 8)

#More expression plots
DefaultAssay(mammal.combined_sub) <- "RNA"
for (i in 1:length( fullmarkerlist ) ) {
  FeaturePlot(mammal.combined_sub,  reduction = "tsne", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/TSNE_", fullmarkerlist[i],"_RNA_sub.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined_sub,  reduction = "umap", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/umap_", fullmarkerlist[i],"_RNA_sub.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined_sub,  reduction = "pca", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/PCA_", fullmarkerlist[i],"_RNA_sub.pdf",sep=""),width = 26, height = 8)
}

DefaultAssay(mammal.combined_sub) <- "integrated"
for (i in 1:length( fullmarkerlist ) ) {
  FeaturePlot(mammal.combined_sub,  reduction = "tsne", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/TSNE_", fullmarkerlist[i],"_Int_sub.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined_sub,  reduction = "umap", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/umap_", fullmarkerlist[i],"_Int_sub.pdf",sep=""),width = 26, height = 8)
  FeaturePlot(mammal.combined_sub,  reduction = "pca", features = fullmarkerlist[i], split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
  ggsave(filename=paste(saveext,"/Markers/OldList/PCA_", fullmarkerlist[i],"_Int_sub.pdf",sep=""),width = 26, height = 8)
}


#Now do some DE analysis 
IntG <- rownames(mammal.combined@assays$integrated) #Integrationi features - these are the important features

#Pull out the marmoset data only and find markers
split <- SplitObject(mammal.combined, split.by = "DataOrder")
AllMarkers <- FindAllMarkers(split$`1) Marmoset`, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

dir.create(paste(saveext,"/EpiCS3/",sep=""))

uGenes1 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Epi_CS3")] )
for (i in 1:length(uGenes1)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes1[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/EpiCS3/Markerscatter_",as.character(uGenes1[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/HypCS3/",sep=""))
uGenes2 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Hyp_CS3")] )
for (i in 1:length(uGenes2)) {
  FeaturePlot(mammal.combined,  reduction = "pca", features = as.character(uGenes2[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/HypCS3/Markerscatter_",as.character(uGenes2[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/TbCS3/",sep=""))
uGenes3 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS3")] )
for (i in 1:length(uGenes3)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes3[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/TbCS3/Markerscatter_",as.character(uGenes3[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/EmDisc_CS5/",sep=""))
uGenes4 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="EmDisc_CS5")] )
for (i in 1:length(uGenes4)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes4[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/EmDisc_CS5/Markerscatter_",as.character(uGenes4[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/SYS_CS5/",sep=""))
uGenes5 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="SYS_CS5")] )
for (i in 1:length(uGenes5)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes5[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/SYS_CS5/Markerscatter_",as.character(uGenes5[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/Tb_CS5/",sep=""))
uGenes6 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS5")] )
for (i in 1:length(uGenes6)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes6[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/Tb_CS5/Markerscatter_",as.character(uGenes6[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/ExMes_CS5/",sep=""))
uGenes7 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="ExMes_CS5")] )
for (i in 1:length(uGenes7)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes7[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/ExMes_CS5/Markerscatter_",as.character(uGenes7[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/EmDisc_CS6/",sep=""))
uGenes4 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="EmDisc_CS6")] )
for (i in 1:length(uGenes4)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes4[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/EmDisc_CS6/Markerscatter_",as.character(uGenes4[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/SYS_CS6/",sep=""))
uGenes5 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="SYS_CS6")] )
for (i in 1:length(uGenes5)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes5[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/SYS_CS6/Markerscatter_",as.character(uGenes5[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/Tb_CS6/",sep=""))
uGenes6 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS6")] )
for (i in 1:length(uGenes6)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes6[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/Tb_CS6/Markerscatter_",as.character(uGenes6[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/ExMes_CS7/",sep=""))
uGenes7 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="ExMes_CS7")] )
for (i in 1:length(uGenes7)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes7[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/ExMes_CS7/Markerscatter_",as.character(uGenes7[i]),".pdf",sep=""),width = 35, height = 8)
}


dir.create(paste(saveext,"/EmDisc_CS7/",sep=""))
uGenes4 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="EmDisc_CS7")] )
for (i in 1:length(uGenes4)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes4[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/EmDisc_CS7/Markerscatter_",as.character(uGenes4[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/SYS_CS7/",sep=""))
uGenes5 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="SYS_CS7")] )
for (i in 1:length(uGenes5)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes5[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/SYS_CS5/Markerscatter_",as.character(uGenes5[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/Tb_CS7/",sep=""))
uGenes6 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="Tb_CS7")] )
for (i in 1:length(uGenes6)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes6[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/Tb_CS7/Markerscatter_",as.character(uGenes6[i]),".pdf",sep=""),width = 35, height = 8)
}

dir.create(paste(saveext,"/ExMes_CS7/",sep=""))
uGenes7 <- intersect(IntG, AllMarkers$gene[which(AllMarkers$cluster=="ExMes_CS5")] )
for (i in 1:length(uGenes7)) {
  FeaturePlot(mammal.combined, reduction = "pca", features = as.character(uGenes7[i]), combine=TRUE, cols =  c("lightgrey", "black"), split.by = "DataOrder", pt.size = 8)
  ggsave(filename=paste(saveext,"/ExMes_CS7/Markerscatter_",as.character(uGenes7[i]),".pdf",sep=""),width = 35, height = 8)
}
