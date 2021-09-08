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

#saveext1 = "~/Desktop/Thorsten/FINAL/AllPlatesCCA_redooriginal/"
#saveext2 = "~/Desktop/Thorsten/FINAL/AllPlatesCCA_3_dataset/"

saveext = "./CCA_4_dataset/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
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
          "cMor_CS3",
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
          "STb_CS7",
          "Blood Progenitors",
          "DE(P)", 
          "Erythroblasts",
          "Hemogenic Endothelium",
          "Non-Neural Ectoderm",
          "Hypoblast",
          "Erythro-Myeloid Progenitors",
          "DE(NP)", 
          "Myeloid Progenitors",
          "YS Endoderm")

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
             "#601A35",
             "#f80000",
             "#BF3C04",
             "#1c0000",
             "#7c0000",
             "#1A0873",
             "#E6B500",
             "#ba0000",
             "#D74404",
             "#3e0000",
             "#BF7104")


mammal.combined <- readRDS("../Data/mammal_combined_3D.rds")
Idents(mammal.combined) <- mammal.combined$Cells


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

#And plot these
DimPlot(mammal.combined, reduction = "pca", cols = coluse, shape.by = "cell.orig", pt.size = 2,  label.size = 2,  split.by = "divergence1", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_4Datasets",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "umap",  cols = coluse, shape.by = "cell.orig",pt.size = 2, label.size = 2,  split.by = "divergence1", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_4Datasets",".pdf",sep=""),width = 26, height = 8)
DimPlot(mammal.combined, reduction = "tsne",  cols = coluse, shape.by = "cell.orig", pt.size = 2, label.size = 2, split.by = "divergence1", label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_4Datasets",".pdf",sep=""),width = 26, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.8)
#mammal.combined$Cl0 <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$Cl9
DimPlot(mammal.combined, pt.size = 4, reduction = "tsne", split.by = "divergence1", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/TSNE_Cl10.pdf",sep=""),width = 32, height = 8)
DimPlot(mammal.combined, pt.size = 4, reduction = "umap", split.by = "divergence1", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl10.pdf",sep=""),width = 32, height = 8)
DimPlot(mammal.combined, pt.size = 4, reduction = "pca", split.by = "divergence1", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/PCA_Cl10.pdf",sep=""),width = 32, height = 8)

DefaultAssay(mammal.combined) <- "RNA"
FeaturePlot(mammal.combined,  reduction = "tsne", features = "TFAP2A", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_TFAP2A_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2A", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_TFAP2A_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2A", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/PCA_TFAP2A_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "VTCN1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_VTCN1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "VTCN1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_VTCN1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "VTCN1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/PCA_VTCN1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "TFAP2C", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_TFAP2C_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "TFAP2C", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_TFAP2C_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "TFAP2C", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/PCA_TFAP2C_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "NANOS3", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_NANOS3_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "NANOS3", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_NANOS3_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "NANOS3", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/PCA_NANOS3_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "SOX17", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/TSNE_SOX17_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "SOX17", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/UMAP_SOX17_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SOX17", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/Markers/PCA_SOX17_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "PDGFRA", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/TSNE_PDGFRA_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "PDGFRA", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_PDGFRA_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "PDGFRA", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/PCA_PDGFRA_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "SNAI2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/TSNE_SNAI2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "SNAI2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_SNAI2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "SNAI2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/PCA_SNAI2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "EOMES", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/TSNE_EOMES_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "EOMES", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_EOMES_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "EOMES", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/PCA_EOMES_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "tsne", features = "FOXA2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/TSNE_FOXA2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "FOXA2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_FOXA2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "pca", features = "FOXA2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/PCA_FOXA2_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "WNT8A", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_WNT8A_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "WNT5A", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_WNT5A_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "BMPER", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_BMPER_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "SFRP1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_SFRP1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "SFRP2", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_SFRP3_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "FST", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_FST_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "TDGF1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_TDGF1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "WNT5B", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_WNT5B_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "CDX1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_CDX1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "MIXL1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_MIXL1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "ISL1", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_ISL1_alignstages.pdf",sep=""),width = 32, height = 8)
FeaturePlot(mammal.combined,  reduction = "umap", features = "GATA4", split.by = "DataOrder", cols =  c("lightgrey", "black"), pt.size = 4)
ggsave(filename=paste(saveext,"/DimRed/Markers/UMAP_GATA4_alignstages.pdf",sep=""),width = 32, height = 8)

#Condense down the names
uID <- as.character(mammal.combined$uData)
uID[which(uID=="1) Marmoset")] <- "Marm1"
uID[which(uID=="2) Cynomolgous")] <- "Cyno1"
uID[which(uID=="3) Human (in vitro)")] <- "Human1"
uID[which(uID=="3) Human (in vivo)")] <- "Human2"
mammal.combined$uData <- uID

#Now do scatter plots on subclustered pseudobulks
uID <- paste(mammal.combined$uData,mammal.combined$Cl9,sep="_")
Idents(mammal.combined) <- uID
DefaultAssay(mammal.combined) <- "RNA"
Ae <- AverageExpression(mammal.combined)
Ae <- Ae$RNA
Ae$gene <- rownames(Ae)
SIGNAL<-read.table("../Data/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("../Data/TF.txt",header = F)
TF <- TF$V1

#And now similarities / differences
#PS vs Pluripotent
M_5 <- FindMarkers(mammal.combined, ident.1 = c("Marm1_5"), ident.2 = c("Marm1_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
H2_5 <- FindMarkers(mammal.combined, ident.1 = c("Human2_5"), ident.2 = c("Human2_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

M_7 <- FindMarkers(mammal.combined, ident.1 = c("Marm1_5"), ident.2 = c("Marm1_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
H2_7 <- FindMarkers(mammal.combined, ident.1 = c("Human2_5"), ident.2 = c("Human2_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#PGC vs PlPo
M_10 <- FindMarkers(mammal.combined, ident.1 = c("Marm1_10"), ident.2 = c("Marm1_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
H2_10 <- FindMarkers(mammal.combined, ident.1 = c("Human2_10"), ident.2 = c("Human2_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

#Am vs PlPo
M_1 <- FindMarkers(mammal.combined, ident.1 = c("Marm1_1"), ident.2 = c("Marm1_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))
H2_1 <- FindMarkers(mammal.combined, ident.1 = c("Human2_1"), ident.2 = c("Human2_0"), verbose = FALSE, test.use = "MAST", only.pos = FALSE, logfc.threshold = log(0))

write.csv(as.data.frame(M_1), file=paste(saveext,"/Marm_Cl1_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(M_5), file=paste(saveext,"/Marm_Cl5_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(M_7), file=paste(saveext,"/Marm_Cl7_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(M_10), file=paste(saveext,"/Marm_Cl10_vs_Cl0.csv",sep=""))

write.csv(as.data.frame(H2_1), file=paste(saveext,"/Human_Cl1_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(H2_5), file=paste(saveext,"/Human_Cl5_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(H2_7), file=paste(saveext,"/Human_Cl7_vs_Cl0.csv",sep=""))
write.csv(as.data.frame(H2_10), file=paste(saveext,"/Human_Cl10_vs_Cl0.csv",sep=""))


##################################
#First compare Marmo vs human1, Cl1  == Amnion
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2))],rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M_1),"X"] <- M_1$avg_logFC
Ae5[rownames(H2_1),"Y"] <- H2_1$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M_1)[which(M_1$p_val_adj<0.1)],rownames(H2_1)[which(H2_1$p_val_adj<0.1)])),"Keep"] <- 1
#Expressed > FCPM in one conditoin
Ae5 <- Ae5[which((Ae5$Marm1_0>0.1 | Ae5$Human2_0>0.1 | Ae5$Marm1_1>0.1 | Ae5$Human2_1>0.1)  ) ,]
genes.to.label1 = c( intersect(intersect(rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2))],rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))]),TF))
genes.to.label2 = c( intersect(intersect(rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2))],rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))]),SIGNAL1))
genes.to.label3 = c( intersect(intersect(rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2))],rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))]),SIGNAL2))
genes.to.label4 = c( intersect(intersect(rownames(M_1)[which(abs(M_1$avg_logFC)>log(1.2))],rownames(H2_1)[which(abs(H2_1$avg_logFC)>log(1.2))]),SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
#Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm1_Human2_Cl1_scatter_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed") +  geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm_Human2_Cl1_scatter_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


##################################
#First compare Marmo vs human1, C1l0 == PGC
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2))],rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M_10),"X"] <- M_10$avg_logFC
Ae5[rownames(H2_10),"Y"] <- H2_10$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M_10)[which(M_10$p_val_adj<0.1)],rownames(H2_10)[which(H2_10$p_val_adj<0.1)])),"Keep"] <- 1
#Expressed > FCPM in one conditoin
Ae5 <- Ae5[which((Ae5$Marm1_0>0.1 | Ae5$Human2_0>0.1 | Ae5$Marm1_10>0.1 | Ae5$Human2_10>0.1)  ) ,]
genes.to.label1 = c( intersect(intersect(rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2))],rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))]),TF))
genes.to.label2 = c( intersect(intersect(rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2))],rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))]),SIGNAL1))
genes.to.label3 = c( intersect(intersect(rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2))],rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))]),SIGNAL2))
genes.to.label4 = c( intersect(intersect(rownames(M_10)[which(abs(M_10$avg_logFC)>log(1.2))],rownames(H2_10)[which(abs(H2_10$avg_logFC)>log(1.2))]),SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
#Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm1_Human2_Cl10_scatter_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed") +  geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm_Human2_Cl10_scatter_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()

############################
###NEXT CLUSTER 5 = PS

Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2))],rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M_5),"X"] <- M_5$avg_logFC
Ae5[rownames(H2_5),"Y"] <- H2_5$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M_5)[which(M_5$p_val_adj<0.1)],rownames(H2_5)[which(H2_5$p_val_adj<0.1)])),"Keep"] <- 1
#Expressed > FCPM in one conditoin
Ae5 <- Ae5[which((Ae5$Marm1_0>0.1 | Ae5$Human2_0>0.1 | Ae5$Marm1_5>0.1 | Ae5$Human2_5>0.1)  ) ,]
genes.to.label1 = c( intersect(intersect(rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2))],rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))]),TF))
genes.to.label2 = c( intersect(intersect(rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2))],rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))]),SIGNAL1))
genes.to.label3 = c( intersect(intersect(rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2))],rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))]),SIGNAL2))
genes.to.label4 = c( intersect(intersect(rownames(M_5)[which(abs(M_5$avg_logFC)>log(1.2))],rownames(H2_5)[which(abs(H2_5$avg_logFC)>log(1.2))]),SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
#Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm1_Human2_Cl5_scatter_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed") +  geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "Delta Expression")
ggsave(filename=paste(saveext,"Marm_Human2_Cl5_scatter_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()



############################
###NEXT CLUSTER 7 = PS
Ae5 <- Ae
Ae5$Pval1 <- 1
Ae5$FC2 <- NA
Ae5$Indi <- 0
Ae5$CG <- 0
Ae5[rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2) )],"Indi"] <- 1
Ae5[rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))],"Indi"] <- 2
Ae5[intersect(rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2))],rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))]),"Indi"] <- 3
Ae5$Indi <- as.factor(Ae5$Indi)
Ae5$X <- 0
Ae5$Y <- 0
Ae5[rownames(M_7),"X"] <- M_7$avg_logFC
Ae5[rownames(H2_7),"Y"] <- H2_7$avg_logFC
Ae5$Keep <- 0
Ae5[unique(c(rownames(M_7)[which(M_7$p_val_adj<0.1)],rownames(H2_7)[which(H2_7$p_val_adj<0.1)])),"Keep"] <- 1
#Expressed > FCPM in one conditoin
Ae5 <- Ae5[which((Ae5$Marm1_0>0.1 | Ae5$Human2_0>0.1 | Ae5$Marm1_7>0.1 | Ae5$Human2_7>0.1)  ) ,]
genes.to.label1 = c( intersect(intersect(rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2))],rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))]),TF))
genes.to.label2 = c( intersect(intersect(rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2))],rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))]),SIGNAL1))
genes.to.label3 = c( intersect(intersect(rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2))],rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))]),SIGNAL2))
genes.to.label4 = c( intersect(intersect(rownames(M_7)[which(abs(M_7$avg_logFC)>log(1.2))],rownames(H2_7)[which(abs(H2_7$avg_logFC)>log(1.2))]),SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
#Ae5 <- Ae5[which(Ae5$Keep==1),]
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed")  + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " ", y = " ")
ggsave(filename=paste(saveext,"Marm1_Human2_Cl7_scatter_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(X,Y)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','#777777','#555555','black')) + geom_hline(yintercept=log(1.2),linetype="dashed") + geom_hline(yintercept=-log(1.2),linetype="dashed") +  geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed") #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = " ", y = " ")
ggsave(filename=paste(saveext,"Marm_Human2_Cl7_scatter_Lig.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()