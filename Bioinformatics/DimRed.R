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
saveext = "~/Desktop/Thorsten/FINAL/Final_AllGoodShots_wCS6rrrrrr/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))
dir.create(paste(saveext,"/Docs/",sep=""))
dir.create(paste(saveext,"/Monocle/",sep=""))
dir.create(paste(saveext,"/GO/",sep=""))

#Load marmoset data key
BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/KeyCorrect_CPAll2.csv",sep=",",header = T, row.names=1)
TF<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Leaving_package/Dimensionality\ reduction\ techniques\ smart-seq2\ -\ Boroviaklab\ data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)

#Load the marmoset data ... to get ordering
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)

#Get logical indexing for different subgroups. In this case we look for a column E25All
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
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20, perplexity = 30)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 10) #k.param reduce for later applications? 
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)

write.csv(as.data.frame( GetAssayData(object = mammal.combined, slot = 'scale.data') ), file=paste(saveext,"/NormData.csv",sep=""))
write.csv(as.data.frame( GetAssayData(object = mammal.combined, slot = 'data') ), file=paste(saveext,"/RawData.csv",sep=""))

write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsCl.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNE.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddings.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCA.csv",sep=""))



#Plot dimensionality reduction
DimPlot(mammal.combined, reduction = "pca", label = TRUE, dim.1 = 1, dim.2 = 2, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_base",".pdf",sep=""),width = 6, height = 4)

FeatureScatter(object = mammal.combined, feature1 = "PC_1", feature2 = "PC_2")
ggsave(filename=paste(saveext,"/DimRed/PCA_1_2",".pdf",sep=""),width = 6, height = 4)
FeatureScatter(object = mammal.combined, feature1 = "PC_1", feature2 = "PC_3")
ggsave(filename=paste(saveext,"/DimRed/PCA_1_3",".pdf",sep=""),width = 6, height = 4)
FeatureScatter(object = mammal.combined, feature1 = "PC_2", feature2 = "PC_3")
ggsave(filename=paste(saveext,"/DimRed/PCA_2_3",".pdf",sep=""),width = 6, height = 4)

ElbowPlot(mammal.combined,ndims = 30)
ggsave(filename=paste(saveext,"/DimRed/Components",".pdf",sep=""),width = 4, height = 4)


DimPlot(mammal.combined, reduction = "umap", label = TRUE, dim.1 = 10, dim.2 = 10, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP",".pdf",sep=""),width = 6, height = 4)
DimPlot(mammal.combined, reduction = "tsne", label = TRUE, dim.1 = 10, dim.2 = 10)
ggsave(filename=paste(saveext,"/DimRed/TSNE",".pdf",sep=""),width = 6, height = 4)

VizDimLoadings(mammal.combined, dims = 1:6, nfeatures = 30, col = "blue",reduction = "pca", projected = FALSE, balanced = FALSE,ncol = NULL, combine = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_loadings",".pdf",sep=""),width = 8, height = 20)

#Identify markers based on cluster
mammal.combined.markerscl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 100,limitsize = FALSE)
write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
#top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMap_100",".pdf",sep=""),width = 5, height = 40)
TFstrue <- merge(x = as.data.frame(mammal.combined.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/Markers.csv",sep=""))

#dbs <- listEnrichrDbs()
#dbs <- c("GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b","KEGG_2016","WikiPathways_2015","Single_Gene_Perturbations_from_GEO_up","Kinase_Perturbations_from_GEO_up","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
#for (i in 1:length(unique(mammal.combined.markerscl$cluster))) {
#enriched <- enrichr(rownames(mammal.combined.markerscl)[which( (mammal.combined.markerscl$cluster==unique(mammal.combined.markerscl$cluster)[i]) & (mammal.combined.markerscl$p_val_adj<0.05) )], dbs)
#write.csv(as.data.frame(enriched[["GO_Biological_Process_2017b"]]), file=paste(saveext,"/GO/GO_Biological_Process_2017b_", i,sep=""))
#write.csv(as.data.frame(enriched[["GO_Cellular_Component_2017b"]]), file=paste(saveext,"/GO/GO_Cellular_Component_2017b_", i,sep=""))
#write.csv(as.data.frame(enriched[["GO_Molecular_Function_2017b"]]), file=paste(saveext,"/GO/GO_Molecular_Function_2017b_", i,sep=""))
#write.csv(as.data.frame(enriched[["WikiPathways_2015"]]), file=paste(saveext,"/GO/WikiPathways_2015_", i,sep=""))
#write.csv(as.data.frame(enriched[["KEGG_2016"]]), file=paste(saveext,"/GO/KEGG_2016_", i,sep=""))
#write.csv(as.data.frame(enriched[["Single_Gene_Perturbations_from_GEO_up"]]), file=paste(saveext,"/GO/Single_Gene_Perturbations_from_GEO_up_", i,sep=""))
#write.csv(as.data.frame(enriched[["Kinase_Perturbations_from_GEO_up"]]), file=paste(saveext,"/GO/Kinase_Perturbations_from_GEO_up_", i,sep=""))
#write.csv(as.data.frame(enriched[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]), file=paste(saveext,"/GO/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_", i,sep=""))
#}

#Specific cluster comparison ...
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)

#Take top 20 for plotting
top20 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/Top20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top20$cluster))) {
  genes<-top20$gene[which(top20$cluster==unique(top20$cluster)[i])]
  VlnPlot(mammal.combined, features = genes)
  ggsave(filename=paste(saveext,"/Markers/Marker_",i,".pdf",sep=""),width = 24, height = 24)
  FeaturePlot(mammal.combined, features = genes, combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
  ggsave(filename=paste(saveext,"/Markers/Markerscatter_",i,".pdf",sep=""),width = 24, height = 24)
}









#marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
#DimPlot(marrow)

#FeaturePlot(mammal.combined, features = top20$gene[which(top20$cluster==unique(top20$cluster)[1])], split.by = "species")
#ggsave(filename=paste(saveext,"UMAP_FeatureCluster1_Split",".pdf",sep=""),width = 24, height = 24)
#dev.off()

#Can do the same for violin plot representations
#VlnPlot(mammal.combined, features = mammal.combined.markerscl$gene)
#ggsave(filename=paste(saveext,"Cluster_Split",".pdf",sep=""),width = 24, height = 24)
#dev.off()

#Manually compare clusters to get markers out for specific clusters rather than 1-vs-all in previous
#mammal.combined.specificmarkers1_6 <- FindMarkers(mammal.combined, ident.1 = "1", ident.2 = "6", only.pos = TRUE)
#mammal.combined.specificmarkers4_6 <- FindMarkers(mammal.combined, ident.1 = "4", ident.2 = "6", only.pos = TRUE)
#mammal.combined.specificmarkers0_6 <- FindMarkers(mammal.combined, ident.1 = "0", ident.2 = "6", only.pos = TRUE)
#write.csv(as.data.frame(mammal.combined.specificmarkers1_6), file=paste(saveext,"Markers_cl_1_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers4_6), file=paste(saveext,"Markers_cl_4_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers0_6), file=paste(saveext,"Markers_cl_0_6.csv",sep=""))

#...
#mammal.combined.specificmarkers6_1 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "1", only.pos = TRUE)
#mammal.combined.specificmarkers6_4 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "4", only.pos = TRUE)
#mammal.combined.specificmarkers6_0 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "0", only.pos = TRUE)
#write.csv(as.data.frame(mammal.combined.specificmarkers1_6), file=paste(saveext,"Markers_cl_1_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers4_6), file=paste(saveext,"Markers_cl_4_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers0_6), file=paste(saveext,"Markers_cl_0_6.csv",sep=""))

#Assign cluster based on annotated cell type (read in from a seperate file)
Cls <- Idents(mammal.combined)
Idents(mammal.combined) <- labs2


mammal.combined.markersncl <- FindMarkers(mammal.combined, ident.1 = c("Gland_CS5","Gland_CS6","Gland_CS7","Myo_CS7","ReGland_CS5","ReGland_CS7") , ident.2 = c("Am_CS5","Am_CS6","Am_CS7","EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","ExMes_CS5","ExMes_CS6","ExMes_CS7","ExMes_stalk_CS7","PGC_CS5","PGC_CS6","SYS_CS5","SYS_CS6","SYS_CS7","Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","VE_CS5","VE_CS6"),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(as.data.frame(mammal.combined.markersncl), file=paste(saveext,"/Maternal.csv",sep=""))



mammal.combined.markersncl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Markers <- c("OOEP", "DNMT1", "PPP1R3A", "WEE2", "TCL1A", "POU5F1", "NANOG", "SOX2", "SPIC", "NOTO", "MIOX", "WWTR1","KLF17", "ANPEP", "LEFTY2", "KHDC3L", "DPPA5", "NLRP7", "DPPA3","DNMT3B", "TDGF1", "SFRP1", "IGFBP2", "SFRP2","APOA1", "RSPO3", "GSG1", "GATM", "OTC", "SPARC", "FSCN1", "CFL1","NODAL", "POSTN", "LHX1", "TRPA1", "FZD5","TTR", "APOB", "SERPINA1", "TF", "CYTL1", "SPINK1", "GATA2", "CGA", "CGB3", "PIEZO1", "PIEZO2", "MLLT1", "CA12")
for (i in 1:length(Markers) ) {
  FeaturePlot(mammal.combined, reduction = "tsne", features = as.character(Markers[i]), combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
  ggsave(filename=paste(saveext,"/Markers/Markerscatter_tsne_",Markers[i],".pdf",sep=""),width = 5, height = 4)
  
}

#for (i in 1:length(unique(mammal.combined.markersncl$cluster))) {
#  enriched <- enrichr(rownames(mammal.combined.markersncl)[which( (mammal.combined.markersncl$cluster==unique(mammal.combined.markersncl$cluster)[i]) & (mammal.combined.markerscl$p_val_adj<0.05) )], dbs)
#  write.csv(as.data.frame(enriched[["GO_Biological_Process_2015"]]), file=paste(saveext,"/GO/GO_Biological_Process_2015_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["GO_Cellular_Component_2015"]]), file=paste(saveext,"/GO/GO_Cellular_Component_2015_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["GO_Molecular_Function_2015"]]), file=paste(saveext,"/GO/GO_Molecular_Function_2015_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["WikiPathways_2015"]]), file=paste(saveext,"/GO/WikiPathways_2015_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["KEGG_2016"]]), file=paste(saveext,"/GO/KEGG_2016_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["Single_Gene_Perturbations_from_GEO_up"]]), file=paste(saveext,"/GO/Single_Gene_Perturbations_from_GEO_up_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["Kinase_Perturbations_from_GEO_up"]]), file=paste(saveext,"/GO/Kinase_Perturbations_from_GEO_up_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#  write.csv(as.data.frame(enriched[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]), file=paste(saveext,"/GO/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_", unique(mammal.combined.markersncl$cluster)[i],sep=""))
#}

write.csv(as.data.frame(table(labs2)), file=paste(saveext,"NumberofCells.csv",sep=""))
cellcount <- as.data.frame(table(labs2))
#p<-ggplot(data=cellcount, aes(x=labs2, y=Freq)) + geom_bar(stat="Identity") + theme_classic()
#ggsave(filename=paste(saveext,"/DimRed/NumberofCells",".pdf",sep=""),width = 60, height = 10,limitsize = FALSE)
cellcount <- cellcount[which(cellcount$Freq>0),]
p<-ggplot(data=cellcount, aes(x=labs2, y=Freq)) + geom_bar(stat="Identity") + theme_classic()
ggsave(filename=paste(saveext,"/DimRed/NumberofCells",".pdf",sep=""),width = 20, height = 10,limitsize = FALSE)


mammal.combined.markersncl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
top100 <- mammal.combined.markersncl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapncl_100",".pdf",sep=""),width = 20, height = 120,limitsize = FALSE)


#Take top 20 for plotting
top20 <- mammal.combined.markersncl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/Top20MarkersNCL.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
genes<-unique(top20$gene)
for (i in 1:length(genes)) {
  VlnPlot(mammal.combined, features = genes[i])
  ggsave(filename=paste(saveext,"/Markers/Marker_",genes[i],".pdf",sep=""),width = 24, height = 24)
  FeaturePlot(mammal.combined, features = genes[i], reduction="tsne", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
  ggsave(filename=paste(saveext,"/Markers/Markerscatter_tsne_",genes[i],".pdf",sep=""),width = 24, height = 24)
}

#PIEZO, IFITM5
#"IFITM5"
Markers <- c("OOEP", "DNMT1", "PPP1R3A", "WEE2", "TCL1A", "POU5F1", "NANOG", "SOX2", "SPIC", "NOTO", "MIOX", "WWTR1","KLF17", "ANPEP", "LEFTY2", "KHDC3L", "DPPA5", "NLRP7", "DPPA3","DNMT3B", "TDGF1", "SFRP1", "IGFBP2", "SFRP2","APOA1", "RSPO3", "GSG1", "GATM", "OTC", "SPARC", "FSCN1", "CFL1","NODAL", "POSTN", "LHX1", "TRPA1", "FZD5","TTR", "APOB", "SERPINA1", "TF", "CYTL1", "SPINK1", "GATA2", "CGA", "CGB3", "PIEZO1", "PIEZO2", "MLLT1", "CA12")
for (i in 1:length(Markers) ) {
  FeaturePlot(mammal.combined, reduction = "tsne", features = as.character(Markers[i]), combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
  ggsave(filename=paste(saveext,"/Markers/Markerscatter_tsne_",Markers[i],".pdf",sep=""),width = 5, height = 4)
  
}



FeaturePlot(mammal.combined, reduction = "tsne", features = "TFAP2C", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste(saveext,"/Markers/Markerscatter_tsne_TFAP2C",sep=""),width = 5, height = 4)

#FeaturePlot(mammal.combined, features = c("OOEP"), combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 3)
#ggsave(filename=paste(saveext,"/Markers/Markerscatter_Z.pdf",sep=""),width = 5, height = 4)



#Label-type plot
#Zygote/4-cell: OOEP, DNMT1, PPP1R3A, WEE2, TCL1A,
#Core pluripotency: POU5F1, NANOG, SOX2   ---    Early ICM: SPIC, NOTO, MIOX, WWTR1,
#Preimplantation EPI: KLF17, ANPEP, LEFTY2, KHDC3L, DPPA5, IFITM5, NLRP7, DPPA3
#Postimplantation EPI: DNMT3B, TDGF1, SFRP1, IGFBP2, SFRP2
#HYPO: APOA1, RSPO3, GSG1, GATM, OTC, SPARC, FSCN1, CFL1
#VE: NODAL, POSTN, LHX1, TRPA1, FZD5
#SYS: TTR, APOB, SERPINA1, TF, CYTL1, SPINK1
#TB: GATA2, CGA, CGB3, PIEZO, MLLT1, CA12


#mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:30)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:30, k.param = 10) #k.param reduce for later applications? 
##mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:30)
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)


pdf('00_qPCR.pdf', height=10, width=10)(heatmap.2(log2(mat+1),(          Rowv=NA, (          # Colv=NA, (          col =  ePalFunc(100), (          scale="row", (          # margins=c(2,2), (          dendrogram = "col",(          trace = "none"(          # labCol = c("ICM", "PreEPI", "PostEPI_early", "PostEPI_late")())

mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20, perplexity = 40)


FeaturePlot(mammal.combined, reduction = "tsne", features = "OTX2", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_tsne_OTX2.pdf",sep=""),width = 5, height = 4)
FeaturePlot(mammal.combined, reduction = "tsne", features = "GATA6", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_tsne_GATA6.pdf",sep=""),width = 5, height = 4)
FeaturePlot(mammal.combined, reduction = "tsne", features = "GATA2", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_tsne_GATA2.pdf",sep=""),width = 5, height = 4)
FeaturePlot(mammal.combined, reduction = "tsne", features = "KRT7", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_tsne_KRT7.pdf",sep=""),width = 5, height = 4)
FeaturePlot(mammal.combined, reduction = "tsne", features = "SOX17", combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
ggsave(filename=paste("~/Desktop/Markerscatter_tsne_SOX17.pdf",sep=""),width = 5, height = 4)


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
DimPlot(mammal.combined, reduction = "umap",  cols = coluse,label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)
DimPlot(mammal.combined, reduction = "tsne", cols = coluse, label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)

DimPlot(mammal.combined, reduction = "tsne", cols = coluse, label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_p40",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)


dir.create(paste(saveext,"/Markers2/",sep=""))

DimPlot(mammal.combined, reduction = "tsne", cols = coluse, label.size = 2, no.legend = TRUE, label = FALSE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/Markers2/TSNE_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)

Markers <- c("OOEP", "DNMT1", "PPP1R3A", "WEE2", "TCL1A", "POU5F1", "NANOG", "SOX2", "SPIC", "NOTO", "MIOX", "WWTR1","KLF17", "ANPEP", "LEFTY2", "KHDC3L", "DPPA5", "NLRP7", "DPPA3","DNMT3B", "TDGF1", "SFRP1", "IGFBP2", "SFRP2","APOA1", "RSPO3", "GSG1", "GATM", "OTC", "SPARC", "FSCN1", "CFL1","NODAL", "POSTN", "LHX1", "TRPA1", "FZD5","TTR", "APOB", "SERPINA1", "TF", "CYTL1", "SPINK1", "GATA2", "CGA", "CGB3", "PIEZO1", "PIEZO2", "MLLT1", "CA12")
for (i in 1:length(Markers) ) {
  FeaturePlot(mammal.combined, reduction = "tsne", features = as.character(Markers[i]), combine=TRUE, cols =  c("lightgrey", "black"), pt.size = 2)
  ggsave(filename=paste(saveext,"/Markers2/Markerscatter_tsne_",Markers[i],".pdf",sep=""),width = 5, height = 4)
  
}


DimPlot(mammal.combined, reduction = "umap",  cols = coluse,label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 15, height = 8, useDingbats=FALSE)


#mammal.combined$CellType2 <- labs3
#
#subsetdata <- subset(mammal.combined, idents = c("earlyTb_CS3","Hyp_CS3","Epi_CS3"), invert = FALSE)
#subsetdata <- RunPCA(subsetdata, npcs = 30, verbose = FALSE) #,features = rownames(Markers))
#subsetdata <- RunUMAP(subsetdata, reduction = "pca", dims = 1:30)

#DimPlot(subsetdata, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_subset",".pdf",sep=""),width = 15, height = 8)

#Idents(subsetdata) <- subsetdata$CellType2
#DimPlot(subsetdata, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_subset2",".pdf",sep=""),width = 15, height = 8)


dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b", "GO_Biological_Process_2017b","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","Gene_Perturbations_from_GEO_up","Kinase_Perturbations_from_GEO_up")

AvExpC <- AverageExpression(object = mammal.combined, use.scale = TRUE)
R3 <- as.data.frame(AvExpC$RNA)
for (i in 1:dim(R3)[2]) {
  Genes <- R3[order(R3[,i],decreasing = TRUE), i, drop = FALSE]
  enriched <- enrichr(as.character(rownames(Genes)[1:1000]), dbs)
  #write.csv(as.data.frame(enriched[["GO_Biological_Process_2017b"]]), file=paste(saveext,"GO/GOBP_", colnames(R3)[i],".csv",sep=""))
  #write.csv(as.data.frame(enriched[["GO_Cellular_Component_2017b"]]), file=paste(saveext,"GO/GOCC_", colnames(R3)[i],".csv",sep=""))
  #write.csv(as.data.frame(enriched[["GO_Molecular_Function_2017b"]]), file=paste(saveext,"GO/GOMF_", colnames(R3)[i],".csv",sep=""))
  #write.csv(as.data.frame(enriched[["Kinase_Perturbations_from_GEO_up"]]), file=paste(saveext,"GO/Kinase_", colnames(R3)[i],".csv",sep=""))
  write.csv(as.data.frame(enriched[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]), file=paste(saveext,"GO/Encode_", colnames(R3)[i],".csv",sep=""))
  write.csv(as.data.frame(enriched[["Gene_Perturbations_from_GEO_up"]]), file=paste(saveext,"GO/GenePert_", colnames(R3)[i],".csv",sep=""))
}



FortyEightPanel <- c("POU5F1","NANOG","SOX2","IL1RN","NOV","ZNF80","SPIC","ESRRB","STAT3","SLC6A8","MT3","ATP12A","GATA2","CDX2","LAMA3","GATA3","NPPA","RSPO3","APOA1","GATA4","GATA6","SOX17","NLRP9","BRDT","CD53","KLF17","TFCP2L1","TFAP2C","SOX15","DDX43","MFAP2","SFRP2","FBXO2","WFDC2","COL1A1","SFRP1","WNT5A","T","APLNR","FGF17","LHX1","KLHL14","LGR5","HGF","CREB3L1","HAND2")




write.csv(as.data.frame(AverageExpression(object = mammal.combined, use.counts = TRUE)), file=paste(saveext,"AvExp_CPM.csv",sep=""))
write.csv(as.data.frame(AverageExpression(object = mammal.combined, use.scale = TRUE)), file=paste(saveext,"AvExp_Scaled.csv",sep=""))

FortyEightPanel <- c("")
DoHeatmap(mammal.combined, features = FortyEightPanel) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/FortyEightPanel",".pdf",sep=""),width = 20, height = 20,limitsize = FALSE)


avexpcount <- as.data.frame(AverageExpression(object = mammal.combined, use.counts = TRUE))
avexpcount1 <- transpose(avexpcount)
avexpcount2 <- as.data.frame(rowMeans(avexpcount1))
colnames(avexpcount2) <- "Freq"
avexpcount2$ID <- colnames(avexpcount)
p<-ggplot(data=avexpcount2, aes(x=ID, y=Freq)) + geom_bar(stat="Identity") + theme_classic()
ggsave(filename=paste(saveext,"/DimRed/AverageExperssionCount",".pdf",sep=""),width = 20, height = 10,limitsize = FALSE)


#avexpscaled <- as.data.frame(AverageExpression(object = mammal.combined, use.scale = TRUE))
#avexpscaled <- colMeans(avexpscaled)




#TFstrue <- merge(x = as.data.frame(mammal.combined.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
#write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Merkers/Markers.csv",sep=""))

AvExp <- AverageExpression(object = mammal.combined, use.counts = TRUE)
AvExpDF <- as.data.frame(AvExp, row.names = rownames(AvExp$RNA))
AvExpDFTF <- merge(x = AvExpDF, y=as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(AvExpDFTF), file=paste(saveext,"AvExp_CPM.csv",sep=""))

AvExp <- AverageExpression(object = mammal.combined, use.scale = TRUE)
AvExpDF <- as.data.frame(AvExp, row.names = rownames(AvExp$RNA))
AvExpDFTF <- merge(x = AvExpDF, y=as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(AvExpDFTF), file=paste(saveext,"AvExp_Scaled.csv",sep=""))


clusterZ4C.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Zy_CS1_"), ident.2 = c("4-cell_CS2_"), min.pct = 0.25) #, only.pos = TRUE)
cluster4C8C.markers <- FindMarkers(object = mammal.combined, ident.1 = c("4-cell_CS2_"), ident.2 = c("8-cell_CS2_"), min.pct = 0.25) #, only.pos = TRUE)
cluster8CMor.markers <- FindMarkers(object = mammal.combined, ident.1 = c("8-cell_CS2_"), ident.2 = c("cMor_CS3_"), min.pct = 0.25) #, only.pos = TRUE)
clusterCMoreICM.markers <- FindMarkers(object = mammal.combined, ident.1 = c("cMor_CS3_"), ident.2 = c("ICM_CS3_"), min.pct = 0.25) #, only.pos = TRUE)
clusterICMpreEpi.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ICM_CS3_"), ident.2 = c("Epi_CS3_"), min.pct = 0.25) #, only.pos = TRUE)

#Early differences
clusterICMHyp.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Epi_CS3_"), ident.2 = c("Hyp_CS3_"), min.pct = 0.25) #, only.pos = TRUE)
clusterICMTb.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Epi_CS3_"), ident.2 = c("earlyTb_CS3_"), min.pct = 0.25) #, only.pos = TRUE)
clusterHypTb.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Hyp_CS3_"), ident.2 = c("earlyTb_CS3_"), min.pct = 0.25) #, only.pos = TRUE)

#Embryonic development
clusterpreEpiEmDisc.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Epi_CS3_"), ident.2 = c("EmDisc_CS5_"), min.pct = 0.25) #, only.pos = TRUE)
clusterEmDisc.markers <- FindMarkers(object = mammal.combined, ident.1 = c("EmDisc_CS5_"), ident.2 = c("EmDisc_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#SYS
clusterHypSYS.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Hyp_CS3_"), ident.2 = c("VE_CS5_"), min.pct = 0.25) #, only.pos = TRUE)
clusterSYS.markers <- FindMarkers(object = mammal.combined, ident.1 = c("VE_CS5_"), ident.2 = c("SYS_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#Tb developmet
clusterTbTb1.markers <- FindMarkers(object = mammal.combined, ident.1 = c("earlyTb_CS3_"), ident.2 = c("Tb_CS5_"), min.pct = 0.25) #, only.pos = TRUE)
clusterTbTb2.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Tb_CS5_"), ident.2 = c("Tb_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#ExMes developmet
clusterHypMes.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Hyp_CS3_"), ident.2 = c("ExMes_CS5_"), min.pct = 0.25) #, only.pos = TRUE)
clusterMes.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ExMes_CS5_"), ident.2 = c("ExMes_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#EmDisc vs Amnion
clusterEmDiscAm.markers <- FindMarkers(object = mammal.combined, ident.1 = c("EmDisc_CS7_"), ident.2 = c("Am_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#SYS Amnion
clusterSYSAm.markers <- FindMarkers(object = mammal.combined, ident.1 = c("SYS_CS7_"), ident.2 = c("Am_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#SYS Tb
clusterSYSTb.markers <- FindMarkers(object = mammal.combined, ident.1 = c("SYS_CS7_"), ident.2 = c("Tb_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#SYS EMDisc
clusterSYSEm.markers <- FindMarkers(object = mammal.combined, ident.1 = c("SYS_CS7_"), ident.2 = c("EmDisc_CS7_"), min.pct = 0.25) #, only.pos = TRUE)


#Amn Tb
clusterAmTb.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Am_CS7_"), ident.2 = c("Tb_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#Mes Tb
clusterMesTb.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ExMes_CS7_"), ident.2 = c("Tb_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#Mes Am
clusterMesAm.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ExMes_CS7_"), ident.2 = c("Am_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#Mes Epi
clusterMesEmDisc.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ExMes_CS7_"), ident.2 = c("EmDisc_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#Mes SYS
clusterMesSYS.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ExMes_CS7_"), ident.2 = c("SYS_CS7_"), min.pct = 0.25) #, only.pos = TRUE)

#
clusterMaternal1.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Gland_CS7_"), ident.2 = c("ReGland_CS7_"), min.pct = 0.25) #, only.pos = TRUE)
clusterMaternal2.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Gland_CS7_"), ident.2 = c("Myo_CS7_"), min.pct = 0.25) #, only.pos = TRUE)
clusterMaternal3.markers <- FindMarkers(object = mammal.combined, ident.1 = c("ReGland_CS7_"), ident.2 = c("Myo_CS7_"), min.pct = 0.25) #, only.pos = TRUE)
clusterMaternal4.markers <- FindMarkers(object = mammal.combined, ident.1 = c("Gland_CS5_","Gland_CS7_"), ident.2 = c("Stroma_CS5_","Stroma_CS7_"), min.pct = 0.25) #, only.pos = TRUE)



clusterZ4C.markers$Comp <- "4C_Z"
cluster4C8C.markers$Comp <- "4C_8C"
cluster8CMor.markers$Comp <- "8C_cMor"
clusterCMoreICM.markers$Comp <- "eICM_cMor"
clusterICMpreEpi.markers$Comp <- "preEpi_eICM"

clusterICMHyp.markers$Comp <- "Hyp_eICM"
clusterICMTb.markers$Comp <- "ICM_Tb"
clusterHypTb.markers$Comp <- "Hyp_Tb"

clusterpreEpiEmDisc.markers$Comp <- "EmDiscE15_preEpi"
clusterEmDisc.markers$Comp <- "EmDiscE25_EmDiscE15" 

clusterHypSYS.markers$Comp <- "SYS_Hyp" 
clusterSYS.markers$Comp <- "SYS"

clusterSYSAm.markers$Comp <- "Am_SYS"
clusterSYSTb.markers$Comp <- "Tb_SYS"
clusterAmTb.markers$Comp <- "Tb_Am"
clusterTbTb1.markers$Comp <- "Tb_Tb1"
clusterTbTb2.markers$Comp <- "Tb_Tb2"
clusterHypMes.markers$Comp <- "HypMes"
clusterMes.markers$Comp <- "Mes"
clustereEmDiscAm.markers$Comp <- "EmDisc_Am"
clusterMesTb.markers$Comp <- "Mes_Tb"
clusterMesAm.markers$Comp <- "Mes_Am"
clusterMesEmDisc.markers$Comp <- "Mes_EmDisc" 
clusterMesSYS.markers$Comp <- "Mes_SYS"
clusterSYSEm.markers$Comp <- "EmDisc_SYS"

clusterMaternal1.markers$Comp <- "Gland_ReGland"
clusterMaternal2.markers$Comp <- "Gland_Myo"
clusterMaternal3.markers$Comp <- "ReGland_Myo"
clusterMaternal4.markers$Comp <- "Gland_Stroma"



Names <- c( rownames(clusterZ4C.markers),
          rownames(cluster4C8C.markers),
          rownames(cluster8CMor.markers),
          rownames(clusterCMoreICM.markers),
          rownames(clustereICMpreEpi.markers),
          rownames(clustereICMHyp.markers),
          rownames(clusterpreEpiEmDisc.markers),
          rownames(clusterEmDisc.markers),
          rownames(clusterHypSYS.markers),
          rownames(clusterSYS.markers),
          rownames(clusterSYSAm.markers),
          rownames(clusterSYSTb.markers),
          rownames(clusterAmTb.markers),
          rownames(clusterTbTb1.markers),
          rownames(clusterTbTb2.markers),
          rownames(clusterHypMes.markers),
          rownames(clusterMes.markers),
          rownames(clustereEmDiscAm.markers),
          rownames(clusterMesTb.markers),
          rownames(clusterMesAm.markers),
          rownames(clusterMesEmDisc.markers),
          rownames(clusterMesSYS.markers),
          rownames(clusterSYSEm.markers),
          rownames(clusterMaternal1),
          rownames(clusterMaternal2),
          rownames(clusterMaternal3),
          rownames(clusterMaternal4))





AllDevmarkers <- bind_rows(clusterZ4C.markers,
cluster4C8C.markers,
cluster8CMor.markers,
clusterCMoreICM.markers,
clusterICMpreEpi.markers,
clusterICMHyp.markers,
clusterpreEpiEmDisc.markers,
clusterEmDisc.markers,
clusterHypSYS.markers,
clusterSYS.markers,
clusterSYSAm.markers,
clusterSYSTb.markers,
clusterAmTb.markers,
clusterTbTb1.markers,
clusterTbTb2.markers,
clusterHypMes.markers,
clusterMes.markers,
clusterEmDiscAm.markers,
clusterMesTb.markers,
clusterMesAm.markers,
clusterMesEmDisc.markers,
clusterMesSYS.markers,
clusterSYSEm.markers,
clusterMaternal1.markers,
clusterMaternal2.markers,
clusterMaternal3.markers,
clusterMaternal4.markers
)

AllDevmarkers$Names <- Names
AllDevmarkers2 <- AllDevmarkers[which(AllDevmarkers$p_val_adj<0.05),]

write.csv(as.data.frame(AllDevmarkers2), file=paste(saveext,"/Markers/AllDevmarkers.csv",sep=""))

Idents(mammal.combined) <- Cls
#Write everythign out to csv
write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed2/EmbeddingsClr.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed2/EmbeddingsTSNEr.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed2/Embeddingsr.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed2/EmbeddingsKey.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed2/EmbeddingsPCAr.csv",sep=""))

#######Now redo with the dynamic gene  lists
#DevMarkersR <-read.table("/Users/christopherpenfold/Desktop/Thorsten/FINAL/Grr/untitled\ folder/ALL_GS_Modelling20k/Markers/AllDevmarkers.csv",sep=",",header = T, row.names=1)
#DevMarkers <- unique(DevMarkersR$Names)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
isinvit <- BS$AllGSNoInVit #$GS #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$GS 
labs <- BS$Type #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$Type
labs2 <- BS$Annotation #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$LABEL_3
labs3 <- BS$LABEL_4 #[order(match(rownames(BS)  , names(Idents(marmoset_data)) )),]$Location
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
labs3 <- labs3[which(isinvit>0)]
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
raw_counts3<- raw_counts3[,which(isinvit>0)]
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "marmoset"
marmoset_data$divergence1 <- "newworld"
marmoset_data$divergence2 <- "primate"
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)
mammal.combined <- marmoset_data
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)

#Dimesionality reduction and clustering
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE, features = unique(AllDevmarkers2$Names) )
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20, k.param = 10) #k.param reduce for later applications? 
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)

Cls <- Idents(mammal.combined)
Idents(mammal.combined) <- labs2

DimPlot(mammal.combined, reduction = "pca", label = TRUE, dim.1 = 1, dim.2 = 2, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_base_Type_DMUsingMarkers",".pdf",sep=""),width = 15, height = 8)
DimPlot(mammal.combined, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_DMUsingMarkers",".pdf",sep=""),width = 15, height = 8)
DimPlot(mammal.combined, reduction = "tsne", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_DMUsingMarkers",".pdf",sep=""),width = 15, height = 8)


#write.csv(as.data.frame( GetAssayData(object = mammal.combined, slot = 'scale.data') ), file=paste(saveext,"/NormData.csv",sep=""))

ND <- as.data.frame( GetAssayData(object = mammal.combined, slot = 'data') )
#mtch <- match(as.character(FortyEightPanel), rownames(ND))
mtch <- match(as.character(FortyEightPanel), rownames(ND))
mtch <- mtch[is.na(mtch)==0]
D <- ND[mtch,isinvit]
colnames(D) <- labs2
pdf(paste(saveext,"/Exp",i ,"_8Cl.pdf",sep=""),width=10,height=siz[count]) #,width=2200, height=2000)
pheatmap(D,cluster_cols=FALSE,cluster_rows = TRUE,fontsize = 6, cutree_rows = 2, cellheight = 5, border_color=FALSE) #, annotation_col = anno)
dev.off()


rownames(ND)
