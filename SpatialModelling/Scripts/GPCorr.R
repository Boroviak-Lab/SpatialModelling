



D1<-read.table("~/Desktop/CS6-update.csv",sep=",", header = T, row.names=1)
D2<-read.table("~/Desktop/CS7-update.csv",sep=",", header = T, row.names=1)
D3<-read.table("~/Desktop/CS5.csv",sep=",", header = T, row.names=1)




D2B <- cbind(D3[which(D1[,6]>0 & D2[,6]>0 & (D1[,5]>1 | D2[,5]>1) ), 7:1006], D1[which(D1[,6]>0 & D2[,6]>0 & (D1[,5]>1 | D2[,5]>1) ), 7:1006], D2[which(D1[,6]>0 & D2[,6]>0 & (D1[,5]>1 | D2[,5]>1) ), 7:1006])

#D4 <- D2[order(D2$BF,decreasing = TRUE),]

#mydata.cor = cor(t(D4[1:500,7:1006]), method = c("spearman"))

mat_breaks <- seq(-2, 5, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
out <- pheatmap(D2B,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, cutree_rows = 20, show_rownames=F, show_colnames=F,  filename = "~/Desktop/GPtest.pdf")



mat_breaks <- seq(-2, 3, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
out <- pheatmap(D2B,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, cutree_rows = 20, show_rownames=F, show_colnames=F,  filename = "~/Desktop/GPtest2.pdf")


mat_breaks <- seq(-2, 4, length.out = 20)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
out <- pheatmap(D2B,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=FALSE, scale ="row" ,cutree_rows = 20, show_rownames=F, show_colnames=F,  filename = "~/Desktop/GPtest3.pdf")


#C1 <- sort(cutree(out$tree_row, k=2))
#C2 <- sort(cutree(out$tree_row, k=3))
#C3 <- sort(cutree(out$tree_row, k=4))
#C4 <- sort(cutree(out$tree_row, k=5))

#write.csv(as.data.frame(C1), file=paste("~/Desktop/CS6_Cl_k=2.csv",sep=""))
#write.csv(as.data.frame(C2), file=paste("~/Desktop/CS6_Cl_k=3.csv",sep=""))
#write.csv(as.data.frame(C3), file=paste("~/Desktop/CS6_Cl_k=4.csv",sep=""))
#write.csv(as.data.frame(C4), file=paste("~/Desktop/CS6_Cl_k=5.csv",sep=""))


#genes1 <- rownames(as.data.frame(C4[which(C4==1)]))
#genes2 <-rownames(as.data.frame(C4[which(C4==2)]))
#genes3 <-rownames(as.data.frame(C4[which(C4==3)]))
#genes4 <-rownames(as.data.frame(C4[which(C4==4)]))
#genes5 <-rownames(as.data.frame(C4[which(C4==5)]))#

#my_genes1 <- D4[row.names(D4) %in% genes1 ,7:1006]
#out <- pheatmap(my_genes1,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/CS6_k=5_Cl1.pdf")

#my_genes2 <- D4[row.names(D4) %in% genes2 ,7:1006]
#out <- pheatmap(my_genes2,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/CS6_k=5_Cl2.pdf")

#my_genes3 <- D4[row.names(D4) %in% genes3 ,7:1006]
#out <- pheatmap(my_genes3,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/CS6_k=5_Cl3.pdf")

#my_genes4 <- D4[row.names(D4) %in% genes4 ,7:1006]
#out <- pheatmap(my_genes4,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/CS6_k=5_Cl4.pdf")

#my_genes5 <- D4[row.names(D4) %in% genes5 ,7:1006]
#out <- pheatmap(my_genes5,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/CS6_k=5_Cl5.pdf")

