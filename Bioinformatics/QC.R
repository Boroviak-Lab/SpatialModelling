library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("tidyverse")
library("DropletUtils")
library(gtools)
set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "QC"
dir.create(saveext)

#Read in the key
BS<-read.table("../Data/marmKey.csv",sep=",",header = T, row.names=1)

CS <- as.character(BS$Carnegie.Stage)
CS[which(CS=="CS1" | CS=="CS2" | CS=="CS3")] <- "CS1-3"
BS$Carnegie.Stage <- as.factor(CS)

stage <- which( (BS$Carnegie.Stage=="CS1-3" | BS$Carnegie.Stage=="CS5" | BS$Carnegie.Stage=="CS6" | BS$Carnegie.Stage=="CS7") & BS$QC==1 )
BS2 <- BS[stage,]

p <- ggplot(data=BS2, aes(x=Carnegie.Stage, y=Unique_Reads)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads", breaks=c(0,500000,1000000,5000000,10000000))
ggsave(filename=paste("QC/ReadsByStage.pdf",sep=""), plot = p, width = 16, height = 16)

p <- ggplot(data=BS2, aes(x=Carnegie.Stage, y=Unique_Reads)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads", breaks=c(0,500000,1000000,5000000,10000000))
ggsave(filename=paste("QC/ReadsByStage_box.pdf",sep=""), plot = p, width = 16, height = 16)

stage <- which( (BS$Carnegie.Stage=="CS1-3" |BS$Carnegie.Stage=="CS5" | BS$Carnegie.Stage=="CS6" | BS$Carnegie.Stage=="CS7") & BS$QC==1 )
BS2 <- BS[stage,]

p <- ggplot(data=BS2, aes(x=Carnegie.Stage, y=NoGenes)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste("QC/NoGenesByStage.pdf",sep=""), plot = p, width = 16, height = 16)

p <- ggplot(data=BS2, aes(x=Carnegie.Stage, y=NoGenes)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()
ggsave(filename=paste("QC/NoGenesByStage_boxplot.pdf",sep=""), plot = p, width = 16, height = 16)


stage2 <- which( (BS$Annotation2=="EmDisc_CS5" | BS$Annotation2=="EmDisc_CS6" | BS$Annotation2=="EmDisc_CS7" | BS$Annotation2=="ExMes_CS5" | BS$Annotation2=="ExMes_CS6" | BS$Annotation2=="ExMes_CS7" | BS$Annotation2=="Gland_CS5" | BS$Annotation2=="Gland_CS6" | BS$Annotation2=="Gland_CS7" | BS$Annotation2=="Tb_CS5_" | BS$Annotation2=="Tb_CS6" | BS$Annotation2=="Tb_CS7" | BS$Annotation2=="SYS_CS5"  | BS$Annotation2=="SYS_CS6" | BS$Annotation2=="SYS_CS7") & BS$QC==1 )
BS3 <- BS[stage2,]

ann <- as.character(BS3$Annotation2)

ann[which(ann=="EmDisc_CS5")] <- "1 EmDisc"
ann[which(ann=="EmDisc_CS6")] <- "1 EmDisc"
ann[which(ann=="EmDisc_CS7")] <- "1 EmDisc"
ann[which(ann=="Tb_CS5")] <- "4 Tb"
ann[which(ann=="Tb_CS6")] <- "4 Tb"
ann[which(ann=="Tb_CS7")] <- "4 Tb"
ann[which(ann=="Gland_CS5")] <- "5 Gland"
ann[which(ann=="Gland_CS6")] <- "5 Gland"
ann[which(ann=="Gland_CS7")] <- "5 Gland"
ann[which(ann=="SYS_CS5")] <- "3 SYS"
ann[which(ann=="SYS_CS6")] <- "3 SYS"
ann[which(ann=="SYS_CS7")] <- "3 SYS"
ann[which(ann=="ExMes_CS5")] <- "2 ExMes"
ann[which(ann=="ExMes_CS6")] <- "2 ExMes"
ann[which(ann=="ExMes_CS7")] <- "2 ExMes"

BS3$ann <- as.factor(ann) 

p <- ggplot(data=BS3, aes(x=ann, y=Unique_Reads)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads", breaks=c(0,500000,1000000,5000000,10000000))
ggsave(filename=paste("QC/GenesMappedByStage.pdf",sep=""), plot = p, width = 16, height = 16)

p <- ggplot(data=BS3, aes(x=ann, y=Unique_Reads)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads", breaks=c(0,500000,1000000,5000000,10000000))
ggsave(filename=paste("QC/GenesMappedByStage_boxplot.pdf",sep=""), plot = p, width = 16, height = 16)

p <- ggplot(data=BS3, aes(x=ann, y=NoGenes)) + geom_violin(trim=FALSE)  + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads")
ggsave(filename=paste("QC/NoGenesByType.pdf",sep=""), plot = p, width = 16, height = 16)

p <- ggplot(data=BS3, aes(x=ann, y=NoGenes)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="Unique reads")
ggsave(filename=paste("QC/NoGenesByType_boxplot.pdf",sep=""), plot = p, width = 16, height = 16)


#Load the marmoset data ... to get ordering
raw_counts3<-read.table("../Data/marmData.csv",sep=",",header = T, row.names=1)
isinvit <- BS$PostImp * BS$QC
raw_counts3<- raw_counts3[,which(isinvit>0)]

D2 <- sign(as.data.frame(raw_counts3) )
write.csv(colSums(D2), file=paste("QC/QC_NumberofGenes.csv",sep=""))

D <- log10( as.data.frame(raw_counts3) +1 )

x2 <- raw_counts3
Counts <- colSums(sign(x2))
C <- as.data.frame(Counts)
colnames(C) <- "No_genes"
C$Frac <- 1

downs <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01)

for (j in 1:length(downs)){

  x2 <- downsampleMatrix(as.matrix(raw_counts3),downs[j], bycol=TRUE)
  Counts <- colSums(sign(x2))
  C2 <- as.data.frame(Counts)
  colnames(C2) <- "No_genes"
  C2$Frac <- downs[j]
  
  C <- rbind(C,C2)
  
}

C$No_genes <- log2(C$No_genes)
C$Fraction <- as.factor(C$Frac)

p <- ggplot(data=C, aes(x=Fraction, y=No_genes)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() + scale_y_continuous(name="No genes",breaks = log2(c(1000,2000,3000,4000,5000,10000,20000)))
ggsave(filename=paste("QC/NoGenesDownSamplelog2_2.pdf",sep=""), plot = p, width = 24, height = 8)


C$Fraction <- as.factor(1-C$Frac)

p <- ggplot(data=C, aes(x=Fraction, y=No_genes)) + geom_violin(trim=FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() + scale_y_continuous(name="No genes",breaks = log2(c(1000,2000,3000,4000,5000,10000,20000)))
ggsave(filename=paste("QC/NoGenesDownSamplelog2_1.pdf",sep=""), plot = p, width = 24, height = 8)


C$Fraction <- as.factor(1-C$Frac)
p <- ggplot(data=C, aes(x=Fraction, y=No_genes)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic() +scale_y_continuous(name="No genes")
ggsave(filename=paste("QC/NoGenesDownSample1_boxplot.pdf",sep=""), plot = p, width = 16, height = 16)


ids = which(BS$ForMapping==1)

Data <- data.frame(cbind( BS$X.Mapped_to_annotation_base[ids],BS$X.Mapped_to_annotation_1kbextended[ids],BS$X.Mapped_to_annotation_5kbextended[ids], apply(cbind(BS$X.Mapped_to_annotation_base[ids],BS$X.Mapped_to_annotation_5kbextended[ids],BS$X.Mapped_to_annotation_1kbextended[ids]), 1, FUN = sum) ) ) #, col.names = c('X', 'y')) )
colnames(Data) <- c("Base","Extended (1kb)","Extended (5kb)","Sum")

x <- seq(from = 1, to = dim(Data)[1], by = 1)

Data1 = Data[order(Data$Base),]
Data1$t = x

D1 <- Data1[c("t","Base")]
D1$Type <- "Base"
colnames(D1) <- c("t","Efficiency","Type")

D2 <- Data1[c("t","Extended (1kb)")]
D2$Type <- "Extended (1kb)"
colnames(D2) <- c("t","Efficiency","Type")

D3 <- Data1[c("t","Extended (5kb)")]
D3$Type <- "Extended (5kb)"
colnames(D3) <- c("t","Efficiency","Type")

library("ggpubr")

Data2 <- rbind(D1,D2,D3)

p <- ggplot(data=Data2, aes(x=t, y=Efficiency,group=Type)) +
     geom_line(aes(color=Type))+
     scale_color_grey() + theme_classic()
ggsave(filename=paste("QC/Mapped_to_gtf",".pdf",sep=""), plot = p, width = 16, height = 16)

Data2 <- rbind(D1,D3)
ggpaired(Data2, x = "Type", y = "Efficiency",
         color = "Type", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)
ggsave(filename=paste("QC/Mapped_to_gtf_Basevs5kb",".pdf",sep=""), width = 16, height = 16)



#compare_means(Efficiency ~ Type,  data = Data2)

compare_means(Efficiency ~ Type,  data = Data2, paired = TRUE)

my_comparisons <- list( c("Base", "Extended (1kb)"), c("Base", "Extended (5kb)"))
p <- ggpaired(Data2, x = "Type", y = "Efficiency", color = "Type", palette = "jco") +geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4) + theme_classic()
ggsave(filename=paste("QC/Mapped_to_gtf_Basevs5kb",".pdf",sep=""), width = 16, height = 16)


ggpaired(Data2, x = "Type", y = "Efficiency",
         color = "Type", line.color = "gray", line.size = 0.4,
         palette = "jco")


#ggpaired(Data2, x = "t", y = "Efficiency",
#         group = "Type", color = "supp", line.color = "gray", line.size = 0.4,palette = "jco")+stat_compare_means(paired = TRUE)


#Now generate the expression densit plots
library(reshape2)

#Now plot QC as a funnction of time ...
Locs <- BS[is.na(BS$Location)==FALSE,]
locs <- Locs$Location
eff <- Locs$X.Mapped


#For the E25A dataset
E25A_ID <- locs[grep("E25A", locs) ]
E25A_Eff <- eff[grep("E25A", locs) ]
E25A_ID <- gsub("_", "-", E25A_ID)
E25A_ID <- as.data.frame(E25A_ID)
colnames(E25A_ID) <- "test"
E25A_ID <- separate(E25A_ID, test, into = c("Stage","Slide","Time"), "-")
E25A_ID$Eff <- E25A_Eff

t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")
for (i in 1:length(t1)){
  E25A_ID$Time <- replace(E25A_ID$Time, E25A_ID$Time == t1[i], t2[i])
}

E25A_Rc = matrix(, nrow = length(unique(E25A_ID$Slide)), ncol = 2)
count <- 0
for (i in unique(E25A_ID$Slide)) {
x <- E25A_ID[which(E25A_ID$Slide==as.double(i)),]

t <- sort(unique(x$Time))
t2 <- seq(1,length(t),by=1)
t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
count <- count +1
E25A_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
E25A_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")

x$T <- t3
p <- ggplot(data=x, aes(x=T, y=Eff)) +
  geom_point() + theme_classic()+ylim(0, 1) + 
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle(paste("R = ", E25A_Rc[count,1],sep=""))

ggsave(filename=paste("QC/E25A_",i,".pdf",sep=""), plot = p, width = 16, height = 16)

}



#E25B
E25B_ID <- locs[grep("E25B", locs) ]
E25B_Eff <- eff[grep("E25B", locs) ]
E25B_ID <- gsub("_", "-", E25B_ID)
E25B_ID <- as.data.frame(E25B_ID)
colnames(E25B_ID) <- "test"
E25B_ID <- separate(E25B_ID, test, into = c("Stage","Slide","Time"), "-")
E25B_ID$Eff <- E25B_Eff

t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")
for (i in 1:length(t1)){
  E25B_ID$Time <- replace(E25B_ID$Time, E25B_ID$Time == t1[i], t2[i])
}

E25B_Rc = matrix(, nrow = length(unique(E25B_ID$Slide)), ncol = 2)
count <- 0
for (i in unique(E25B_ID$Slide)) {
  x <- E25B_ID[which(E25B_ID$Slide==as.double(i)),]

  t <- sort(unique(x$Time))
  t2 <- seq(1,length(t),by=1)
  t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
  count <- count +1
  E25B_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
  E25B_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")
  
  x$T <- t3
  p <- ggplot(data=x, aes(x=T, y=Eff)) +
    geom_point() + theme_classic()+ylim(0, 1) + 
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle(paste("R = ", E25B_Rc[count,1],sep=""))
  
  ggsave(filename=paste("QC/E25B_",i,".pdf",sep=""), plot = p, width = 16, height = 16)
  
  
}

#E25D
E25D_ID <- locs[grep("E25D", locs) ]
E25D_Eff <- eff[grep("E25D", locs) ]
E25D_ID <- gsub("_", "-", E25D_ID)
E25D_ID <- as.data.frame(E25D_ID)
colnames(E25D_ID) <- "test"
E25D_ID <- separate(E25D_ID, test, into = c("Stage","Slide","Time"), "-")
E25D_ID$Eff <- E25D_Eff
t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")

for (i in 1:length(t1)){
  E25D_ID$Time <- replace(E25D_ID$Time, E25D_ID$Time == t1[i], t2[i])
}

E25D_Rc = matrix(, nrow = length(unique(E25D_ID$Slide)), ncol = 2)
count <- 0
for (i in unique(E25D_ID$Slide)) {
  x <- E25D_ID[which(E25D_ID$Slide==as.double(i)),]

  t <- sort(unique(x$Time))
  t2 <- seq(1,length(t),by=1)
  t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
  count <- count +1
  E25D_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
  E25D_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")

  x$T <- t3
  p <- ggplot(data=x, aes(x=T, y=Eff)) +
    geom_point() + theme_classic()+ylim(0, 1) + 
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle(paste("R = ", E25D_Rc[count,1],sep=""))
  
  ggsave(filename=paste("QC/E25D_",i,".pdf",sep=""), plot = p, width = 16, height = 16)
  
  
}


#E15C
E15C_ID <- locs[grep("E15C", locs) ]
E15C_Eff <- eff[grep("E15C", locs) ]
E15C_ID <- gsub("_", "-", E15C_ID)
E15C_ID <- as.data.frame(E15C_ID)
colnames(E15C_ID) <- "test"
E15C_ID <- separate(E15C_ID, test, into = c("Stage","Slide","Time"), "-")
E15C_ID$Eff <- E15C_Eff
t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")

for (i in 1:length(t1)){
  E15C_ID$Time <- replace(E15C_ID$Time, E15C_ID$Time == t1[i], t2[i])
}

E15C_Rc = matrix(, nrow = length(unique(E15C_ID$Slide)), ncol = 2)
count <- 0
for (i in unique(E15C_ID$Slide)) {
  x <- E15C_ID[which(E15C_ID$Slide==as.double(i)),]

  t <- sort(unique(x$Time))
  t2 <- seq(1,length(t),by=1)
  t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
  count <- count +1
  E15C_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
  E15C_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")
  
  x$T <- t3
  p <- ggplot(data=x, aes(x=T, y=Eff)) +
    geom_point() + theme_classic()+ylim(0, 1) + 
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle(paste("R = ", E15C_Rc[count,1],sep=""))
  
  ggsave(filename=paste("QC/E15C_",i,".pdf",sep=""), plot = p, width = 16, height = 16)
  
    
}


#E15B1
E15B1_ID <- locs[grep("P1_E15B", locs) ]
E15B1_Eff <- eff[grep("P1_E15B", locs) ]
E15B1_ID <- gsub("_", "-", E15B1_ID)
E15B1_ID <- as.data.frame(E15B1_ID)
colnames(E15B1_ID) <- "test"
E15B1_ID <- separate(E15B1_ID, test, into = c("Plate","Stage","Slide","Time"), "-")
E15B1_ID$Eff <- E15B1_Eff
t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")


for (i in 1:length(t1)){
  E15B1_ID$Time <- replace(E15B1_ID$Time, E15B1_ID$Time == t1[i], t2[i])
}


E15B1_Rc = matrix(, nrow = length(unique(E15B1_ID$Slide)), ncol = 2)
count <- 0
for (i in unique(E15B1_ID$Slide)) {
  x <- E15B1_ID[which(E15B1_ID$Slide==as.double(i)),]

  t <- sort(unique(x$Time))
  t2 <- seq(1,length(t),by=1)
  t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
  count <- count +1
  E15B1_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
  E15B1_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")
  
  x$T <- t3
  p <- ggplot(data=x, aes(x=T, y=Eff)) +
    geom_point() + theme_classic()+ylim(0, 1) + 
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle(paste("R = ", E15B1_Rc[count,1],sep=""))
  
  ggsave(filename=paste("QC/E15B1_",i,".pdf",sep=""), plot = p, width = 16, height = 16)
  
  
}


#E15B2
E15B2_ID <- locs[grep("P2_E15B", locs) ]
E15B2_Eff <- eff[grep("P2_E15B", locs) ]
E15B2_ID <- gsub("_", "-", E15B2_ID)
E15B2_ID <- as.data.frame(E15B2_ID)
colnames(E15B2_ID) <- "test"
E15B2_ID <- separate(E15B2_ID, test, into = c("Plate","Stage","Slide","Time"), "-")
E15B2_ID$Eff <- E15B2_Eff
t1 <- c("AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD")
t2 <- c("ZA","ZB","ZC","ZD","ZE","ZF","ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV","ZW","ZX","ZY","ZZ","ZZA","ZZB","ZZC","ZZD")

for (i in 1:length(t1)){
  E15B2_ID$Time <- replace(E15B2_ID$Time, E15B2_ID$Time == t1[i], t2[i])
}

E15B2_Rc = matrix(, nrow = length(unique(E15B2_ID$Slide)), ncol = 2)
count <- 0

for (i in unique(E15B2_ID$Slide)) {
  x <- E15B2_ID[which(E15B2_ID$Slide==as.double(i)),]
 
  t <- sort(unique(x$Time))
  t2 <- seq(1,length(t),by=1)
  t3 <- seq(1,length(t),by=1)
  for (j in 1:length(t3)) {
    inds <- which(x$Time==t[j])
    t3[j] <- t2[inds]
  }
  count <- count +1
  E15B2_Rc[count,1] <- cor(t3,x$Eff,method = "spearman")
  E15B2_Rc[count,2] <- cor(permute(t3),x$Eff,method = "spearman")
  
  
  
  x$T <- t3
  p <- ggplot(data=x, aes(x=T, y=Eff)) +
    geom_point() + theme_classic()+ylim(0, 1) + 
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle(paste("R = ", E15B2_Rc[count,1],sep=""))
  
  ggsave(filename=paste("QC/E15B2_",i,".pdf",sep=""), plot = p, width = 16, height = 16)
  
}

AllCorrelations <- as.data.frame(rbind(E25A_Rc,E25B_Rc,E25D_Rc,E15C_Rc,E15B1_Rc,E15B2_Rc))
colnames(AllCorrelations) <- c("Pearson","Pearson (permuted)")


AllCorrelations <- na.omit(AllCorrelations)
AllCorrelations2 <- melt(AllCorrelations)
AllCorrelations2 <- AllCorrelations2[which(AllCorrelations2$variable=="Pearson"),]
p <- ggplot(AllCorrelations2, aes(x=variable, y=value)) + geom_boxplot(width=0.1) + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()

ggsave(filename=paste("QC/E15E25Correlations.pdf",sep=""), plot = p, width = 16, height = 16)