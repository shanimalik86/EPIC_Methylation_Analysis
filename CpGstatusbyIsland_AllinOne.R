library(limma)
library(minfi)
#library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(WGCNA)
library(topGO)
library(ggplot2)
library(qqman)
library(CMplot)  
library(plyr)
library(devtools)

B3 <- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/DMPs_CpGsOnly/DMPs_B3vsEP_Beta0.2_0.05FDR.txt", quote="", comment="", sep="\t", header=TRUE, na.strings=c("NA","NaN", " ", "?"))
C1 <- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/DMPs_CpGsOnly/DMPs_C1vsEP_Beta0.2_0.05FDR.txt", quote="", comment="", sep="\t", header=TRUE, na.strings=c("NA","NaN", " ", "?"))
#A4 <- read.table("/Users/fazal2/Desktop/MethylationData_Updated/DMPs_A4vsNT21_Beta0.3_0.01FDR.txt", quote="", comment="", sep="\t", header=TRUE, na.strings=c("NA","NaN", " ", "?"))
C2 <- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/DMPs_CpGsOnly/DMPs_C2vsNT21.txt", quote="", comment="", sep="\t", header=TRUE, na.strings=c("NA","NaN", " ", "?"))
B4 <- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/DMPs_CpGsOnly/DMPs_B4vsK33.txt", quote="", comment="", sep="\t", header=TRUE, na.strings=c("NA","NaN", " ", "?"))

B3$Relation_to_Island <- as.character(B3$Relation_to_Island)
B3$Relation_to_Island[B3$Relation_to_Island == "N_Shore" | B3$Relation_to_Island == "S_Shore" ] <- "Shore"
B3$Relation_to_Island[B3$Relation_to_Island == "N_Shelf" | B3$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(B3)

C1$Relation_to_Island <- as.character(C1$Relation_to_Island)
C1$Relation_to_Island[C1$Relation_to_Island == "N_Shore" | C1$Relation_to_Island == "S_Shore" ] <- "Shore"
C1$Relation_to_Island[C1$Relation_to_Island == "N_Shelf" | C1$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(C1)

#A4$Relation_to_Island <- as.character(A4$Relation_to_Island)
#A4$Relation_to_Island[A4$Relation_to_Island == "N_Shore" | A4$Relation_to_Island == "S_Shore" ] <- "Shore"
#A4$Relation_to_Island[A4$Relation_to_Island == "N_Shelf" | A4$Relation_to_Island == "S_Shelf" ] <- "Shelf"
#head(A4)

C2$Relation_to_Island <- as.character(C2$Relation_to_Island)
C2$Relation_to_Island[C2$Relation_to_Island == "N_Shore" | C2$Relation_to_Island == "S_Shore" ] <- "Shore"
C2$Relation_to_Island[C2$Relation_to_Island == "N_Shelf" | C2$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(C2)

B4$Relation_to_Island <- as.character(B4$Relation_to_Island)
B4$Relation_to_Island[B4$Relation_to_Island == "N_Shore" | B4$Relation_to_Island == "S_Shore" ] <- "Shore"
B4$Relation_to_Island[B4$Relation_to_Island == "N_Shelf" | B4$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(B4)

AllB3=count(B3$Relation_to_Island)
AllC1=count(C1$Relation_to_Island)
#AllA4=count(A4$Relation_to_Island)
AllC2=count(C2$Relation_to_Island)
AllB4=count(B4$Relation_to_Island)

# Get the significant CpG sites at less than 5% FDR
B4_sigCpGs = B4[B4$adj.P.Val<0.05, ] 
head(B4_sigCpGs)
B4_Hyper = B4_sigCpGs[B4_sigCpGs$DeltaBeta >= 0.2,] 
head(B4_Hyper)
B4_Hypo = B4_sigCpGs[B4_sigCpGs$DeltaBeta <= -0.2,]

B4_hyper_count=count(B4_Hyper$Relation_to_Island)
B4_hypo_count=count(B4_Hypo$Relation_to_Island)
#B4_hypo_count=rbind(B4_hypo_count,list("ExonBnd",0))
#B4_hypo_count <- B4_hypo_count[order(B4_hypo_count$x),]
features <- AllB4$x

df <- data.frame(AllB4$freq,B4_hyper_count$freq,B4_hypo_count$freq)

colnames(df)=c("All Significant", "Hyper","Hypo")
rownames(df)=AllB4$x
key=as.matrix(rownames(df))
coul = brewer.pal(9, "Paired")
#coul = brewer.pal(7,"BrBG")
barplot(as.matrix(df), col=coul, space=0.6, font.axis=2, ylab="portion", xpd = FALSE, main="CpGs Status -By Island", cex.axis=1.0, cex.names=0.8)
