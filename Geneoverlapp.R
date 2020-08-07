library(GeneOverlap)
library(tidyverse)
library(reshape)
library("ggpubr")
library(dplyr)
library(msigdbr)
case1 <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/C2_Case1.txt")
case1=as.character(case1$V1)
case2 <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/C2_Case2.txt")
case2=as.character(case2$V1)
case3 <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/C2_Case3.txt")
case3=as.character(case3$V1)
case4 <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/C2_Case4.txt")
case4=as.character(case4$V1)

NonCpGs <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/C1_NonCpGs.txt")
NonCpGs=as.character(NonCpGs$V1)
C2 = msigdbr(species = "Homo sapiens", category = "C2")
C2_Symbols = C2 %>% split(x = .$gene_symbol, f = .$gs_name)

background <- read.table("/Users/fazal2/Desktop/Cases_Enrichment/AllGenes_EPIC.txt")
back=as.character(background$V1)

for(i in 1:length(C2_Symbols)){
go.obj <- newGeneOverlap(C2_Symbols[[i]],
                         NonCpGs,
                         genome.size=length(back))
enrich<- testGeneOverlap(go.obj)
write.table(cbind(names(C2_Symbols[i]), getPval(enrich),getOddsRatio(enrich)), file="/Users/fazal2/Desktop/Cases_Enrichment/C1_NonCpGs_C2_Enrichment.txt", sep="\t", row.names=FALSE,append=TRUE,quote = FALSE,col.names = FALSE)
}
