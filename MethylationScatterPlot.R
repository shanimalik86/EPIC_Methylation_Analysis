library(pheatmap)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(corrplot)
library(reshape)
library("ggpubr")
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
source("http://www.sthda.com/upload/rquery_cormat.r")
#exp<- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/NT21_C2_Geneexpression_Common_Reps_CpGsID.txt",header= TRUE)
#meth<- read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/NT21_C2_Methylation_Common_reps.txt",header= TRUE)
#exp=melt(exp)
#meth=melt(meth)
#write.table(exp, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/NT21_C2_Geneexpression_Common_Reps_scatterformated.txt", sep="\t", row.names=TRUE)
#write.table(meth, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/NT21_C2_Methylation_Common_Reps_scatterformated.txt", sep="\t", row.names=TRUE)

res=read.table("/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/K33_B4_Correlation_Matrix_Figure5.txt",header= TRUE)
#genes <- as.factor(data$Gene)
#cpgs <- as.factor(data$Probe)


#Cor = res %>% 
  #group_by(Gene, Probe) %>% 
  #mutate(corr = cor(Geneexpression, Beta))%>%
  #arrange(Gene, Probe)
#write.table(Cor, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/K33_B4_Correlation_Matrix_MultipleGenes.txt", sep="\t", row.names=TRUE)

cor2=res %>%
  group_by(Gene, Probe) %>%
  mutate(Cor = cor.test(Geneexpression, Beta)$estimate) %>%
  mutate(p.value = cor.test(Geneexpression, Beta)$p.value) %>%
  arrange(Gene, Probe)
#write.table(cor2, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/K33_B4_Correlation_Matrix_Pvalue_MultipleGenes.txt", sep="\t", row.names=TRUE)

fit=group_split(cor2, .keep = TRUE)
length(fit)

scaleFUN2 <- function(x) sprintf("%.2f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)
ppi = 200

exp_vs_meth_plot <- function(fit){
  scaleFUN2 <- function(x) sprintf("%.2f", x)
  scaleFUN1 <- function(x) sprintf("%.1f", x)
  ppi = 200
  #cor.coef=TRUE, cor.coeff.args = list(method = "pearson", label.x.npc = "middle", label.y.npc = "top")
  ggscatter(fit[[i]], x = "Geneexpression" ,  y="Beta", size= 3.5, color="Type") +
    stat_smooth(method= "lm") + 
    stat_cor(method="pearson", label.x.npc= "middle", size = 3)+
    labs(title = fit[[i]]$Gene) +
    xlab("Geneexpression") + 
    ylab("Beta") +
    #scale_x_continuous(trans = "log")
    #facet_wrap(Gene~Probe, scales= "free")+
    theme(strip.text = element_text(color= "black"), 
          panel.spacing= unit(1.5, "lines"),
          panel.grid.major = element_blank(), 
          plot.title = element_text(color = "red", size = 12, 
                                    face = "bold", hjust = 0.5),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8, linetype = "solid"), 
          panel.background = element_blank(),
          legend.position="right") +
    scale_x_continuous(labels= scaleFUN2) + 
    scale_y_continuous(labels= scaleFUN1)
}

for(i in 1:length(fit)){
p1 <- exp_vs_meth_plot(fit)
  if(as.numeric(fit[[i]]$Cor)>0 & as.numeric(fit[[i]]$FC>0))
  {
    png(paste0("/Users/fazal2/Desktop/Verify/", fit[[i]]$Gene, "vs", fit[[i]]$Probe,"_",i, ".png"), width=6*ppi, height=4*ppi, res=ppi)
    pdf("/Users/fazal2/Desktop/Verify/B4_Figure5_4.pdf") 
    print(p1)
    dev.off()
    #write.table(fit[[i]], file="/Users/fazal2/Desktop/Verify/Case4.txt", sep="\t", row.names=TRUE,append=TRUE)
  }


  if(as.numeric(fit[[i]]$Cor)>0 & as.numeric(fit[[i]]$FC<0))
  {
   png(paste0("/Users/fazal2/Desktop/Verify/", fit[[i]]$Gene, "vs", fit[[i]]$Probe,"_",i, ".png"), width=6*ppi, height=4*ppi, res=ppi)
    pdf("/Users/fazal2/Desktop/Verify/B4_Figure5_3.pdf") 
   print(p1)
   dev.off()
   #write.table(fit[[i]], file="/Users/fazal2/Desktop/Verify/Case3.txt", sep="\t", row.names=TRUE,append=TRUE)
 }

 if(as.numeric(fit[[i]]$Cor)<0 & as.numeric(fit[[i]]$FC>0))
 {
  png(paste0("/Users/fazal2/Desktop/Verify/", fit[[i]]$Gene, "vs", fit[[i]]$Probe,"_",i, ".png"), width=6*ppi, height=4*ppi, res=ppi)
   pdf("/Users/fazal2/Desktop/Verify/B4_Figure5_2.pdf") 
  print(p1)
  dev.off()
  #write.table(fit[[i]], file="/Users/fazal2/Desktop/Verify/Case2.txt", sep="\t", row.names=TRUE,append=TRUE)
 }

 if(as.numeric(fit[[i]]$Cor)<0 & as.numeric(fit[[i]]$FC<0))
 {
  png(paste0("/Users/fazal2/Desktop/Verify/", fit[[i]]$Gene, "vs", fit[[i]]$Probe,"_",i, ".png"), width=6*ppi, height=4*ppi, res=ppi)
  pdf("/Users/fazal2/Desktop/Verify/B4_Figure5_1.pdf") 
  print(p1)
  dev.off()
  #write.table(fit[[i]], file="/Users/fazal2/Desktop/Verify/Case1.txt", sep="\t", row.names=TRUE,append=TRUE)
  
 }

}

EPICanno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
head(EPICanno)
EPICanno=as.data.frame(EPICanno)
#write.table(EPICanno, file="/Users/fazal2/Desktop/MethylationData_Updated/Scatter_Data/.txt", sep="\t", row.names=TRUE)
CpG_Gene_Interest = EPICanno %>%
  dplyr::select(Name, UCSC_RefGene_Name) %>%
  tidyr::separate_rows(UCSC_RefGene_Name, sep = ";")%>%
  dplyr::distinct()

#write.table(CpG_Gene_Interest, file="/Users/fazal2/Desktop/MethylationData_Updated/AllCpGs_withMissings.txt", sep="\t", row.names=TRUE)