source("/Users/fazal2/Desktop/Rscripts/summarizeFit.R")
dataDirectory <- "/Users/fazal2/Documents/MethylationData"
# list the files
list.files(dataDirectory, recursive=TRUE)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(org.Hs.eg.db)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(DMRcate)
library(stringr)
library(ggplot2)
library(gplots) 
library(msigdbr)
library(gplots)
library(rgl)
library(Glimma)

# get the EPIC annotation data
EPICanno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(EPICanno)
EPICanno=as.data.frame(EPICanno)

#write.table(EPICanno2, file="/Users/fazal2/Desktop/EPICanno2.txt", sep="\t", row.names=TRUE)

targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet_working.csv")
targets

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet
head(rgSet)

targets$Label <- targets$Sample_Name
rownames(targets) <- targets$Label
sampleNames(rgSet) <- targets$Label
rgSet

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
dim(detP)
#write.table(detP, file="/Users/fazal2/Desktop/detp.txt", sep="\t", row.names=TRUE)

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(12,"Set3")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8,ylab="Mean detection p-values")
abline(h=0.01,col="red")
legend("topright", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylim = c(0,0.002), ylab="Mean detection p-values")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

coul = brewer.pal(8, "Dark2") 
# I can add more tones to this palette :
#coul = colorRampPalette(coul)(25)
densityPlot(rgSet, sampGroups=targets$Sample_Type,main="Density Plot by Cell Type", pal= coul, legend=FALSE)
legend("topleft", legend = levels(factor(targets$Sample_Type)), text.col=coul)

densityBeanPlot(rgSet, sampGroups = targets$Sample_Type, main="Bean Plot by Cell Type", pal= coul)
legend("topleft", legend = levels(factor(targets$Sample_Type)), text.col=coul)

#qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Type,
         #pdf="qcReport_All.pdf")

# remove poor quality samples
keep <- colMeans(detP) < 0.01
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:7]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessFunnorm(rgSet)
#mSetSq_N <- preprocessNoob(rgSet)
# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(getBeta(mSetRaw), sampGroups=targets$Sample_Type,main="Raw", pal= coul, legend=FALSE)
legend("topleft", legend = levels(factor(targets$Sample_Type)), text.col=coul)

densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Type,
            main="Funnorm_Normalized", legend=FALSE)
legend("topleft", legend = levels(factor(targets$Sample_Type)), text.col=coul)

#densityPlot(getBeta(mSetSq_N), sampGroups=targets$Sample_Type,
            #main="Normalized_Noob", legend=FALSE)
#legend("topleft", legend = levels(factor(targets$Sample_Type)), text.col=coul)

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
#keep <- rowSums(detP < 0.01) >= 18
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
snps <- getSnpInfo(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)#, snps = c("SBE","CpG","Probe"), maf  = 0.05)
mSetSqFlt

# exclude cross reactive probes of CpG Sites
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "13059_2016_1066_MOESM1_ESM.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt

coul = brewer.pal(12, "Set3") 
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:6])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:6])
# get the table of results for the first contrast (naive - rTreg)
EPICannoSub <- EPICanno[match(rownames(bVals),EPICanno$Name),
                        c(1:4,12:19,22:ncol(EPICanno))]
#write.table(bVals, file="/Users/fazal2/Desktop/Normalized_Filtered_BetaValue_test.txt", sep="\t", row.names=TRUE)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))

coul = brewer.pal(8, "Set2") 
par(mfrow=c(1,1))
plotMDS(getBeta(mSetSqFlt), gene.selection="common",
        col=pal[factor(targets$Sample_Name)],dim=c(1,2), cex=0.8)
legend("top", legend=levels(factor(targets$Sample_Name)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getBeta(mSetSqFlt), gene.selection="common",
        col=pal[factor(targets$Sample_Group)],dim=c(1,2))
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getBeta(mSetSqFlt), gene.selection="common",
        col=pal[factor(targets$Sample_Type)],dim=c(1,2),pch=20)
legend("right", legend=levels(factor(targets$Sample_Type)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getBeta(mSetSqFlt), gene.selection="common",
        col=pal[factor(targets$Sample_Type)],dim=c(1,2))
legend("right", legend=levels(factor(targets$Sample_Type)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getBeta(mSetSqFlt), gene.selection="common",
        col=pal[factor(targets$Source)],dim=c(1,2))
legend("right", legend=levels(factor(targets$Source)), text.col=pal,
       cex=0.7, bg="white")

# this is the factor of interest
cellType <- factor(targets$Sample_Group)
design <- model.matrix(~0+cellType)
colnames(design) <- levels(cellType)
design

#create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(C2vsNT21 = C2-NT21,
                            #A4vsNT21 = A4-NT21,
                            B3vsEP   = B3-EP,
                            C1vsEP   = C1-EP,
                            B4vsK33  = B4-K33,
                            levels=design)
contMatrix

# fit the linear model
fit <- lmFit(bVals, design)
fit2 <- eBayes(contrasts.fit(fit, contMatrix),trend=T)
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
EPICannoSub <- EPICanno[match(rownames(bVals),EPICanno$Name),
                        c(1:4,12:19,22:ncol(EPICanno))]
#out2 <- summarizeFit(fit2)
#out2 <- cbind(EPICannoSub, out2)
#head(out2, n=100)
#write.table(out2, file="/Users/fazal2/Desktop/MPs_All.txt", sep="\t", row.names=TRUE)

#A4_DMPs <- topTable(fit2,  num=Inf, coef='A4vsNT21',genelist=EPICannoSub)
#C2_DMPs <- topTable(fit2,  num=Inf, coef='C2vsNT21', genelist=EPICannoSub)
B3_DMPs <- topTable(fit2,  num=Inf, coef='B3vsEP', genelist=EPICannoSub)
#C1_DMPs <- topTable(fit2,  num=Inf, coef='C1vsEP', genelist=EPICannoSub)
#B4_DMPs <- topTable(fit2,  num=Inf, coef='B4vsK33', genelist=EPICannoSub)

#write.table(A4_DMPs, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_AllCellines/MPs_A4vsNT21.txt", sep="\t", row.names=TRUE, quote = FALSE)
#write.table(C2_DMPs, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/MPs_C2vsNT21.txt", sep="\t", row.names=TRUE, quote = FALSE)
write.table(B3_DMPs, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/MPs_B3vsEP_b4ill.txt", sep="\t", row.names=TRUE, quote = FALSE)
#write.table(C1_DMPs, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/MPs_C1vsEP.txt", sep="\t", row.names=TRUE, quote = FALSE)
#write.table(B4_DMPs, file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_Allrelaxed/MPs_B4vsK33.txt", sep="\t", row.names=TRUE, quote = FALSE)

#DMRs
myAnnotation <- cpg.annotate(object = bVals, datatype = "array", what = "Beta",
                             analysis.type= "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "C2vsNT21", arraytype = "EPIC")

str(myAnnotation)

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# convert the regions to annotated genomic ranges
#data(dmrcatedata)
results.ranges <- extractRanges(DMRs, genome = "hg19")
head(results.ranges)
write.table(results.ranges,file="/Users/fazal2/Desktop/MethylationData_Updated/Methylation_CorelationAnalysis_Selected/DMRs_C2vsNT21.txt", sep="\t", row.names=TRUE, quote = FALSE )


#groups <- pal[1:length(unique(targets$Sample_Group))]
#names(groups) <- levels(factor(targets$Sample_Group))
#cols <- groups[as.character(factor(targets$Sample_Group))]
#samps <- 1:nrow(targets)



#par(mfrow=c(1,1))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=bVals, phen.col=cols, what = "Beta",
         #arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         #genome="hg19", samps=samps)



B3=B3_DMPs[B3_DMPs$adj.P.Val<0.01 & (B3_DMPs$logFC >= 0.3 |  B3_DMPs$logFC <= -0.3),] 
C1=C1_DMPs[C1_DMPs$adj.P.Val<0.01 & (C1_DMPs$logFC >= 0.3 |  C1_DMPs$logFC <= -0.3),] 
A4=A4_DMPs[A4_DMPs$adj.P.Val<0.01 & (A4_DMPs$logFC >= 0.3 |  A4_DMPs$logFC <= -0.3),] 
C2=C2_DMPs[C2_DMPs$adj.P.Val<0.05 & (C2_DMPs$logFC >= 0.2 |  C2_DMPs$logFC <= -0.2),] 
B4=B4_DMPs[B4_DMPs$adj.P.Val<0.05 & (B4_DMPs$logFC >= 0.2 |  B4_DMPs$logFC <= -0.2),] 

B3$Relation_to_Island <- as.character(B3$Relation_to_Island)
B3$Relation_to_Island[B3$Relation_to_Island == "N_Shore" | B3$Relation_to_Island == "S_Shore" ] <- "Shore"
B3$Relation_to_Island[B3$Relation_to_Island == "N_Shelf" | B3$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(B3)

#C1$Relation_to_Island <- as.character(C1$Relation_to_Island)
C1$Relation_to_Island[C1$Relation_to_Island == "N_Shore" | C1$Relation_to_Island == "S_Shore" ] <- "Shore"
C1$Relation_to_Island[C1$Relation_to_Island == "N_Shelf" | C1$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(C1)

#A4$Relation_to_Island <- as.character(A4$Relation_to_Island)
A4$Relation_to_Island[A4$Relation_to_Island == "N_Shore" | A4$Relation_to_Island == "S_Shore" ] <- "Shore"
A4$Relation_to_Island[A4$Relation_to_Island == "N_Shelf" | A4$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(A4)

#C2$Relation_to_Island <- as.character(C2$Relation_to_Island)
C2$Relation_to_Island[C2$Relation_to_Island == "N_Shore" | C2$Relation_to_Island == "S_Shore" ] <- "Shore"
C2$Relation_to_Island[C2$Relation_to_Island == "N_Shelf" | C2$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(C2)

#B4$Relation_to_Island <- as.character(B4$Relation_to_Island)
B4$Relation_to_Island[B4$Relation_to_Island == "N_Shore" | B4$Relation_to_Island == "S_Shore" ] <- "Shore"
B4$Relation_to_Island[B4$Relation_to_Island == "N_Shelf" | B4$Relation_to_Island == "S_Shelf" ] <- "Shelf"
head(B4)

AllB3=count(B3$Relation_to_Island)
AllC1=count(C1$Relation_to_Island)
AllA4=count(A4$Relation_to_Island)
AllC2=count(C2$Relation_to_Island)
AllB4=count(B4$Relation_to_Island)



# Get the significant CpG sites at less than 5% FDR
#sigCpGs=B3_DMPs$Name[B3_DMPs$adj.P.Val<0.05 & (B3_DMPs$logFC >= 0.2 |  B3_DMPs$logFC <= -0.2)] 
#sigCpGs=B3_DMPs[B3_DMPs$adj.P.Val<0.05,] 
#head(sigCpGs)
#sigCpGs_hypo=sigCpGs[sigCpGs$logFC <= -0.2,]
#sigCpGs_hyper=sigCpGs[sigCpGs$logFC >= 0.2,] 
#write.table(sigCpGs_hypo,file="/Users/fazal2/Desktop/Hypo_B3vsEP.txt", sep="\t", row.names=TRUE, quote = FALSE )
#write.table(sigCpGs_hyper,file="/Users/fazal2/Desktop/Hyper_B3vsEP.txt", sep="\t", row.names=TRUE, quote = FALSE )

#sigCpGs_hypo=sigCpGs$Name[sigCpGs$logFC <= -0.2]
#sigCpGs_hyper=sigCpGs$Name[sigCpGs$logFC >= 0.2]
all <- B3_DMPs$Name

#B3_HypermappedEz <- getMappedEntrezIDs(sigCpGs_hyper,all,array.type="EPIC")
#B3_HypomappedEz <- getMappedEntrezIDs(sigCpGs_hypo,all,array.type="EPIC")

#All_freq <- mapIds(org.Hs.eg.db, as.character(B3_HypermappedEz$freq), 'SYMBOL', 'ENTREZID')
#write.table(B3_HypermappedEz$freq, file="/Users/fazal2/Desktop/All_freq_Chip.txt", sep="\t", row.names=TRUE)
#B3_HypermappedSymbol <- mapIds(org.Hs.eg.db, as.character(B3_HypermappedEz$sig.eg), 'SYMBOL', 'ENTREZID')
#B3_HypomappedSymbol <- mapIds(org.Hs.eg.db, as.character(B3_HypomappedEz$sig.eg), 'SYMBOL', 'ENTREZID')

#write.table(B3_HypermappedSymbol, file="/Users/fazal2/Desktop/B3_HypermappedSymbol.txt", sep="\t", row.names=TRUE)
#write.table(B3_HypomappedSymbol, file="/Users/fazal2/Desktop/B4_HypomappedSymbol.txt", sep="\t", row.names=TRUE)


# Get all the CpG sites used in the analysis to form the background
#all <- B3_DMPs$Name
#Total number of CpG sites tested
length(all)
#set <- read.table("/Users/fazal2/Desktop/UpdatedC2_GSEA/c2.all.v6.2.entrez.txt", quote="", comment="", fill = TRUE, sep="\t", header=FALSE, na.strings=c("NA","NaN", " ", "?"))
load(paste(dataDirectory,"human_c2_v5.rdata",sep="/"))

C2 = msigdbr(species = "Homo sapiens", category = "C2")
C2_Symbols = C2 %>% split(x = .$gene_symbol, f = .$gs_name)
C2_enterz = C2 %>% split(x = .$entrez_gene, f = .$gs_name)
#cat(sapply(C2_enterz, toString), file="test.txt", sep="\n", row.names(C2_enterz))
gsa_full <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=Hs.c2, plot.bias = TRUE, prior.prob = TRUE, array.type="EPIC", anno = EPICanno, fract.counts = TRUE)
gsa_hypo <- gsameth(sig.cpg=sigCpGs_hypo, all.cpg=all, collection=C2_enterz, plot.bias = TRUE, prior.prob = TRUE, array.type="EPIC", anno = EPICanno, fract.counts = TRUE)
gsa_hyper <- gsameth(sig.cpg=sigCpGs_hyper, all.cpg=all, collection=C2_enterz, plot.bias = TRUE, prior.prob = TRUE, array.type="EPIC",  anno = EPICanno, fract.counts = TRUE )
# top 10 gene sets
topGSA(gsa_full)
topGSA(gsa_hypo)
topGSA(gsa_hyper)
#full=topGSA(gsa_full,number=30)
#Hypo=topGSA(gsa_hypo)
#Hyper=topGSA(gsa_hyper)
#write.table(Hypo, file="/Users/fazal2/Desktop/MethylationData_Updated/Hypo_B4_C6_GSEA.txt", sep="\t", row.names=TRUE)
#write.table(Hyper, file="/Users/fazal2/Desktop/MethylationData_Updated/Hyper_B4_C6_GSEA.txt", sep="\t", row.names=TRUE)

library("ChAMP")
dataDirectory2 <- "./chAMP"
myLoad <- champ.load(dataDirectory2, arraytype="EPIC", filterXY=FALSE, method="minfi", filterNoCG=FALSE)
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
myebGSEA<- champ.ebGSEA(beta=bVals,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")

