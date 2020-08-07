library(EnhancedVolcano)
library(magrittr)
library(RColorBrewer)
res <- read.table("/Users/fazal2/Desktop/MethylationData_updated/Methylation_CorelationAnalysis_Selected/REMP/MPs_B3vsEP_REMPProfile_L1Volcano.txt", header=TRUE, na.strings=c("NA","NaN", " ", "?"))
keyvals <- ifelse(
  res$logFC <= -0.2 & res$padj<0.05 , 'darkgoldenrod2',
  ifelse(res$logFC >= 0.2 & res$padj<0.05, 'royalblue',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'royalblue'] <- 'Hypermythylated'
names(keyvals)[keyvals == 'black'] <- 'No Change'
names(keyvals)[keyvals == 'darkgoldenrod2'] <- 'Hypomethylated'

EnhancedVolcano(res,
                lab = "",
                x = "logFC",
                y = "padj",
                xlab = bquote("Delta Beta"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                FCcutoff = 0.2,
                pCutoff = 0.05,
                xlim = c(-1.0, 1.0),
                ylim=c(0,20),
                pointSize = 0.08,
                #pointSize = c(ifelse(res$logFC<0.4, 0.7, 1)),
                colCustom = keyvals,
                legendPosition = "right",
                legendLabSize = 12,
                legendIconSize = 3.0,
                colAlpha = 0.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = "B3 vs EP L1")
ggsave("B3vsEP_L1_SuppFigure3.eps", device=cairo_ps,width = 8, height = 6, dpi = 300)
dev.off()

down=res[res$logFC <= -0.2 & res$padj<0.05,]
up=res[res$logFC >= 0.2 & res$padj<0.05,]
dim(down)
dim(up)
