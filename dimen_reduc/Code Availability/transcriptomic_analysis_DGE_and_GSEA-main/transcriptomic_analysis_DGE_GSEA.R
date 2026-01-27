library(Biobase)
library(GEOquery)
library(limma)

library(RColorBrewer)
library(EnhancedVolcano)
library(RNAseq123)
library(VennDiagram)
library(gplots)
library(EGSEA)

#######################################################
# load series and platform data from GEO
#######################################################
gset <- getGEO("GSE121239", GSEMatrix =TRUE)
# if (length(gset) > 1) idx <- grep("GPL11532", attr(gset, "names")) else idx <- 1
if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group
gsms <-  paste0(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 3, 2,
                2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 3, 3, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                1, 1, 1, 1, 2, 3, 3, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 2, 3,
                3, 2, 2, 3, 3, 3, 3, 2, 2, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 2, 1, 1,
                1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 2, 3, 2, 2,
                2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3,
                3, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2,
                2, 2, 2, 2, 2, 1, 1, 1, 3, 2, 3, 2, 1, 1, 1, 2, 2, 3, 3, 2, 2, 2,
                2, 3, 3, 3, 3, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 2, 2, 1, 1, 1, 2, 3,
                3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 3, 3,
                2, 2, 1, 1, 1, 2, 2, 1, 2, 2, 2, 3, 2, 3, 2, 2, 1, 1, 2, 2, 2, 2,
                2, 1, 1, 1, 1, 1, 2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 2, 3, 2, 3, 3,
                3, 2, 2, 3, 1, 1, 1, 2, 2, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1,
                1, 1, 2, 2, 2, 3, 3, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1,
                2, 2, 1, 1)

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

#######################################################
# log2 transform
#######################################################
ex <- exprs(gset)
# rownames(ex) <- id_ref
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#######################################################
# DGE by limma
#######################################################
# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)

# design matrix
gset$description <- fl
design <- model.matrix(~0+description, gset )
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)

# set contrast 
cont.matrix <- makeContrasts(G1-G2, G3-G2, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tfit <- treat(fit2)

# extract top differentially expressed genes
tT_1_vs_2 <- topTreat(tfit, coef = 1, n = Inf)
tT_3_vs_2 <- topTreat(tfit, coef = 2, n = Inf)

#######################################################
# heatmap
#######################################################
ex <- as.matrix(ex)
tT_1_vs_2.topgenes <- rownames(tT_1_vs_2)[1:200] 
i <- which(rownames(ex) %in% tT_1_vs_2.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
colSide <- brewer.pal(4,'Set1')[as.numeric(group)]

heatmap.2(ex[i,], scale="row", 
          #labRow=genes_tb$Gene.symbol[i], 
          srtRow = 0, 
          labCol=sml, cexCol = 0.6,
          col=mycol, trace="none", density.info="none", 
          margin=c(2,4), keysize = 1.5, dendrogram="both",
          ColSideColors = colSide )
dev.off()

#######################################################
# Gene set enrichment analyses
# perform the EGSEA analysis set report = TRUE to generate
# the EGSEA interactive report
#######################################################

# dataset annotation
entrezID <- gset@featureData@data[["ENTREZ_GENE_ID"]]
entrezID <- gsub(" ///.*", "", entrezID)

genes_symbol <- gset@featureData@data[["Gene.Symbol"]]
genes_symbol <- gsub(" ///.*", "", genes_symbol)

count = ex
rownames(count) <- entrezID

# remove duplicates
duplicate_pos <- duplicated(entrezID)
na_pos <- is.na(entrezID)
genes_symbol <- genes_symbol[!(duplicate_pos | na_pos)]
count <- count[!(duplicate_pos | na_pos),] # remove genes in correpsonding postions

symbol <- mapIds(org.Hs.eg.db, keys=rownames(count), keytype="ENTREZID", column="SYMBOL")
duplicate_pos_symbol <- duplicated(symbol)
na_pos_symbol <- is.na(symbol)
genes_symbol <- genes_symbol[!(duplicate_pos_symbol | na_pos_symbol)]
count <- count[!(duplicate_pos_symbol | na_pos_symbol),] # remove genes in correpsonding postions

# GSEA by package 'EGSEA'
gs.annots <- buildIdx(entrezIDs = rownames(count), species = "human", gsdb.gsets = "all", 
                      kegg.exclude = c("Metabolism"))
names(gs.annots)

# GSEA with report generated
gsa = egsea.cnt(counts= count, group=group, design=design, 
                contrasts=cont.matrix, gs.annots=gs.annots[c("h")], 
                symbolsMap=cbind(rownames(count), genes_symbol), 
                baseGSEAs=baseMethods, display.top = 20, sort.by="avg.rank", 
                report.dir="./SLE_egsea_avgRank_j-cnt-report", 
                num.threads = 4, report = T)
summary(gsa)
