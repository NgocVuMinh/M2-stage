

library(biomaRt)
library(plotly)


count <- count0[,sel$Run]

colData <- sel
rownames(colData) <- colData$Run
colnames(colData)
colData <- colData[,-c(3,8,10)]
colnames(colData) <- c("id_sample", "batch", "id_individual", "karyotype", "tissue",
                       "site", "stage", "multiple_batches")


colData <- colData %>%
  mutate(batch = str_remove(batch, "sequencing batch\\s*:\\s*")) %>%
  group_by(id_individual, batch) %>% 
  mutate(
    batch = mapply(function(s, i) {
      # a comma or dot followed specifically by a space
      if (str_detect(s, "[,.]\\s+")) {
        str_split_i(s, "[,.]\\s+", i)
      } else {
        s
      }
    }, batch, row_number())
  ) %>%
  ungroup()

length(unique(colData$batch)) # 29
batch.counts <- as.data.frame(table(colData$batch))


colData$stage <- factor(colData$stage, 
                        levels=c("4 pcw", "5 pcw", "6 pcw", "7 pcw", "8 pcw","9 pcw",
                                 "10 pcw", "11 pcw", "12 pcw", "13 pcw", "14 pcw", "15 pcw", "16 pcw",
                                 "17 pcw"))
colData$site <- factor(colData$site, 
                        levels=c("right", "left", "unknown"))


gene.ids <- data.frame(ensemble = rownames(count),
                       entrez = NA,
                       symbol = NA,
                       biotype = NA,
                       stringsAsFactors = F)
rownames(gene.ids) <- gene.ids$ensemble

# -------------- get gene info
keytypes(org.Hs.eg.db)
gene.ids$symbol <- mapIds(org.Hs.eg.db, keys=rownames(gene.ids),keytype = "ENSEMBL", column = "SYMBOL", multiVals="first")
gene.ids$name <- mapIds(org.Hs.eg.db, keys=rownames(gene.ids), keytype = "ENSEMBL", column = "GENENAME", multiVals="first")
gene.ids$entrez <- mapIds(org.Hs.eg.db, keys=rownames(gene.ids), keytype = "ENSEMBL", column = "ENTREZID", multiVals="first")
gene.ids$biotype <- mapIds(org.Hs.eg.db, keys=rownames(gene.ids), keytype = "ENSEMBL", column = "GENETYPE", multiVals="first")
sum(is.na(gene.ids$entrez)) # 26130




all(rownames(colData) == colnames(counts)) # true

#  Keep genes expressed >10 in at least 80% of samples
keep.perc <- 0.8
keep <- rowSums(count > 10) >= ceiling(keep.perc * ncol(count))
sum(keep) # 16535

count <- count[keep,]

logCPM <- cpm(count, prior.count=2, log=TRUE)
rownames(logCPM) <- rownames(count)
# hist(logCPM[,1])


logCPM <- logCPM[head(order(matrixStats::rowVars(logCPM), decreasing=T), 2000), ]
pca <- prcomp(t(logCPM), scale=T)


# -------------- batch effect correction
log.cpm.corrected <- ComBat(dat=logCPM, batch=colData$batch, mod=NULL)


pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, digits=2)
pca.dat <- data.frame(pca$x)

pca.dat$tissue <- colData$tissue
pca.dat$stage <- colData$stage
pca.dat$site <- colData$site
pca.dat$batch <- colData$batch
pca.dat$karyotype <- colData$karyotype
 

ggplot(pca.dat, aes(x = PC1, y = PC2, fill=batch)) + # , shape=Time
  geom_point(size = 2, stroke = 0.6, alpha=1, shape = 21) +  
  # geom_mark_ellipse(aes(fill= Group), show.legend=F) + 
  labs(x = paste0('PC1: ', pca.var.per[1], '%'),
       y = paste0('PC2: ', pca.var.per[2], '%'),
       title="") +
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21), color = "black", stroke = 0.6))) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        #panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold", color = "black", hjust=0.5),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank()
        # legend.title = element_text(size = 12, face = "bold", color = "black")
  )




plot_ly(data = pca.dat, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~batch,
        # colors = c(),
        # = ~Time,
        text = ~batch,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 3, line = list(color = "black", width = 1), opacity = 1)) %>%
  layout(title = " ",
         scene = list(xaxis = list(title = paste0("PC1: ", pca.var.per[1], "%")),
                      yaxis = list(title = paste0("PC2: ", pca.var.per[2], "%")),
                      zaxis = list(title = paste0("PC3: ", pca.var.per[3], "%"))
         )
  )

















