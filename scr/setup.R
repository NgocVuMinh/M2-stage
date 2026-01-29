
BiocManager::install("biomaRt")
BiocManager::install("edgeR")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("sva")
install.packages("Rtsne")
install.packages("umap")
install.packages("gridExtra")


library(biomaRt)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(gridExtra)


theme_custom <- function() {
  theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 16, face = "bold", color = "black", hjust=0.5),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 16, face = "bold", color = "black"),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_blank()
          # legend.title = element_text(size = 12, face = "bold", color = "black")
    )
}




mart <- useMart("ensembl","hsapiens_gene_ensembl")
biotype <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values = gene.ids$ensembl, 
                 mart = mart)
sum(is.na(biotype$gene_biotype)) # 0
gene.ids <- merge(gene.ids, biotype, by.x = "ensembl", by.y="ensembl_gene_id", all.x=T, all.y=T)
