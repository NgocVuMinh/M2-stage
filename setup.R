
BiocManager::install("biomaRt")
BiocManager::install("edgeR")
BiocManager::install("org.Hs.eg.db")


library(biomaRt)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)



mart <- useMart("ensembl","hsapiens_gene_ensembl")
biotype <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values = gene.ids$ensembl, 
                 mart = mart)
sum(is.na(biotype$gene_biotype)) # 0
gene.ids <- merge(gene.ids, biotype, by.x = "ensembl", by.y="ensembl_gene_id", all.x=T, all.y=T)
