library(sva)           # For ComBat and ComBat-seq
library(limma)         # For MDS
library(Rtsne)         # For t-SNE
library(umap)          # For UMAP
library(edgeR)         # For cpm

# ============================================
# BATCH CORRECTION METHODS
# ============================================

# 1. ComBat (applied to logCPM)
# --------------------------------------------
logCPM <- cpm(count, prior.count=2, log=TRUE)
rownames(logCPM) <- rownames(count)

# Create model matrix (batch + covariates)
mod <- model.matrix(~1, data=colData)

# Apply ComBat
logCPM_combat <- ComBat(dat=logCPM, 
                        batch=colData$batch, 
                        mod=mod,
                        par.prior=TRUE, 
                        prior.plots=FALSE)

# 2. ComBat-seq (applied to raw counts)
# --------------------------------------------
# ComBat-seq requires integer counts
count_combatseq <- ComBat_seq(counts=as.matrix(count), 
                              batch=colData$batch, 
                              group=NULL)  # or specify group if needed

# Convert ComBat-seq output to logCPM for consistency
logCPM_combatseq <- cpm(count_combatseq, prior.count=2, log=TRUE)


# ============================================
# DIMENSION REDUCTION METHODS
# ============================================

# Function to select top variable genes
select_top_genes <- function(log_matrix, n_genes=2000) {
  top_var_genes <- head(order(matrixStats::rowVars(log_matrix), decreasing=TRUE), n_genes)
  return(log_matrix[top_var_genes, ])
}

logCPM_combat_filtered <- select_top_genes(logCPM_combat, 2000)
logCPM_combatseq_filtered <- select_top_genes(logCPM_combatseq, 2000)

# --------------------------------------------
# 1. PCA 
# --------------------------------------------
perform_pca <- function(log_matrix) {
  pca <- prcomp(t(log_matrix), scale=TRUE)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, digits=2)
  
  # extract 2D coordinates
  pca.coords <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    var_explained_PC1 = pca.var.per[1],
    var_explained_PC2 = pca.var.per[2]
  )
  
  return(list(pca_object = pca, 
              coords = pca.coords,
              var_per = pca.var.per))
}

pca_combat <- perform_pca(logCPM_combat_filtered)
pca_combatseq <- perform_pca(logCPM_combatseq_filtered)

# --------------------------------------------
# 2. MDS (Multidimensional Scaling)
# --------------------------------------------
perform_mds <- function(log_matrix) {
  mds <- plotMDS(log_matrix, plot=FALSE, top=nrow(log_matrix))
  
  # Extract 2D coordinates
  mds.coords <- data.frame(
    sample = colnames(log_matrix),
    MDS1 = mds$x,
    MDS2 = mds$y,
    var_explained_MDS1 = mds$var.explained[1],
    var_explained_MDS2 = mds$var.explained[2]
  )
  
  return(list(mds_object = mds,
              coords = mds.coords))
}

mds_combat <- perform_mds(logCPM_combat_filtered)
mds_combatseq <- perform_mds(logCPM_combatseq_filtered)

# --------------------------------------------
# 3. t-SNE
# --------------------------------------------
perform_tsne <- function(log_matrix, perplexity=30, seed=42) {
  set.seed(seed)
  
  # t-SNE requires samples as rows
  tsne <- Rtsne(t(log_matrix), 
                dims=2, 
                perplexity=perplexity,
                verbose=TRUE,
                max_iter=1000,
                check_duplicates=FALSE)
  
  # extract 2D coordinates
  tsne.coords <- data.frame(
    sample = colnames(log_matrix),
    tSNE1 = tsne$Y[, 1],
    tSNE2 = tsne$Y[, 2]
  )
  
  return(list(tsne_object = tsne,
              coords = tsne.coords))
}

# apply t-SNE (adjust perplexity based on sample size)
# rule of thumb: perplexity should be less than n_samples/3
n_samples <- ncol(logCPM_combat_filtered)
perp <- min(30, floor(n_samples/3))

tsne_combat <- perform_tsne(logCPM_combat_filtered, perplexity=perp)
tsne_combatseq <- perform_tsne(logCPM_combatseq_filtered, perplexity=perp)

# --------------------------------------------
# 4. UMAP
# --------------------------------------------
perform_umap <- function(log_matrix, n_neighbors=15, min_dist=0.1, seed=42) {
  set.seed(seed)
  
  # UMAP configuration
  custom.config <- umap.defaults
  custom.config$n_neighbors <- n_neighbors
  custom.config$min_dist <- min_dist
  custom.config$random_state <- seed
  
  # UMAP requires samples as rows
  umap_result <- umap(t(log_matrix), config=custom.config)
  
  # extract 2D coordinates
  umap.coords <- data.frame(
    sample = colnames(log_matrix),
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )
  
  return(list(umap_object = umap_result,
              coords = umap.coords))
}

umap_combat <- perform_umap(logCPM_combat_filtered)
umap_combatseq <- perform_umap(logCPM_combatseq_filtered)


# ============================================
# Accessing coordinates for plotting
# ============================================

# PCA coordinates (ComBat)
head(pca_combat$coords)

# MDS coordinates (ComBat-seq)
head(mds_combatseq$coords)

# t-SNE coordinates (ComBat)
head(tsne_combat$coords)

# UMAP coordinates (ComBat-seq)
head(umap_combatseq$coords)

# ============================================
# Merge with colData for plotting
# ============================================

# For PCA + ComBat
pca_combat_data <- merge(pca_combat$coords, colData, 
                         by.x="sample", by.y="id_sample")

# For UMAP + ComBat-seq
umap_combatseq_data <- merge(umap_combatseq$coords, colData,
                             by.x="sample", by.y="id_sample")


# 1. PCA + ComBat: pca_combat$coords
# 2. PCA + ComBat-seq: pca_combatseq$coords
# 3. MDS + ComBat: mds_combat$coords
# 4. MDS + ComBat-seq: mds_combatseq$coords
# 5. t-SNE + ComBat: tsne_combat$coords
# 6. t-SNE + ComBat-seq: tsne_combatseq$coords
# 7. UMAP + ComBat: umap_combat$coords
# 8. UMAP + ComBat-seq: umap_combatseq$coords