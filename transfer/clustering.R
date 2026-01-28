
load("C:/Users/Equipe_Saulnier/OneDrive - INSTITUT CURIE/env/meta.RData")

library(sva)           # For ComBat and ComBat-seq
library(limma)         # For MDS
library(Rtsne)         # For t-SNE
library(umap)          # For UMAP
library(edgeR)         # For cpm



colData <- as.data.frame(colData)
rownames(colData) <- colData$id_sample
logCPM <- cpm(count, prior.count=2, log=TRUE)
rownames(logCPM) <- rownames(count)
# hist(logCPM[,1])


# ============================================
# BATCH CORRECTION METHODS
# ============================================

# 1. ComBat
# --------------------------------------------
logCPM <- cpm(count, prior.count=2, log=TRUE)
rownames(logCPM) <- rownames(count)

# model matrix (covariates?)
mod <- model.matrix(~1+tissue + stage + karyotype, data=colData[,-1])

logCPM_combat <- ComBat(dat=logCPM, 
                        batch=colData$batch, 
                        mod=mod, 
                        par.prior=T, 
                        prior.plots=F)
# !!!!! one batch FC2014_053 has only one sample, setting mean.only=TRUE



# 2. ComBat-seq
# --------------------------------------------
count_combatseq <- ComBat_seq(counts=as.matrix(count), 
                              batch=colData$batch, 
                              covar_mod = mod)
# ComBat-seq doesn't support 1 sample per batch yet
logCPM_combatseq <- cpm(count_combatseq, prior.count=2, log=TRUE)







# ============================================
# DIMENSION REDUCTION METHODS
# ============================================

# select top variable genes: 2000 vs 5000?
select_top_genes <- function(log_matrix, n_genes=2000) {
  top_var_genes <- head(order(matrixStats::rowVars(log_matrix), decreasing=TRUE), n_genes)
  return(log_matrix[top_var_genes, ])
}

logCPM_top <- select_top_genes(logCPM, 2000)
logCPM_combat_top <- select_top_genes(logCPM_combat, 2000)
logCPM_combat_top <- select_top_genes(logCPM_combat, 5000)
# logCPM_combatseq_top <- select_top_genes(logCPM_combatseq, 2000)



# --------------------------------------------
# 1. PCA 
# --------------------------------------------
perform_pca <- function(log_matrix) {
  pca <- prcomp(t(log_matrix), scale=T)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, digits=2)
  
  # extract 2D coordinates
  pca.coords <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    var_explained_PC1 = pca.var.per[1],
    var_explained_PC2 = pca.var.per[2],
    var_explained_PC3 = pca.var.per[3]
  )
  
  return(list(pca_object = pca, 
              coords = pca.coords,
              var_per = pca.var.per))
}

pca_ <- perform_pca(logCPM_top)
pca_combat <- perform_pca(logCPM_combat_top)
# pca_combatseq <- perform_pca(logCPM_combatseq_top)



# --------------------------------------------
# 2. MDS
# --------------------------------------------
perform_mds <- function(log_matrix) {
  mds <- plotMDS(log_matrix, plot=FALSE, top=nrow(log_matrix), dim.plot = c(1,3))
  
  # extract 2D coordinates
  mds.coords <- data.frame(
    sample = colnames(log_matrix),
    MDS1 = mds$x,
    MDS2 = mds$y,
    var_explained_MDS1 = mds$var.explained[1],
    var_explained_MDS2 = mds$var.explained[3]
  )
  
  return(list(mds_object = mds,
              coords = mds.coords))
}

mds_ <- perform_mds(logCPM_top)
mds_combat <- perform_mds(logCPM_combat_top)
# mds_combatseq <- perform_mds(logCPM_combatseq_top)




# --------------------------------------------
# 3. t-SNE
# --------------------------------------------
perform_tsne <- function(log_matrix, perplexity=30, seed=42) {
  set.seed(seed)
  
  # t-SNE requires samples as rows
  tsne <- Rtsne(t(log_matrix), 
                dims=3, 
                perplexity=perplexity,
                theta=0.1,
                verbose=TRUE,
                max_iter=2000,
                check_duplicates=FALSE)
  
  # extract 2D coordinates
  tsne.coords <- data.frame(
    sample = colnames(log_matrix),
    tSNE1 = tsne$Y[, 1],
    tSNE2 = tsne$Y[, 2],
    tSNE3 = tsne$Y[, 3]
  )
  
  return(list(tsne_object = tsne,
              coords = tsne.coords))
}

# apply t-SNE (adjust perplexity based on sample size)
# rule of thumb: perplexity should be less than n_samples/3
n_samples <- ncol(logCPM_combat_top)
perp <- min(30, floor(n_samples/3))

tsne_ <- perform_tsne(logCPM_top, perplexity=30)
tsne_combat <- perform_tsne(logCPM_combat_top, perplexity=30)
# tsne_combatseq <- perform_tsne(logCPM_combatseq_top, perplexity=perp)



# --------------------------------------------
# 4. UMAP
# --------------------------------------------
perform_umap <- function(log_matrix, n_neighbors=15, min_dist=0.1, n_components=3, seed=42) {
  set.seed(seed)
  
  # UMAP configuration
  custom.config <- umap.defaults
  custom.config$n_neighbors <- n_neighbors
  custom.config$min_dist <- min_dist
  custom.config$random_state <- seed
  custom.config$n_components <- n_components
  
  # UMAP requires samples as rows
  umap_result <- umap(t(log_matrix), config=custom.config)
  
  # extract 2D coordinates
  umap.coords <- data.frame(
    sample = colnames(log_matrix),
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    UMAP3 = umap_result$layout[, 3]
  )
  
  return(list(umap_object = umap_result,
              coords = umap.coords))
}

umap_ <- perform_umap(logCPM_top)
umap_combat <- perform_umap(logCPM_combat_top)
# umap_combatseq <- perform_umap(logCPM_combatseq_top)



# ============================================
# Accessing coordinates for plotting
# ============================================

head(pca_combat$coords)
head(mds_combat$coords)
head(tsne_combat$coords)
head(umap_combat$coords)


# ============================================
# Plotting
# ============================================

clus_input <- merge(pca_combat$coords, colData, 
                         by.x="sample", by.y="id_sample")
clus_input_ <- merge(pca_$coords, colData,
                     by.x="sample", by.y="id_sample")


clus_input <- merge(mds_combat$coords, colData,
                          by.x="sample", by.y="id_sample")
clus_input_ <- merge(mds_$coords, colData,
                     by.x="sample", by.y="id_sample")


clus_input <- merge(tsne_combat$coords, colData, 
                         by.x="sample", by.y="id_sample")
clus_input_ <- merge(tsne_$coords, colData,
                     by.x="sample", by.y="id_sample")


clus_input <- merge(umap_combat$coords, colData,
                          by.x="sample", by.y="id_sample")
clus_input_ <- merge(umap_$coords, colData,
                    by.x="sample", by.y="id_sample")


# ------------------------
FILL="tissue_level1"
SHAPE=NULL
X="UMAP2"
Y="UMAP1"
COLORS="colors_tissue1"
ggplot(clus_input, aes(x = get(X), y = get(Y), fill=get(FILL))) +
  geom_point(size = 3, stroke = 0.6, alpha=1, shape=21) +  
  labs(x = X,#paste0('PC1: ', clus_input$var_explained_PC1[1], '%'),
       y = Y, #paste0('PC2: ', clus_input$var_explained_PC2[1], '%'),
       title="ComBat") +
  scale_shape_manual(values = c(21:25, 0:6)) +
  #scale_fill_manual(values = get(COLORS)) +
  guides(fill = guide_legend(override.aes = list(shape = c(21), color = "black", stroke = 0.6))) +
  theme_custom()+ theme(legend.position ="bottom")


# ------------------------ no stage
FILL="batch"
SHAPE="tissue_level1"
X="UMAP1"
Y="UMAP2"
COLORS="colors_stage"
p1 <- ggplot(clus_input_, aes(x = .data[[X]], y = .data[[Y]], 
                              color = .data[[FILL]], shape = .data[[SHAPE]])) +
  geom_point(size = 3, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(16, 17, 18, 15, 0, 1, 2, 3, 4, 5, 6, 8)) +
  #scale_color_manual(values = get(COLORS)) +
  labs(title = "No corrections") +
  # labs(x = paste0('PC1: ', clus_input_$var_explained_PC1[1], '%'),
  #      y = paste0('PC2: ', clus_input_$var_explained_PC2[1], '%'),
  #      title="No corrections") +
  # labs(x = paste0('D1: ', round(clus_input_$var_explained_MDS1[1]*100, 2), '%'),
  #      y = paste0('D3: ', round(clus_input_$var_explained_MDS2[1]*100, 2), '%'),
  #      title="No corrections") +
  theme_custom() + 
  theme(legend.position = "bottom")
p2 <- ggplot(clus_input, aes(x = .data[[X]], y = .data[[Y]], 
                             color = .data[[FILL]], shape = .data[[SHAPE]])) +
  geom_point(size = 3, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(16, 17, 18, 15, 0, 1, 2, 3, 4, 5, 6, 8)) +
  #scale_color_manual(values = get(COLORS)) +
  labs(title = "ComBat") +
  # labs(x = paste0('PC1: ', clus_input$var_explained_PC1[1], '%'),
  #      y = paste0('PC2: ', clus_input$var_explained_PC2[1], '%'),
  #      title="ComBat") +
  # labs(x = paste0('D1: ', round(clus_input$var_explained_MDS1[1]*100, 2), '%'),
  #      y = paste0('D3: ', round(clus_input$var_explained_MDS2[1]*100, 2), '%'),
  #      title="ComBat") +
  theme_custom() + 
  theme(legend.position = "bottom")
grid.arrange(p1, p2, ncol = 2)


# ------------------------ with stage
FILL="stage"
SHAPE="tissue_level1"
X="UMAP2"
Y="UMAP1"
COLORS="colors_stage"
p1 <- ggplot(clus_input_, aes(x = .data[[X]], y = .data[[Y]], color = .data[[FILL]], shape = .data[[SHAPE]])) +
  geom_point(size = 3, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(16, 17, 18, 15, 0, 1, 2, 3, 4, 5, 6, 8)) +
  scale_color_manual(values = get(COLORS)) +
  labs(title = "No corrections") +
  theme_custom() + 
  theme(legend.position = "bottom")
p2 <- ggplot(clus_input, aes(x = .data[[X]], y = .data[[Y]], color = .data[[FILL]], shape = .data[[SHAPE]])) +
  geom_point(size = 3, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 8)) +
  scale_color_manual(values = get(COLORS)) +
  labs(title = "ComBat") +
  theme_custom() + 
  theme(legend.position = "bottom")
grid.arrange(p1, p2, ncol = 2)




plot_ly(data = clus_input, 
        #x = ~PC1, y = ~PC2, z = ~PC3, 
        #x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, 
        x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
        color = ~stage, #batch,
        colors = colors_stage, #colors_tissue2,
        symbol = ~tissue_level1,
        symbols = c("circle", "square", "diamond"),
        # = ~Time,
        text = ~batch,
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 4, line = list(color = "black", width = 1), opacity = 1)) %>%
  layout(title = " ",
         # scene = list(xaxis = list(title = paste0("PC1: ", clus_input$var_explained_PC1[1], "%")),
         #              yaxis = list(title = paste0("PC2: ", clus_input$var_explained_PC2[1], "%")),
         #              zaxis = list(title = paste0("PC3: ", clus_input$var_explained_PC3[1], "%"))),
         legend = list(font = list(size = 10),
                       itemsizing = "constant")
  )






# 1. PCA + ComBat: pca_combat$coords
# 2. PCA + ComBat-seq: pca_combatseq$coords
# 3. MDS + ComBat: mds_combat$coords
# 4. MDS + ComBat-seq: mds_combatseq$coords
# 5. t-SNE + ComBat: tsne_combat$coords
# 6. t-SNE + ComBat-seq: tsne_combatseq$coords
# 7. UMAP + ComBat: umap_combat$coords
# 8. UMAP + ComBat-seq: umap_combatseq$coords