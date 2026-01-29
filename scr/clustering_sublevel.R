

TOP_VAR=250

logCPM_top <- select_top_genes(logCPM, TOP_VAR)
logCPM_combat_top <- select_top_genes(logCPM_combat, TOP_VAR)


# ----- top var genes
# FORE_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
# MID_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
# HIND_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="hindbrain")%>%pull(id_sample)], TOP_VAR)
# dim(FORE_logCPM_top)
# dim(MID_logCPM_top)
# dim(HIND_logCPM_top)

FORE_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
MID_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
HIND_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="hindbrain" & !tissue_level2%in%c("pons", "whole hindbrain"))%>%pull(id_sample)], TOP_VAR)



input="FORE"
input="MID"
input="HIND"

# pca_ <- perform_pca(get(paste0(input, "_logCPM_top")))
pca_combat <- perform_pca(get(paste0(input, "_logCPM_combat_top")))

# tsne_ <- perform_tsne(get(paste0(input, "_logCPM_top")), perplexity=30)
tsne_combat <- perform_tsne(get(paste0(input, "_logCPM_combat_top")), perplexity=15)

# umap_ <- perform_umap(get(paste0(input, "_logCPM_top")))
umap_combat <- perform_umap(get(paste0(input, "_logCPM_combat_top")))



# ============================================
# Plotting
# ============================================

clus_input <- merge(pca_combat$coords, colData, by.x="sample", by.y="id_sample")
# clus_input_ <- merge(pca_$coords, colData, by.x="sample", by.y="id_sample")

clus_input <- merge(tsne_combat$coords, colData, by.x="sample", by.y="id_sample")
# clus_input_ <- merge(tsne_$coords, colData, by.x="sample", by.y="id_sample")

clus_input <- merge(umap_combat$coords, colData, by.x="sample", by.y="id_sample")
# clus_input_ <- merge(umap_$coords, colData, by.x="sample", by.y="id_sample")


FILL1="stage"
COLORS1="colors_stage"
FILL2="tissue_level2"
COLORS2="colors_tissue2"
FILL3="tissue"
COLORS3="colors_tissue"
SHAPE= "tissue_level2"
POINT_SIZE=3
X="UMAP1"
Y="UMAP2"

p1 <- ggplot(clus_input, aes(x = .data[[X]], y = .data[[Y]], color = .data[[FILL1]], shape = .data[[SHAPE]])) +
  geom_point(size = POINT_SIZE, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(16, 17, 15, 0, 1, 2, 3, 4, 5, 6, 8)) +
  scale_color_manual(values = get(COLORS1)) +
  labs(title = " ") +
  theme_custom() # + theme(legend.position = "bottom")
p2 <- ggplot(clus_input, aes(x = .data[[X]], y = .data[[Y]], color = .data[[FILL2]], shape = .data[[SHAPE]])) +
  geom_point(size = POINT_SIZE, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 8)) +
  scale_color_manual(values = get(COLORS2)) +
  labs(title = " ") +
  theme_custom() + 
  theme(legend.position = "bottom")
p3 <- ggplot(clus_input, aes(x = .data[[X]], y = .data[[Y]], color = .data[[FILL3]], shape = .data[[SHAPE]])) +
  geom_point(size = POINT_SIZE, stroke = 0.8, alpha = 1) +  
  scale_shape_manual(values = c(15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 8)) +
  scale_color_manual(values = get(COLORS3)) +
  labs(title = " ") +
  theme_custom() + 
  theme(legend.position = "bottom")
grid.arrange(p1, p2, p3, ncol = 3)


plot_ly(data = clus_input,
        #x = ~PC1, y = ~PC2, z = ~PC3, 
        #x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, 
        x = ~UMAP2, y = ~UMAP1, z = ~UMAP3,
        color = ~stage,
        colors = colors_stage,
        symbol = ~tissue_level2,
        symbols = c("circle", "diamond", "square", "x", "cross", "triangle-down"),
        text = ~stage,
        type = "scatter3d", mode = "markers",
        marker = list(size = 4, line = list(color = "black", width = 1), opacity = 1)) %>%
  layout(title = " ",
         # scene = list(xaxis = list(title = paste0("PC1: ", clus_input$var_explained_PC1[1], "%")),
         #              yaxis = list(title = paste0("PC2: ", clus_input$var_explained_PC2[1], "%")),
         #              zaxis = list(title = paste0("PC3: ", clus_input$var_explained_PC3[1], "%"))),
         legend = list(font = list(size = 10),
                       itemsizing = "constant")
  )










