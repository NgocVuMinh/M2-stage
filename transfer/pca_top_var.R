

PCA_TOP_VARS <- data.frame(combination = c(),
                           sumPC2 = c(),
                           sumPC3 = c(),
                           stringsAsFactors = F)


for (TOP_VAR in c(500, 1000, 2000)) {
  
  logCPM_top <- select_top_genes(logCPM, TOP_VAR)
  logCPM_combat_top <- select_top_genes(logCPM_combat, TOP_VAR)
  
  FORE_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
  FORE_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
  
  MID_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
  MID_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
  
  HIND_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="hindbrain")%>%pull(id_sample)], TOP_VAR)
  HIND_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="hindbrain")%>%pull(id_sample)], TOP_VAR)
  
  for (input in c("FORE_", "MID_", "HIND_" , "")) {
    for (i in c("logCPM_combat_top", "logCPM_top")) {
      pca_combat <- perform_pca(get(paste0(input, i)))
      clus_input <- merge(pca_combat$coords, colData, by.x="sample", by.y="id_sample")
      
      # sum of %var explained by 2PC, 3PC
      sum2PC <- clus_input$var_explained_PC1[1] + clus_input$var_explained_PC2[1]
      sum3PC <- sum2PC + clus_input$var_explained_PC3[1]
      
      temp_df <- data.frame(
        combination = paste0(var_name, "_", TOP_VAR),
        sumPC2 = sum2PC,
        sumPC3 = sum3PC,
        stringsAsFactors = FALSE
      )
      PCA_TOP_VARS[[length(PCA_TOP_VARS) + 1]] <- temp_df
      
      combination <- paste0(input, i, "_", TOP_VAR)
      PCA_TOP_VARS <- rbind(PCA_TOP_VARS,
                            data.frame(combination = combination,
                                       sumPC2 = sum2PC,
                                       sumPC3 = sum3PC,
                                       stringsAsFactors = F))
    }
    
  }
}



PCA_TOP_VARS <- list()

for (TOP_VAR in c(250, 500, 750, 1000, 2000)) {
  
  logCPM_top <- select_top_genes(logCPM, TOP_VAR)
  logCPM_combat_top <- select_top_genes(logCPM_combat, TOP_VAR)
  
  FORE_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
  FORE_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="forebrain")%>%pull(id_sample)], TOP_VAR)
  
  MID_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
  MID_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="midbrain")%>%pull(id_sample)], TOP_VAR)
  
  HIND_logCPM_top <- select_top_genes(logCPM[,colData%>%filter(tissue_level1=="hindbrain")%>%pull(id_sample)], TOP_VAR)
  HIND_logCPM_combat_top <- select_top_genes(logCPM_combat[,colData%>%filter(tissue_level1=="hindbrain")%>%pull(id_sample)], TOP_VAR)
  
  for (region in c("FORE_", "MID_", "HIND_", "")) {
    for (method in c("logCPM_combat_top", "logCPM_top")) {
      
      var_name <- paste0(region, method)
      pca_combat <- perform_pca(get(var_name))
      clus_input <- merge(pca_combat$coords, colData, by.x = "sample", by.y = "id_sample")
      
      v1 <- clus_input$var_explained_PC1[1]
      v2 <- clus_input$var_explained_PC2[1]
      v3 <- clus_input$var_explained_PC3[1]
      
      sum2PC <- v1 + v2
      sum3PC <- sum2PC + v3
      
      temp_df <- data.frame(
        combination = paste0(var_name, "_", TOP_VAR),
        sumPC2 = sum2PC,
        sumPC3 = sum3PC,
        stringsAsFactors = FALSE
      )
      PCA_TOP_VARS[[length(PCA_TOP_VARS) + 1]] <- temp_df
    }
  }
}

PCA_TOP_VARS <- do.call(rbind, PCA_TOP_VARS)

PCA_TOP_VARS <- PCA_TOP_VARS %>%
  mutate(
    TopVar = as.character(str_extract(combination, "\\d+$")),
    BatchCorrected = ifelse(str_detect(combination, "combat"), "ComBat", "Raw"),
    Tissue = case_when(
      str_detect(combination, "FORE") ~ "Forebrain",
      str_detect(combination, "MID")  ~ "Midbrain",
      str_detect(combination, "HIND") ~ "Hindbrain",
      TRUE                            ~ "All regions"
    )
  ) %>%
  mutate(TopVar = factor(TopVar, levels=c("250", "500", "750", "1000", "2000")))


ggplot(PCA_TOP_VARS, aes(x = TopVar, y = sumPC3, fill = BatchCorrected)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, padding=0.0), 
           color = "black", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = round(sumPC3, 1), group = BatchCorrected), 
            position = position_dodge2(width = 0.9, padding=0.0), 
            vjust = -0.5, size = 4) + 
  facet_wrap(~Tissue, scales = "free_y") +
  scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(ylim = c(50, 75)) + 
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  scale_fill_manual(values = c("#6caed6", "#cbaa7a", "grey")) +
  labs(x = "Number of top variable genes",
       y = "Total % variance explained by PC1-3") +
  theme_custom()



ggplot(PCA_TOP_VARS %>% filter(BatchCorrected == "ComBat"), 
       aes(x = Tissue, y = sumPC2, fill = factor(TopVar))) +
  geom_bar(stat = "identity", 
           position = position_dodge2(width = 0.9, padding = 0), 
           color = "black", linewidth = 0.3, width = 0.85) +
  geom_text(aes(label = round(sumPC2, 1)), 
            position = position_dodge2(width = 0.9, padding = 0), 
            vjust = -0.5, size = 5) + 
  coord_cartesian(ylim = c(40, 75)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_brewer(palette = "Blues", direction=-1) +
  labs(x = " ",
       y = "Total % variance\nexplained by PC1-2",
       fill = "Number of top\nvariable genes") +
  theme_custom()







