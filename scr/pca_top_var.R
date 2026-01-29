
# ============================================
# PCA % variance explained FOR BATCH CORRECTION:
# - ComBat on 464 samples
# - ComBat on 450 samples (removed those from low-count batches)
# - ComBat-seq 450 samples 
# ============================================

logCPM <- cpm(count, prior.count=2, log=TRUE)

# model matrix (covariates only)
mod <- model.matrix(~ tissue + stage + karyotype, data=colData)

# ComBat on all 464 samples
logCPM_combat <- ComBat(dat = logCPM,
                        batch = colData$batch,
                        mod = mod,
                        par.prior = TRUE, prior.plots = FALSE)
# ComBat on 450 samples
logCPM_combat_filt <- ComBat(dat = as.matrix(logCPM[, !colnames(logCPM) %in% batch.low.samples]),
                             batch = colData$batch[!colnames(logCPM) %in% batch.low.samples],
                             mod = as.matrix(mod[!colnames(logCPM) %in% batch.low.samples, , drop=FALSE]),
                             par.prior  = TRUE, prior.plots = FALSE)
# ComBat-seq on 450 samples - long as hell omg
count_combatseq <- ComBat_seq(counts= as.matrix(count[, !colnames(count) %in% batch.low.samples]),
                              batch = colData$batch[!colnames(count) %in% batch.low.samples],
                              covar_mod = as.matrix(mod[!colnames(logCPM) %in% batch.low.samples, , drop=FALSE]))
logCPM_combatseq_filt <- cpm(count_combatseq, prior.count=2, log=TRUE)



PCA_COMBAT_COMP <- list()

for (TOP_VAR in c(250, 500, 1000, 2000)) {
  
  logCPM_top <- select_top_genes(logCPM, TOP_VAR)
  logCPM_combat_top <- select_top_genes(logCPM_combat, TOP_VAR)
  logCPM_combat_filt_top <- select_top_genes(logCPM_combat_filt, TOP_VAR)
  logCPM_combatseq_filt_top <- select_top_genes(logCPM_combatseq_filt, TOP_VAR)
  
  for (i in c("logCPM_top", "logCPM_combat_top", "logCPM_combat_filt_top", "logCPM_combatseq_filt_top")) 
    {
      pca_combat <- perform_pca(get(i))
      clus_input <- merge(pca_combat$coords, colData, by.x="sample", by.y="id_sample")
      
      # sum of %var explained by 2PC, 3PC
      sum2PC <- clus_input$var_explained_PC1[1] + clus_input$var_explained_PC2[1]
      sum3PC <- sum2PC + clus_input$var_explained_PC3[1]
      
      temp_df <- data.frame(method = i,
                            top_var = TOP_VAR,
                            sumPC2 = sum2PC,
                            sumPC3 = sum3PC,
                            stringsAsFactors = FALSE)
      PCA_COMBAT_COMP[[length(PCA_COMBAT_COMP) + 1]] <- temp_df
    }
  }
PCA_COMBAT_COMP <- do.call(rbind, PCA_COMBAT_COMP)
PCA_COMBAT_COMP <- PCA_COMBAT_COMP%>%
  mutate(top_var = factor(as.character(top_var), levels=c("250", "500", "1000", "2000")),
         method = case_when(
           method == "logCPM_top" ~ "No corrections",
           method == "logCPM_combat_top" ~ "ComBat",
           method == "logCPM_combat_filt_top" ~ "ComBat (filtered small batches)",
           method == "logCPM_combatseq_filt_top" ~ "ComBat-seq (filtered small batches)"
         ))%>%
  mutate(method = factor(method, levels=c("No corrections", "ComBat", "ComBat (filtered small batches)", "ComBat-seq (filtered small batches)")))


ggplot(PCA_COMBAT_COMP, aes(x = top_var, y = sumPC3, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, padding=0.0), 
           color = "black", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = round(sumPC3, 1), group = method), 
            position = position_dodge2(width = 0.8, padding=0.0), 
            vjust = -0.5, size = 4) + 
  # facet_wrap(~Tissue, scales = "free_y") +
  scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(ylim = c(50, 70)) + 
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  scale_fill_manual(values = c("#6caed6", "#cbaa7a", "grey", "green")) +
  labs(x = "Number of top variable genes",
       y = "Total % variance explained by PC1-3") +
  theme_custom()









PCA_TOP_VARS <- list()

for (TOP_VAR in c(250, 500, 1000, 2000)) {
  
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
      TRUE~ "All regions"
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
   aes(x = Tissue, y = sumPC3, fill = factor(TopVar))) +
  geom_bar(stat = "identity", 
   position = position_dodge2(width = 0.9, padding = 0), 
   color = "black", linewidth = 0.3, width = 0.85) +
  geom_text(aes(label = round(sumPC3, 1)), 
position = position_dodge2(width = 0.9, padding = 0), 
vjust = -0.5, size = 5) + 
  coord_cartesian(ylim = c(40, 75)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_brewer(palette = "Blues", direction=-1) +
  labs(x = " ",
   y = "Total % variance\nexplained by PC1-3",
   fill = "Number of top\nvariable genes") +
  theme_custom()







