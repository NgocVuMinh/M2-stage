
install.packages("paletteer")
library(paletteer)


as.character(paletteer_d("colorBlindness::SteppedSequential5Steps"))

c(list(c("a", "b"),
       c(1, 2)))


colors_tissue <- setNames(c("#E57E7EFF", "red3", "#FFB2B2FF", 
                            "#99540FFF", "#E5B17EFF", "#FFD8B2FF",
                            "#465e17", "#85B22CFF", "#E5FFB2FF", 
                            "#6551CCFF", "#5ba3c9", "#c9e6f5"),
                          unique(colData$tissue))
colors_tissue1 <- setNames(c("#E57E7EFF", "#85B22CFF", "#5ba3c9"),
                          unique(colData$tissue_level1))
colors_tissue2 <- setNames(c("#E57E7EFF", "red3", "#85B22CFF", "#E5FFB2FF", 
                             "#5ba3c9", "#c9e6f5", "#E5B17EFF", "#BFB2FFFF"),
                           unique(colData$tissue_level2))


colors_stage <- setNames(as.character(paletteer_c("ggthemes::Blue", 14)),
                         levels(colData$stage))
