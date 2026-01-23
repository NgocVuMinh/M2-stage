
library(dplyr)

dt <- read.csv("metadata.txt", sep="\t")
colnames(dt) <- c("id_hdbr", "cs_stage", "id_sample", "tissue", "left_right", "karyotype", "time_in_transit_mins")

table(dt$tissue)

dt$tissue_clean <- sapply(dt$tissue, function(x) {
  if (grepl("slice", x)) {
    strsplit(x, split=" slice")[[1]][1]
  } else {x}
})
table(dt$tissue_clean)


dt$fragment <- ifelse(grepl("fragment|fragmenet", dt$tissue_clean), 1, 0)
table(dt$fragment)




table(dt$karyotype)
# 46, XX  46, XY UNKNOWN 
# 332     258      38 
# remove?

table(dt$cs_stage)
