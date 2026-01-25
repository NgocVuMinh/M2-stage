
setwd("workdir/")

library(dplyr)


# ----- Supplementary table 1 from the HDBR website
dt.supp <- read.csv("data/hdbr/metadata/metadata.txt", sep="\t")
colnames(dt.supp) <- c("id_hdbr", "cs_stage", "id_sample", "tissue", "left_right", 
                  "karyotype", "time_in_transit_mins")

# ----- Clean tissue morphology
table(dt.supp$tissue)
dt.supp$tissue_clean <- sapply(dt.supp$tissue, function(x) {
  if (grepl("slice", x)) {
    strsplit(x, split=" slice")[[1]][1]
  } else {x}
})
table(dt.supp$tissue_clean)


# ----- Fragment status
dt.supp$fragment <- ifelse(grepl("fragment|fragmenet", dt.supp$tissue_clean), 1, 0)
table(dt.supp$fragment)


# ----- Gender?
table(dt.supp$karyotype)
# 46, XX  46, XY UNKNOWN 
# 332     258      38 
# remove?


# ----- Stage? 
table(dt.supp$cs_stage)



# ----- EBI metadata
dt.ebi <- read.csv("data/hdbr/metadata/E-MTAB-4840-experiment-design.tsv", sep="\t")
colnames(dt.ebi)
View(dt.ebi[,c(16:19)])

dt.ebi2 <- dt.ebi[, !grepl("Ontology", names(dt.ebi))]
colnames(dt.ebi2)


# ----- Duplicated columns --> remove
identical(dt.ebi2$Sample.Characteristic.developmental.stage., dt.ebi2$Factor.Value.developmental.stage.)
identical(dt.ebi2$Sample.Characteristic.organism.part., dt.ebi2$Factor.Value.organism.part.)
dt.ebi2 <- dt.ebi2[,c(1:5, 7,8, 11)]

# all lowercase
dt.ebi2 <- dt.ebi2 %>%
  mutate(across(c(Sample.Characteristic.developmental.stage., Sample.Characteristic.karyotype.,
                  Sample.Characteristic.organism.part., Sample.Characteristic.sampling.site.), tolower))


table(dt.ebi2$Sample.Characteristic.sampling.site.)
table(dt.ebi2$Sample.Characteristic.karyotype.)
table(dt.ebi2$Sample.Characteristic.developmental.stage.)

# this is important: why 100 fragments, why stomach?
table(dt.ebi2$Sample.Characteristic.organism.part.)










