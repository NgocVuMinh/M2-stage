
setwd("workdir/")
load("C:/Users/Equipe_Saulnier/OneDrive - INSTITUT CURIE/env/meta.RData")

library(dplyr)
library(stringr)


# ----- Supplementary table 1 from the HDBR website
dt.supp <- read.csv("data/hdbr/metadata/metadata.txt", sep="\t")
colnames(dt.supp) <- c("id_hdbr", "cs_stage", "id_sample", "tissue", "left_right", 
                  "karyotype", "time_in_transit_mins")
dt.supp <- dt.supp %>%
  mutate(across(c(tissue, karyotype, left_right), tolower))%>%
  mutate(across(everything(), str_trim))
# typo
dt.supp <- dt.supp%>% 
  mutate(tissue = case_when(
    tissue == "hindbrain fragmenet" ~ "hindbrain fragment",
    TRUE ~ tissue))

# ----- Clean tissue morphology
table(dt.supp$tissue)
dt.supp$tissue_clean <- sapply(dt.supp$tissue, function(x) {
  if (grepl("slice", x)) {
    strsplit(x, split=" slice")[[1]][1]
  } else {x}
})
dt.supp <- dt.supp%>% 
  mutate(tissue_clean = case_when(
    tissue_clean == "temporal lobe (hippocampus)" ~ "hippocampus",
    tissue_clean == "cortex" ~ "cerebral cortex", 
    tissue_clean == "basal ganglia" ~ "basal ganglion", 
    tissue_clean == "forebrain- frontal" ~ "forebrain-frontal", 
    tissue_clean == "diencephalon and pituitary" ~ "pituitary and diencephalon",
    TRUE ~ tissue_clean))
table(dt.supp$tissue_clean)
tissue.count.supp <- as.data.frame(table(dt.supp$tissue_clean))
colnames(tissue.count.supp) <- c("tissue", "count")


# ----- Fragment status
dt.supp$fragment <- ifelse(grepl("fragment", dt.supp$tissue_clean), 1, 0)
table(dt.supp$fragment)


# ----- Gender?
table(dt.supp$karyotype)
# 46, XX  46, XY UNKNOWN 
# 332     258      38 
# remove?


# ----- Stage? 
table(dt.supp$cs_stage)
# converting CS to PCW

dt.supp <- dt.supp%>% 
  mutate(pcw_stage = case_when(
    cs_stage == "CS 13" ~ "4 pcw",
    cs_stage %in% c("CS 14", "CS 15") ~ "5 pcw",
    cs_stage %in% c("CS 16", "CS 17") ~ "6 pcw",
    cs_stage %in% c("CS 18", "CS 19") ~ "7 pcw",
    cs_stage %in% c("CS 20", "CS 21", "CS 22", "CS 23", "Late 8 pcw") ~ "8 pcw",
    TRUE ~ cs_stage))

# left right
dt.supp <- dt.supp%>% 
  mutate(left_right = case_when(
    left_right=="" ~ "unknown",
    TRUE ~ left_right))



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

# stage
table(dt.ebi2$Sample.Characteristic.developmental.stage.)
dt.ebi2 <- dt.ebi2%>% 
  mutate(pcw_stage = case_when(
    Sample.Characteristic.developmental.stage. == "carnegie stage 13" ~ "4 post conception weeks",
    Sample.Characteristic.developmental.stage. %in% c("carnegie stage 14", "carnegie stage 15") ~ "5 post conception weeks",
    Sample.Characteristic.developmental.stage. %in% c("carnegie stage 16", "carnegie stage 17") ~ "6 post conception weeks",
    Sample.Characteristic.developmental.stage. %in% c("carnegie stage 18", "carnegie stage 19") ~ "7 post conception weeks",
    Sample.Characteristic.developmental.stage. %in% c("carnegie stage 20", "carnegie stage 21", "carnegie stage 22", "carnegie stage 23", "late 8 post conception weeks") ~ "8 post conception weeks",
    TRUE ~ Sample.Characteristic.developmental.stage.))
dt.ebi2$pcw_stage <- gsub('post conception weeks', 'pcw', dt.ebi2$pcw_stage)

table(dt.ebi2$Sample.Characteristic.sampling.site.)
table(dt.ebi2$Sample.Characteristic.karyotype.)
table(dt.ebi2$Sample.Characteristic.developmental.stage.)

# left right
dt.ebi2 <- dt.ebi2%>% 
  mutate(Sample.Characteristic.sampling.site. = case_when(
    Sample.Characteristic.sampling.site.=="" ~ "unknown",
    TRUE ~ Sample.Characteristic.sampling.site.))




# why 100 fragments, why stomach?
table(dt.ebi2$Sample.Characteristic.organism.part.)

tissue.count <- as.data.frame(table(dt.ebi2$Sample.Characteristic.organism.part.))
colnames(tissue.count) <- c("tissue", "count")

View(dt.ebi2%>%filter(Sample.Characteristic.organism.part.=="basal ganglion"))

d <- dt.ebi2%>%filter(Sample.Characteristic.organism.part.=="brain fragment")
table(d$pcw_stage)


write.table(dt.ebi2%>%filter(!Sample.Characteristic.organism.part. %in% c("choroid plexus", "spinal cord", "stomach")),
            "metadata_ebi_processed.txt", quote = F, sep = "\t")




# ----- consistency check
dt.ebi2.p <- dt.ebi2%>%filter(!Sample.Characteristic.organism.part. %in% c("choroid plexus", "spinal cord", "stomach"))
dt.supp.p <- dt.supp%>%filter(!tissue_clean %in% c("choroid plexus", "spinal cord", "stomach"))
dt.supp.p$id_hdbr <- as.integer(dt.supp.p$id_hdbr)

tissue.count.supp <- as.data.frame(table(dt.supp.p$tissue_clean))
colnames(tissue.count.supp) <- c("tissue", "count")
tissue.count <- as.data.frame(table(dt.ebi2.p$Sample.Characteristic.organism.part.))
colnames(tissue.count) <- c("tissue", "count")


# all patients exist?
length(unique(dt.ebi2.p$Sample.Characteristic.individual.)) # 174
length(unique(dt.supp.p$id_hdbr)) # 171

all(dt.supp.p$id_hdbr %in% dt.ebi2.p$Sample.Characteristic.individual.) 
all(dt.ebi2.p$Sample.Characteristic.individual. %in% dt.supp.p$id_hdbr) # no
which(!dt.ebi2.p$Sample.Characteristic.individual. %in% dt.supp.p$id_hdbr)
View[dt.ebi2.p[c(123, 391, 392),]]


# ID for matching
dt.ebi2.p$ID_EBI <- paste0(dt.ebi2.p$Sample.Characteristic.individual., "_",
                           dt.ebi2.p$Sample.Characteristic.organism.part., "_",
                           dt.ebi2.p$pcw_stage, "_",
                           dt.ebi2.p$Sample.Characteristic.sampling.site.)
dt.supp.p$ID_SUPP <- paste0(dt.supp.p$id_hdbr, "_",
                            dt.supp.p$tissue_clean, "_",
                            dt.supp.p$pcw_stage, "_",
                            dt.supp.p$left_right)

length(dt.ebi2.p$ID_EBI)
length(unique(dt.ebi2.p$ID_EBI))
length(dt.supp.p$ID_SUPP)
length(unique(dt.supp.p$ID_SUPP))




# count data







