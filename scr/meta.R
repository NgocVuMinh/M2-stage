
setwd("workdir/")
load("C:/Users/Equipe_Saulnier/OneDrive - INSTITUT CURIE/env/meta.RData")

library(dplyr)
library(stringr)
library(ggplot2)


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
    tissue_clean == "forebrain-frontal" ~ "forebrain",
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

all(dt.supp.p$id_hdbr %in% dt.ebi2.p$Sample.Characteristic.individual.) # yes 
all(dt.ebi2.p$Sample.Characteristic.individual. %in% dt.supp.p$id_hdbr) # no
which(!dt.ebi2.p$Sample.Characteristic.individual. %in% dt.supp.p$id_hdbr)
View(dt.ebi2.p[c(123, 391, 392),]) # 3 more samples than supp, check EBI if there are authors' notes


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
count0 <- read.table("C:/Users/Equipe_Saulnier/OneDrive\ -\ INSTITUT\ CURIE/env/STARcounts_hg38_HDBR.txt",
                     sep="\t", head=TRUE, check.names=FALSE)
count0[1:5,1:5]
colnames(dt.ebi2.p)
colnames(dt.supp.p)




sel <- dt.ebi2.p%>%
  filter(!Sample.Characteristic.organism.part. %in% c("brain fragment",
                                                      "diencephalon and midbrain",
                                                      "forebrain and midbrain",
                                                      "forebrain fragment",
                                                      "hindbrain fragment",
                                                      "hindbrain without cerebellum",
                                                      "pituitary and diencephalon"))
sel.supp <- dt.supp.p%>%
  filter(!tissue_clean %in% c("brain fragment",
                              "diencephalon and midbrain",
                              "forebrain and midbrain",
                              "forebrain fragment",
                              "hindbrain fragment",
                              "hindbrain without cerebellum",
                              "pituitary and diencephalon"))

length(unique(sel$Sample.Characteristic.individual.)) # 84


ggplot(as.data.frame(table(sel$Sample.Characteristic.karyotype.)), 
       aes(x = Var1, y = Freq)) +
  geom_col(aes(fill = Var1), position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Freq, group = Var1),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 300),
                     expand = c(0, 0)) + # labels = comma, 
  scale_fill_manual(values = c("#6caed6", "#cbaa7a", "grey")) +
  labs(title="Gender distribution", x="", y="") +
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        axis.line.x= element_line(colour="black", linewidth=0.5),
        axis.line.y= element_line(colour="black", linewidth=0.5),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.text=element_text(size=13, color="black"),
        # axis.text.x=element_text(size=11, color="black", angle=0, hjust=1),
        legend.title=element_text(size=11, color="black", face="bold"),
        legend.text=element_text(size=11, color="black"))



d <- sel %>%
  group_by(Sample.Characteristic.individual.) %>%
  summarise(n_samples = n())
d2 <- sel.supp %>%
  group_by(id_hdbr) %>%
  summarise(n_samples = n())
cmp <- d %>%
  rename(id_hdbr = Sample.Characteristic.individual.) %>%
  full_join(d2, by = "id_hdbr", suffix = c("_d", "_d2")) %>%
  mutate(
    same_n_samples = n_samples_d == n_samples_d2
  )
cmp <- cmp %>% filter(same_n_samples == FALSE)

x <- sel.supp%>%filter(id_hdbr==1923)%>%pull(ID_SUPP)
x2 <- sel%>%filter(Sample.Characteristic.individual.==1923)%>%pull(ID_EBI)%>%unique()
which(x2 %in% x)


sel <- sel%>%
  mutate(multiple_batches = grepl(",", Sample.Characteristic.block.))
sel.multibatch <- sel %>%
  filter(multiple_batches)

dt.ebi2.p <- dt.ebi2.p%>%
  mutate(multiple_batches = grepl(",", Sample.Characteristic.block.))
dim(dt.ebi2.p%>%filter(multiple_batches)) # 34



ggplot(sel %>%
         group_by(Sample.Characteristic.individual.) %>%
         summarise(n_samples = n()), 
       aes(x=n_samples)) + 
  geom_histogram(bins=10, fill="#D2D0A0", color="black") + 
  scale_y_continuous(expand=expansion(mult=c(0.01, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) + 
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        # axis.line.x= element_blank(), axis.line.y= element_blank(),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=13, color="black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=1),
        legend.title=element_text(size=12, color="black", face="bold"),
        legend.text=element_text(size=12, color="black")) + 
  labs(title="Number of samples per individual", x="", y="# samples")












