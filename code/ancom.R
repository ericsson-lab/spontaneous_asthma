library(ALDEx2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtext)
library(cowplot)
library(ggrepel)
library(EnhancedVolcano)
library(ggthemes)
library(scales)

# get metadata
metadata <- read_delim("data/metadata.txt") %>% 
  as.data.frame()

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt") %>% 
  as.data.frame(row.names = featureid)

table <- orig.table %>% 
  select(contains(metadata$sampleid), Taxon) %>% 
  rownames_to_column(var = "ASV") %>% 
  mutate(ASV = paste0("ASV", ASV))
table


table
head(table)

taxonomy <- table %>% dplyr::select(ASV, Taxon) %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  as.data.frame()

table.tax <- table %>% dplyr::select(-Taxon) %>% 
  pivot_longer(-ASV) %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
                        "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  as.data.frame()

family.table <- table.tax %>% 
  filter(level == "Family") %>% 
  group_by(name, taxon) %>% 
  summarise(family_count = sum(value)) %>% 
  ungroup() %>% 
  mutate(taxon = gsub("D_4__", "", taxon))

phylum.level <- table.tax %>% 
  filter(level == "Phylum") %>% 
  group_by(name, taxon) %>% 
  summarise(family_count = sum(value)) %>% 
  ungroup() %>% 
  mutate(taxon = gsub("D_1__", "", taxon))

##### ASV

# Removes any NA row names and sets to "Undefined Family"
# Move rowname to col 1
## psuedocount 1 to every cell in table
table.pseudo.asv<- table %>% column_to_rownames(var = "ASV") %>% 
  select(-Taxon)
table.pseudo.asv <- table.pseudo.asv +1
table

table
## Folowing ANCOM Standard Analysis
# https://github.com/FrederickHuangLin/ANCOM#standard-analysis 
otu_data.asv = table.pseudo.asv

meta_data = metadata
meta_data= meta_data %>% rename(Sample.ID = sampleid)

source("code/ANCOM-master/scripts/ancom_v2.1.R")

# Step 1: Data preprocessing

feature_table = otu_data.asv; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Group"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
taxonomy

ancom_sig_ASV_res <- res$out %>% 
  rename(ASV = taxa_id) %>% 
  inner_join(taxonomy, ., by = "ASV") %>% 
  filter(detected_0.9 == TRUE)

write.csv(ancom_sig_ASV_res, "data/ancom_sig_ASV_res.csv")

inner_join(taxonomy, res$out, by = ASV)
res$fig

res$out <- res$out %>% rename(family = taxa_id)


log2FC <- table.pseudo.asv %>% 
  # move family col back to rowname
  rownames_to_column(var = "ASV") %>% 
  # lengthen table to perform mutations
  pivot_longer(-ASV, names_to = "sampleid") %>% 
  left_join(., metadata, by = "sampleid") %>% 
  # Group by family and time point to find mean within each family at
  # both time points
  group_by(ASV, Group) %>% 
  summarise(mean_fam = mean(value), .groups = "drop") %>% 
  ## widen table to split up time points
  pivot_wider(names_from = Group, values_from = mean_fam)  %>% 
  # group by family for future plotting and perform log2 FC calculation
  # + FC = up in Chronic (right); - FC = up in Healthy (left)
  group_by(ASV) %>% 
  summarise(log2FC = log2(`SpA`/`Healthy`)) %>% 
  arrange(desc(log2FC))

log2FC

ancom_res_asv <- read_csv("data/ancom_sig_ASV_res.csv", ) %>% 
  left_join(., log2FC, by = "ASV")
write_csv(ancom_res_asv, "data/ancom_ASV_log2FC_h_v_spa.csv")

### Family

# Removes any NA row names and sets to "Undefined Family"
# Move rowname to col 1
## psuedocount 1 to every cell in table
table.pseudo.family <- family.table %>% pivot_wider(names_from = name,
                             values_from = family_count)  
table.pseudo.family[is.na(table.pseudo.family)]<-"Undefined Family"
table.pseudo.family
table.pseudo.family <- table.pseudo.family %>% column_to_rownames(var = "taxon")
table.pseudo.family <- table.pseudo.family + 1

table
## Folowing ANCOM Standard Analysis
# https://github.com/FrederickHuangLin/ANCOM#standard-analysis 
otu_data.family = table.pseudo.family

meta_data = metadata
meta_data= meta_data %>% rename(Sample.ID = sampleid)

# Step 1: Data preprocessing

feature_table = otu_data.family; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Group"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL

res.family = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
re

res.family$fig
write.csv(res.family$out, "data/ancom_sig_family_res.csv")

inner_join(taxonomy, res$out, by = ASV)
res$fig

res$out <- res$out %>% rename(family = taxa_id)








write_csv(res$out, "data/ANCOM/ancom_spa.csv")


