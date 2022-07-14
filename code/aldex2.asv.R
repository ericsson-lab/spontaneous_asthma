
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
table %>% select(-Taxon) %>% 
  write_csv("data/metabo_upload.csv")


taxonomy <- table %>% dplyr::select(ASV, Taxon) %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  as.data.frame()

# table.tax <- table %>% dplyr::select(-c(Taxon)) %>% 
#   pivot_longer(-ASV) %>% 
#   inner_join(., taxonomy, by = "ASV") %>% 
#   pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
#                         "Order", "Family", "Genus"),
#                names_to = "level",
#                values_to = "taxon") %>% 
#   as.data.frame()
# 
# table.genus <- table.tax %>% 
#   filter(level == "Genus") %>%
#   group_by(name, taxon) %>%
#   summarise(count = sum(value)) %>%
#   ungroup() %>%
#   mutate(taxon = gsub("D_5__", "", taxon)) %>% 
#   pivot_wider(names_from = "name",
#               values_from = "count")
# 
# table.genus[is.na(table.genus)]<-"Undefined" 
# 
# table.genus <- table.genus %>% 
#   column_to_rownames(var = "taxon")

table.asv <- table %>% 
  select(-Taxon) %>% 
  column_to_rownames(var = "ASV")
#######################################

## psuedocount 1 to every cell in table
table.asv <- table.asv +1

## Mutation to calculate log2 fold change
log2FC <- table.asv %>% 
  # move family col back to rowname
  rownames_to_column(var = "ASV") %>% 
  # lengthen table to perform mutations
  pivot_longer(-ASV, names_to = "sampleid") %>% 
  left_join(., metadata, by = "sampleid") %>% 
  # Group by family and time point to find mean within each family at
  # both time points
  group_by(ASV, Group) %>% 
  summarise(mean = mean(value), .groups = "drop") %>% 
  ## widen table to split up time points
  pivot_wider(names_from = Group, values_from = mean)  %>% 
  # group by family for future plotting and perform log2 FC calculation
  # + FC = up in Chronic (right); - FC = up in Healthy (left)
  group_by(ASV) %>% 
  summarise(log2FC = log2(`SpA`/`Healthy`)) %>% 
  arrange(desc(log2FC))

log2FC
## DESEq2 Analysis, following 
# https://microbiome.github.io/OMA/differential-abundance.html#aldex2
# Build condition table for DESEq2

table.family
metadata$Group
conds <- c( rep("Healthy",11), rep("SpA",26))

x.asv <- aldex.clr(
  reads = table.asv,
  conds = conds, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)

x_tt.asv <- aldex.ttest(
  x.asv, 
  paired.test = FALSE, 
  verbose = FALSE)


aldex_out.asv <- data.frame(x_tt.asv) %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(., log2FC)
aldex_out.asv

## Plots generated using Benjami-Hochberg corrected Welch's T test
## and 20 fold change in mean family counts between groups
asv.plot <- EnhancedVolcano(aldex_out.asv,
                              lab = aldex_out.asv$ASV,
                              x = 'log2FC',
                              y = 'we.ep',
                              pCutoff = 0.05,
                              FCcutoff = 1,
                              drawConnectors = TRUE,
                              title = NULL,
                              subtitle = NULL,
                              caption = NULL,
                              axisLabSize = 12,
                              legendPosition = "bottom",
                              legendLabels = c("NS", 
                                               "Log2 FC", 
                                               "p < 0.05",
                                               "p < 0.05 and Log2 FC"),
                              legendLabSize = 10,
                              labSize = 4,
                              labFace = "bold",
                              colAlpha = 0.9) +
  theme_set(
    theme(axis.text = element_text(color = "black",
                                   face = "bold"),
          axis.title = element_text(color = "black",
                                    face = "bold"),
          legend.text = element_text(color = "black",
                                     face = "bold"),
          legend.text.align	= 0)
  )
asv.plot

aldex_out.asv %>% 
  filter()




