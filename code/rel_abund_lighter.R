# library(ggtext)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tidyr)
library(stringr)
library(forcats)
library(ggrepel)



library(viridis)
library(pals)
library(ggthemes)


## get metadata
metadata <- read_delim("data/metadata.txt") 
  # mutate(sampleid = str_replace_all(sampleid, "Fe-BAL-", "FeBAL")) %>% 
  # separate(sampleid, sep = "-", into = c("sampleid")) 

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt", col_types = NULL) %>% 
  as.data.frame()
head(orig.table)
# colnames(orig.table) = gsub("Fe-BAL-", "FeBAL", colnames(orig.table))
# colnames(orig.table) = gsub("-.*", "", colnames(orig.table))

# Add ASV column to parse taxonomy out
table <- orig.table %>% 
  rownames_to_column("ASV") %>% 
  mutate(ASV = paste0("ASV", ASV))

# Generate taxonomy and separate identifiers
taxonomy <-table %>% select(ASV, Taxon) %>% 
  separate(Taxon,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";",
           fill = "right")

setdiff(metadata$sampleid, colnames(orig.table))

# Make feature table longer
table.counts <- table %>% 
  select(ASV, contains(metadata$sampleid)) %>% 
  pivot_longer(-c(ASV),
               names_to = "sampleid",
               values_to = "count")

# Join metadata to feature table and taxonomy
table.rel.abund <- inner_join(metadata, 
                              table.counts, 
                              by = "sampleid") %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  group_by(sampleid) %>% 
  # calculate % relative abundance
  mutate(rel_abund = count / sum(count)) %>% 
  ungroup() %>% 
  # Break up taxonomy
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  # Select only family level
  filter(level == "Family",) %>%
  # Remove string header from silva/gg(?)
  mutate(taxon = gsub("D_4__", "", taxon)) %>% 
  ungroup() %>% 
  group_by(sampleid, taxon) %>% 
  # sum up mean rel abund for each taxa within a sample
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>% 
  inner_join(., metadata, by = "sampleid")

# Generate condition for Other pool.
# Logical condition evaluates if the max mean relabund in any sample
# in any sample is < 1. If yes, taxa is grouped into other.
taxon.pool <- table.rel.abund %>% 
  group_by(taxon) %>% 
  mutate(pool = max(rel_abund) < 0.05, .groups = "drop") %>% 
  mutate(abund.40 = max(rel_abund) >= 0.40, .groups = "drop") %>% 
  mutate(alpha = case_when(abund.40 == TRUE ~ 1,
                           abund.40 == FALSE ~ 0.5))


# Count # of max(mean(rel_abund)) taxa > 1% to set color pallete in plot
taxon.pool %>% count(pool)

# Join pool table to rel.abund and rename family to other is mean(max(rel_abund)) < 1%
table 
inner_join(table.rel.abund, taxon.pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  inner_join(., abund.40, by = "taxon") 

phylums <- c('Actinobacteria','Bacteroidetes',
             'Firmicutes',"Fusobacteria", 
             "Proteobacteria","Tenericutes")

actinobacteria <- c("Corynebacteriaceae",
                    "Microbacteriaceae",
                    "Propionibacteriaceae")

bacteroidetes <- c("Bacteroidaceae",
                   "Chitinophagaceae",
                   "Flavobacteriaceae",
                   "Muribaculaceae",
                   "Porphyromonadaceae")

firmicutes <- c("Bacillaceae",
                "Lachnospiraceae",
                "Ruminococcaceae",
                "Staphylococcaceae",
                "Streptococcaceae")

fusobacteria <- c("Fusobacteriaceae")

proteobacteria <- c("Caulobacteraceae",
                    "Moraxellaceae",
                    "Neisseriaceae",
                    "Nitrosomonadaceae",
                    "Pasteurellaceae",
                    "Pseudomonadaceae",
                    "Rhizobiaceae",
                    "Sphingomonadaceae",
                    "Xanthomonadaceae")

tenericutes <- ("Mycoplasmataceae")

all_family <- c(actinobacteria,
                bacteroidetes,
                firmicutes,
                fusobacteria,
                proteobacteria,
                tenericutes)


reorder 
taxon.pool.ordered <-taxon.pool %>%  ungroup() %>% 
  mutate(phylum = case_when(taxon %in% proteobacteria ~ "Proteobacteria",
                            taxon %in% bacteroidetes ~ "Bacteroidetes",
                            taxon %in% fusobacteria ~ "Fusobacteria",
                            taxon %in% actinobacteria ~ "Actinobacteria",
                            taxon %in% firmicutes ~ "Firmicutes",
                            taxon %in% tenericutes ~ "Tenericutes",
                            TRUE ~ "Other"),
         taxon = if_else(pool, "Other", taxon)) %>% 
  mutate(phylum=factor(phylum, levels=c(phylums, "Other")),
         taxon=fct_reorder(taxon, 10*as.integer(phylum) + grepl("Other", taxon)))

# Full gdocs_pal() pallete. 
pal <- c("#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6",
         "#dd4477", "#66aa00", "#b82e2e", "#316395", "#994499", "#22AA99",
         "#AAAA11", "#6633CC", "#E67300", "#8B0707", "#651067", "#329262", 
         "#5574A6", "#3B3EAC")
palette = colorRampPalette(pal)
#order <- sample(palette(24))
#order
# [1] "#3E8779" "#0692C4" "#39683F" "#506AA7" "#3B3EAC" "#AC5688" "#81409F"
# [8] "#E26E00" "#897413" "#8F4698" "#943B48" "#4A5364" "#94822E" "#418F99"
# [15] "#8F7B5A" "#F27706" "#85139E" "#3A6095" "#700D49" "#6EAA4C" "#8C970B"
# [22] "#BE4032" "#961506" "#3366CC"

taxon.pool.ordered %>% 
  ggplot(aes(x = sampleid,
             y = rel_abund,
             fill = taxon,
            alpha = as.factor(alpha))) +
  facet_grid(~Group, scales = "free_x",
             space = "free_x",
             switch = "both"
  ) +
  
  scale_alpha_discrete(range = c(0.4, 1)) +
  labs(y = "Family Relative Abundance") +
  geom_col()  +
  theme_cowplot() +
  
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent) +
  
  # scale_fill_viridis(discrete = T) +
   scale_fill_manual(values = as.vector(c(order, "#808080"))) +
  # scale_fill_manual(values = as.vector(c(gdocs_pal()))) +
  
  
  
  guides(fill=guide_legend(ncol=1,
                           override.aes = list(size = 0.2)),
         alpha = "none"
         )+
  
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        
        strip.background = element_rect(color = NA,
                                        fill = NA),
        strip.placement = "outside",
        strip.text = element_text(face = "bold",
                                  size = 14),
        
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",
                                   size = 10))


#order <- sample(palette(24))
#order
# [1] "#3E8779" "#0692C4" "#39683F" "#506AA7" "#3B3EAC" "#AC5688" "#81409F"
# [8] "#E26E00" "#897413" "#8F4698" "#943B48" "#4A5364" "#94822E" "#418F99"
# [15] "#8F7B5A" "#F27706" "#85139E" "#3A6095" "#700D49" "#6EAA4C" "#8C970B"
# [22] "#BE4032" "#961506" "#3366CC"
getwd()
ggsave("plots/rel_abund_5pAbund_40p_highlited_lighter.png",
       dpi = 600,
       width = 8,
       height = 7,
       units = c("in"),
       bg = "white")


taxon.pool.ordered %>% 
  filter(rel_abund >= 0.4) %>% 
  count(taxon)
arrange(rel_abund)
  
  
