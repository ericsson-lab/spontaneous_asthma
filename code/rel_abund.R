library(tidyverse)
library(cowplot)
library(viridis)  # Load viridis for the viridis palette
library(stringr)
library(forcats)
library(ggrepel)
library(pals)



## get metadata
metadata <- read_delim("data/metadata.txt") %>% 
  as.data.frame()

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt") %>% 
  as.data.frame(row.names = featureid)

table <- orig.table %>% 
  select(contains(metadata$sampleid), Taxon) %>% 
  rownames_to_column(var = "ASV") %>% 
  mutate(ASV = paste0("ASV", ASV))

taxonomy <- table %>% dplyr::select(ASV, Taxon) %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  as.data.frame()


table.tax <- table %>% dplyr::select(-c(Taxon)) %>% 
  pivot_longer(-ASV) %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
                        "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon") %>% 
  as.data.frame()

table.family <-table.tax %>% 
  filter(level == "Family") %>%
  group_by(name, taxon) %>%
  summarise(count = sum(value)) %>%
  ungroup() %>%
  mutate(taxon = gsub("D_4__", "", taxon)) %>% 
  rename(sampleid = name) %>% 
  ungroup() %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) %>% 
  as.data.frame()

other.list <- table.family %>% 
  ungroup() %>%
  group_by(taxon) %>% 
  dplyr::summarise(pool = rel_abund < 0.05, .groups = "drop")
other.list %>% count(pool)

abund.list <-table.family %>% 
  group_by(taxon) %>% 
  summarise(abund.40 = rel_abund >= 0.4, .groups = "drop") 

table.family
other.list

table <- inner_join(table.family, other.list, by = "taxon")  %>% 
  # inner_join(., abund.list, by = "taxon") %>% 
  # mutate(abund.40 = case_when(abund.40 == TRUE ~ 1,
  #                             abund.40 ==FALSE ~ 0.6,
  #                             TRUE ~ 0.6)) %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  left_join(., metadata, by = "sampleid")
table



# Gather all of the family to be plotted and group by phylum
#https://stackoverflow.com/questions/62627480/how-to-creat-a-bar-graph-of-microbiota-data-with-one-color-for-higher-taxonomic
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
                   "Porphyromonadaceae",
                   "Prevotellaceae",
                   "Weeksellaceae")

firmicutes <- c("Bacillaceae",
                "Lachnospiraceae",
                "Ruminococcaceae",
                "Staphylococcaceae",
                "Streptococcaceae")

fusobacteria <- c("Fusobacteriaceae")

proteobacteria <- c("Burkholderiaceae",
                    "Caulobacteraceae",
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

c('Actinobacteria','Bacteroidetes',
  'Firmicutes',"Fusobacteria", 
  "Proteobacteria","Tenericutes")
reorder 
table.family %>%  
  mutate(phylum = case_when(taxon %in% proteobacteria ~ "Proteobacteria",
                            taxon %in% bacteroidetes ~ "Bacteroidetes",
                            taxon %in% fusobacteria ~ "Fusobacteria",
                            taxon %in% actinobacteria ~ "Actinobacteria",
                            taxon %in% firmicutes ~ "Firmicutes",
                            taxon %in% tenericutes ~ "Tenericutes",
                            TRUE ~ "Other"),
         taxon = case_when(taxon %in% all_family ~ taxon,
                            TRUE ~ "Other")) %>% 
  mutate(phylum=factor(phylum, levels=c(phylums, "Other")),
                taxon=fct_reorder(taxon, 10*as.integer(phylum) + grepl("Other", taxon))) %>% 

  left_join(., metadata, by = "sampleid") %>% 
  group_by(sampleid) %>% 
  mutate(rel_abund = count/sum(count)) 


table %>% 
  ggplot(aes(x = sampleid,
             y = rel_abund,
             fill = taxon)) +
  facet_grid(~Group, scales = "free_x",
             space = "free_x",
             switch = "both"
             ) +
  labs(y = "Family Relative Abundance") +
  geom_col() +
  
  scale_fill_manual(values = as.vector(c("#809580",alphabet2(), "#808080"))) +
  
  
  theme_cowplot() +

  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent) +

  guides(fill=guide_legend(ncol=1,
                           override.aes = list(size = 0.2))) +
  
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

ggsave("plots/family_rel_abund.png",
       bg = "white",
       height = 7,
       width = 7,
       units = c("in"),
       dpi = 600)

as.array(c("#F8A19G",alphabet2(), "#808080"))

#Chitino
#Moraxe
# Myco
# psueuo

# %>% 
#   pivot_wider(names_from = "name",
#               values_from = "count")






# ggsave("plots/alpha_diversities.png",
#        width = 8,
#        height = 6,
#        units = c("in"),
#        bg = "white",
#        dpi = 600)

getwd()
write_csv(alpha_metrics, "data/alpha_diversity.csv")
pivot_longer(-sampleid, names_to = "alpha_metric", values_to = "value") %>% 
  left_join(metadata, .)






featur


