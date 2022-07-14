library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)
library(dplyr)


## get metadata
metadata <- read_delim("data/metadata.txt") %>% 
  as.data.frame()

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt") %>% 
  as.data.frame(row.names = featureid)

table <- orig.table %>% 
  dplyr::select(contains(metadata$sampleid), Taxon) %>% 
  rownames_to_column(var = "ASV") %>% 
  mutate(ASV = paste0("ASV", ASV))
  


table
head(table)

taxonomy <- table %>% dplyr::select(ASV, Taxon) %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  as.data.frame()


alpha_metrics <- table %>% dplyr::select(-Taxon) %>% 
  pivot_longer(-ASV,
               names_to = "sampleid",
               values_to = "counts") %>% ungroup() %>% 
 # select(-ASV) %>% 
  group_by(sampleid) %>%  
  summarize(obs_richness = specnumber(counts),
            shannon_e = diversity(counts, 
                                  index = "shannon",
                                  base = exp(1)),
            simpson = diversity(counts, 
                                index = "simpson")) %>% 
  left_join(metadata, ., by = "sampleid")


alpha_metrics %>% 
  dplyr::select(Group, read_count, obs_richness, 
                         shannon_e, simpson) %>% 
  
  pivot_longer(-Group, names_to = "alpha_diversity",
               values_to = "value") %>% 
  group_by(Group,alpha_diversity) %>% 
  summarize(min = min(value),
            q1 = quantile(value, 0.25),
            median = median(value),
            mean = mean(value),
            q3 = quantile(value, 0.75),
            max = max(value),
            sd = sd(value)) %>% 
  arrange(alpha_diversity) %>% 
  write_csv("data/alpha_diversity_sum_stats.csv")



alpha_metrics %>% 
  select(Group, read_count, 
         obs_richness, shannon_e, 
         simpson) %>% 
  pivot_longer(-Group, names_to = "alpha_diversity",
               values_to = "value") %>% 
  mutate(alpha_diversity = case_when(alpha_diversity == "obs_richness" ~ "Observed ASVs",
                                     alpha_diversity == "shannon_e" ~ "Shannon Index",
                                     alpha_diversity == "simpson" ~ "Simpson Index",
                                     alpha_diversity == "read_count" ~ "Feature Count")) %>% 
  ggplot(aes(x= Group, y = value)) +
  
  geom_boxplot() +
  
  geom_dotplot(aes(fill = Group),
               binaxis = "y",
               stackdir = "center",
               show.legend = F) +
  
  facet_wrap(~alpha_diversity,
             scales = "free",
             strip.position = "left") +
  
  scale_fill_brewer(palette = "Dark2") +
  
  theme_cowplot() +
  
  stat_compare_means(method = "wilcox",
                     label.x = 0.75) +
  
  scale_x_discrete(labels = c("Health", "Asthma")) +
  
  theme(strip.background = element_rect(color = NA,
                                        fill = NA),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        
        axis.text = element_text(face = "bold"),
        
        axis.title = element_blank(),
        text = element_text(face =))

ggsave("plots/alpha_diversities.png",
        width = 8,
        height = 6,
        units = c("in"),
        bg = "white",
        dpi = 600)
  
getwd()
write_csv(alpha_metrics, "data/alpha_diversity.csv")
  pivot_longer(-sampleid, names_to = "alpha_metric", values_to = "value") %>% 
  left_join(metadata, .)






featur


