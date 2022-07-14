library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(plyr)
library(vegan) #2.5-7
library(ape)
library(ggthemes)
library(scales)


set.seed(1851)
permutations = 9999

## get metadata
metadata <- read_delim("data/metadata.txt") %>% 
  as.data.frame()

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt") %>% 
  as.data.frame(row.names = featureid)

table <- orig.table %>% 
  select(contains(metadata$sampleid), Taxon) %>% 
  rownames_to_column(var = "ASV") %>% 
  mutate(ASV = paste0("ASV", ASV)) %>% 
  column_to_rownames(var = "ASV") %>% 
  select(-Taxon)


# Transpose table
table <- table %>% t()


## PCoA plots for all sample types

#quarter root transformation of table
table.transformed <- table^1/4

# Calculate jaccard & bray-curtis distances
bc_dist <- vegdist(table.transformed, method= "bray")
j_dist <- vegdist(table.transformed, method= "jaccard")

# Calculate PCoA with BC/J distances
bc_pcoa <- pcoa(bc_dist, correction = "cailliez") %>% 
  mutate(percent_var = (.$values$Eigenvalues/.$trace)*100)
j_pcoa <- pcoa(j_dist, correction = "cailliez") %>% 
  mutate(percent_var = (.$values$Eigenvalues/.$trace)*100)
## Pull PCo1/PCo2 vectors for plotting

# Rename values
bc_pcoa_vectors <- bc_pcoa$vectors %>% as_tibble(rownames = "sampleid") %>% 
  select(sampleid, Axis.1,Axis.2)
colnames(bc_pcoa_vectors) <- c("sampleid", "bc_PCo1", "bc_PCo2")

j_pcoa_vectors <- j_pcoa$vectors %>% as_tibble(rownames = "sampleid") %>% 
  select(sampleid, Axis.1,Axis.2)
colnames(j_pcoa_vectors) <- c("sampleid", "j_PCo1", "j_PCo2")


bc_hull <- as.data.frame(bc_pcoa_vectors) %>% 
  left_join(metadata) %>% 
  group_by(Group) %>% 
  slice(chull(bc_PCo1, bc_PCo2))

j <- as.data.frame(j_pcoa_vectors) %>% 
  left_join(metadata) %>% 
  group_by(Group) %>% 
  slice(chull(j_PCo1, j_PCo2))

## Pull variance represented for axis lables
bc_variance_rep <- round(bc_pcoa$percent_var[1:2],2)
j_variance_rep <- round(j_pcoa$percent_var[1:2],2)



#Join with metadata for plotting
pcoa.vectors.metadata <- inner_join(metadata, bc_pcoa_vectors, by = "sampleid") %>% 
  inner_join(., j_pcoa_vectors)

## Colors used in following PCoA plot
display.brewer.pal(3, "Set2")
# plot bc pcoa for all samples

show_col(gdocs_pal()(3))



####################################33
# ADONIS
######################################

bc_adonis <- adonis(bc_dist ~ Group , 
                            permutations = permutations, data = metadata)


j_adonis <- adonis(j_dist~Group,
                   permutations = permutations, data = metadata)

permanova_res <- as.numeric()
permanova_res["bc_F"] <- bc_adonis$aov.tab$F.Model[1]
permanova_res["bc_p"] <- bc_adonis$aov.tab$`Pr(>F)`[1]
permanova_res["j_F"] <- j_adonis$aov.tab$F.Model[1]
permanova_res["j_p"] <- j_adonis$aov.tab$`Pr(>F)`[1]




pcoa.vectors.metadata

pcoa.vectors.metadata %>% 
  ggplot(aes(x = bc_PCo1, y = bc_PCo2,
             color = Group)) +
  
  # geom_polygon(data = bc_hull,
  #              aes(fill = Group,
  #                  color = Group),
  #              alpha = 0.1,
  #              show.legend = FALSE) +
  stat_ellipse(show.legend = F) +
  
  # stat_ellipse(show.legend = F) +
  geom_point(aes(size = 3,
                 shape = Group)) +
  
  labs(x = paste("PCo1 - ",bc_variance_rep[1], "%", sep = ""),
       y = paste("PCo2 - ",bc_variance_rep[2], "%", sep = "")) +
  
  scale_color_brewer(palette = "Dark2") +
  
  theme_cowplot() +
  

  
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold",
                                   size = 14),
        legend.position = c(0.05,0.95),
        axis.title = element_text(face = "bold"),
        axis.text = element_text (face = "bold")) +
  
  guides(size = "none",
         color = guide_legend(override.aes = list(size = 3)))


# ggsave("plots/bc_pcoa.png",
#        width = 5,
#        height = 5,
#        dpi = 600,
#        units = c("in"),
#        bg = "white") 


pcoa.vectors.metadata %>% 
  ggplot(aes(x = j_PCo1, y = j_PCo2,
             color = Group)) +
  
  # geom_polygon(data = bc_hull,
  #              aes(fill = Group,
  #                  color = Group),
  #              alpha = 0.1,
  #              show.legend = FALSE) +
  stat_ellipse(show.legend = F) +
  
  # stat_ellipse(show.legend = F) +
  geom_point(aes(size = 3,
                 shape = Group)) +
  
  labs(x = paste("PCo1 - ",j_variance_rep[1], "%", sep = ""),
       y = paste("PCo2 - ",j_variance_rep[2], "%", sep = "")) +
  
  scale_color_brewer(palette = "Dark2",
                     labels = c("Health", "Asthma")) +
  
  scale_shape(labels = c("Health", "Asthma")) +
  
  theme_cowplot() +
  
  
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "bold",
                                   size = 14),
        legend.position = c(0.05,0.95),
        axis.title = element_text(face = "bold"),
        axis.text = element_text (face = "bold")) +
  
  guides(size = "none",
         color = guide_legend(override.aes = list(size = 3)))

# 
# ggsave("plots/j_pcoa.png",
#        width = 7,
#        height = 5,
#        dpi = 600,
#        units = c("in"),
#        bg = "white")  

permanova_res


