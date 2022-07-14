
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
table.tax

table.family <- table.tax %>% 
    filter(level == "Family") %>% 
    group_by(name, taxon) %>% 
    summarise(count = sum(value)) %>% 
    ungroup() %>%
    mutate(taxon = gsub("D_4__", "", taxon)) %>% 
  pivot_wider(names_from = "name",
              values_from = "count")

table.family[is.na(table.family)]<-"Undefined" 

table.family <- table.family %>% 
  column_to_rownames(var = "taxon")
table.family[rowSums(table.family[])>0,]


##Go through each row and determine if a value is zero

#######################################
  
## psuedocount 1 to every cell in table
table.family <- table.family
table.family.pseudo <- table.family + 1
## Mutation to calculate log2 fold change
log2FC <- table.family.pseudo %>% 
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
# table.family %>% 
#   rownames_to_column(var = "Family") %>% 
#   write_csv("data/metabo_genus.csv")
# metadata$Group
conds <- c(rep("SpA",26), rep("Healthy",11))

x.family <- aldex.clr(
  reads = table.family,
  conds = conds, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)

x_tt.family <- aldex.ttest(
  x.family, 
  paired.test = FALSE, 
  verbose = FALSE)

aldex_out.family <- data.frame(x_tt.family) %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(., log2FC)

aldex_out.family %>% 
  filter(log2FC > 1)


## Plots generated using Benjami-Hochberg corrected Welch's T test
## and 20 fold change in mean family counts between groups
family <- EnhancedVolcano(aldex_out.family,
                               lab = aldex_out.family$ASV,
                               x = 'log2FC',
                               y = 'wi.ep',
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               drawConnectors = TRUE,
                               typeConnectors = "open",
                               arrowheads = FALSE,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               ylim = c(0,10),
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
          legend.text.align	= 0 ,
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  )

family


aldex_out.genus 
  

# ggsave("plots/aldex2/family_aldex2.png",
#        dpi = 600,
#        width = 10,
#        height = 7,
#        units = c("in"),
#        bg = "white")


save <- aldex_out.family %>% 
  filter(-log10(wi.ep) > -log10(0.05))

# write.csv(save, "data/sig_families.csv")

# fam.drop <- c("Burkholderiaceae", "Actinomycetaceae", "Sphingobacteriaceae",
#               "Fusobacteriaceae", "Clostridiaceae 1", "Weeksellaceae")
sav
sig_family <- save %>% 
  rownames_to_column(var = "family") %>% 
  arrange(log2FC) %>% 
  select(ASV) %>% 
  # filter(!ASV %in% fam.drop) %>% 
  as.list()



plot.data <- table.family %>% 
  rownames_to_column(var = "family") %>% 
  pivot_longer(-family,
               names_to = "sampleid",
               values_to = "count") %>% 
  filter(family %in% sig_asv$ASV) %>% 
  
  left_join(., metadata) %>% 
  select(family, count, Group) %>% 
  group_by(family)
plot.data

plot.data$family <- factor(plot.data$family, levels = sig_family$ASV)

plot.data %>% 
  mutate(rel_abund = count / sum(count)) %>% 
  ggplot(aes(x = Group, y = rel_abund,
             fill = Group)) +
  
  geom_boxplot(aes(fill = NULL),
               show.legend = F) +
  
  geom_dotplot(binaxis = 'y',
               stackdir = "center",
               dotsize = 1.25,
               show.legend = F) +
  facet_wrap(~family, scales = "free") +
  
  facet_wrap(~factor(family, levels = c( "Pseudomonadaceae",     "Xanthobacteraceae",    "Sphingomonadaceae",   
                                         "Burkholderiaceae" ,    "Propionibacteriaceae", "Actinomycetaceae" ,   
                                         "Sphingobacteriaceae",  "Staphylococcaceae",    "Corynebacteriaceae",  
                                         "Fusobacteriaceae",     "Clostridiaceae 1",     "Weeksellaceae"   ,    
                                         "Mycoplasmataceae" , "Xanthomonadaceae", "Streptococcaceae", "Undefined")), 
             scales = "free",
             ncol = 4) +
  
  labs(y = "Family Relative Abundance") +
  
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = c("Health", "Asthma")) +
  
  theme_cowplot() +
  
  theme(axis.title.x = element_blank(),
        axis.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        
        strip.background = element_rect(color = NA,
                                    fill = NA),
        strip.text = element_text(size = 12,
                                  face = "bold")
  )

ggsave("plots/sig_families.png",
       width = 10,
       height = 9,
       dpi = 600,
       bg = "white",
       units = c("in"))


# plot.data.rel.abund <- table.family %>% 
#   rownames_to_column(var = "family") %>% 
#   pivot_longer(-family,
#                names_to = "sampleid",
#                values_to = "count") %>% 
#   group_by(sampleid) %>% 
#   mutate(rel_abund = count / sum(count)) %>% 
#   filter(family %in% sig_asv$ASV) %>% 
#   
#   left_join(., metadata) %>% 
#   group_by(family)
# 
# 
# plot.data.rel.abund %>% 
#   ggplot(aes(x = Group, y = rel_abund,
#              fill = Group)) +
#   
#   geom_boxplot(aes(fill = NULL),
#                show.legend = F) +
#   
#   geom_dotplot(binaxis = 'y',
#                stackdir = "center",
#                dotsize = 2,
#                show.legend = F) +
#   
#   facet_wrap(~factor(family, levels = c( "Pseudomonadaceae",     "Xanthobacteraceae",    "Sphingomonadaceae",   
#                                          "Burkholderiaceae" ,    "Propionibacteriaceae", "Actinomycetaceae" ,   
#                                          "Sphingobacteriaceae",  "Staphylococcaceae",    "Corynebacteriaceae",  
#                                          "Fusobacteriaceae",     "Clostridiaceae 1",     "Weeksellaceae"   ,    
#                                          "Mycoplasmataceae" , "Xanthomonadaceae", "Streptococcaceae", "Undefined")), 
#              scales = "free",
#              ncol = 4) +
#   
#   labs(y = "Family Count") +
#   
#   scale_fill_brewer(palette = "Dark2") +
#   
#   scale_y_continuous(labels = scales::percent) +
#   
#   theme_cowplot() +
#   
#   theme(axis.title.x = element_blank(),
#         axis.text = element_text(face = "bold"),
#         axis.title.y = element_text(face = "bold"),
#         
#         strip.background = element_rect(color = NA,
#                                         fill = NA),
#         strip.text = element_text(size = 12,
#                                   face = "bold")
#   )
# 
# ggsave("plots/sig_families_rel_abund.png",
#        width = 8.5,
#        height = 11,
#        dpi = 600,
#        bg = "white",
#        units = c("in"))
# 
# 
# write_csv(aldex_out.family, "data/aldex2_families.csv")
# 





