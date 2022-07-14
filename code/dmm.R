library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyverse)


set.seed(1851)
permutations = 9999

## get metadata
metadata <- read_delim("data/metadata.txt") 

## get feature table and assign taxonomy as rownames
orig.table <- read_delim("data/feature.table.txt") %>% 
  as.data.frame(row.names = featureid)


table <- orig.table %>% 
  select(contains(metadata$sampleid), Taxon)
head(table)

taxonomy <- table %>% dplyr::select(Taxon) %>% 
  separate(col = "Taxon", 
           into = c("Kingdom", "Phylum", "Class", 
                    "Order", "Family", "Genus"), 
           sep = ";",
           fill = "right") %>% 
  rownames_to_column(var = "ASV") %>% 
  mutate(ASV = paste0("ASV", ASV))
  

write_csv(taxonomy, "data/taxonomy.csv")
table

table.tax <- table %>% dplyr::select(-Taxon) %>% 
  rownames_to_column(var = "ASV") %>% 
  pivot_longer(-ASV) %>% 
  inner_join(., taxonomy, by = "ASV") %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", 
                        "Order", "Family", "Genus"),
               names_to = "level",
               values_to = "taxon")

table.tax
table.asv <- table %>% 
  rownames_to_column(var = "ASV") %>% 
  select(-Taxon) %>% 
  mutate(ASV = paste0("ASV", ASV)) %>% 
  column_to_rownames(var = "ASV")
table.asv <- table.asv[rowSums(table.asv[])>0,] %>% t() %>% as.matrix()


# family.table <- table.tax %>%
#   filter(level == "Family") %>%
#   group_by(name, taxon) %>%
#   summarise(family_count = sum(value)) %>%
#   ungroup() %>%
#   mutate(taxon = gsub("D_4__", "", taxon))%>%
#   pivot_wider(names_from = "taxon",
#               values_from = "family_count") %>%
#   column_to_rownames(var = "name") %>%
#   as.matrix()


row_number(metadata)
table.asv

fit <- lapply(1:10, dmn, count = table.asv, verbose=TRUE)
# fit.phylum <- lapply(1:10, dmn, count = phylum.table, verbose=TRUE)

lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

# lplc <- base::sapply(fit.phylum, DirichletMultinomial::laplace) # AIC / BIC / Laplace
# aic  <- base::sapply(fit.phylum, DirichletMultinomial::AIC) # AIC / BIC / Laplace
# bic  <- base::sapply(fit.phylum, DirichletMultinomial::BIC) # AIC / BIC / Laplace
# 


best <- fit[[which.min(unlist(lplc))]]
mixturewt(best)


assignments.asv <- apply(mixture(best), 1, which.max) %>% as.data.frame()
colnames(assignments.family) <- c("dmm_group")

colnames(assignments.asv) <- c("Assignment")

best

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  
  ggsave(paste0("plots/com",k,"_asv.png"),
         height = 20,
         width = 5,
         units = c("in"),
         bg = "white",
         dpi = 600)
}



vignette("DirichletMultinomial")


