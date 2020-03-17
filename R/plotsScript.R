## Making Plots for USDA Grant

require(dplyr)
require(tidyr)
require(stringr)
require(ggplot2)

# set ggplot theme
theme_set(theme_bw())

# read in data
counts <- read.table("https://github.com/EmilyB17/amr-brazil/blob/master/data/parsedGeneCounts.txt?raw=true", sep = "\t", header = TRUE)

# summarize data for each gene and make wider by farm
sumfarm <- counts %>% 
  select(-nhits) %>% 
  group_by(class, type, protein, gene, farm, site) %>% 
  summarize(count = sum(count))  %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = count, names_prefix = "farm_", values_fill = list(count = 0))

# make vertical and add nicer column names
sumv <- sumfarm %>% 
  pivot_longer(cols = starts_with("farm_"), names_to = "farm", names_prefix = "farm_", values_to = "count") %>% 
  mutate(logcount = log1p(count)) %>% 
  # fix names of everything for plots
  mutate(`Resistance Type` = gsub("_", " ", protein),
         `Resistance Class` = gsub("_", " ", type),
         `General Type` = gsub("_", " ", class)) %>% 
  rename(Gene = gene) %>% 
  # add log transformation of gene counts
  mutate(logcount = log1p(count))

## ---- Plots ----

## BARPLOT
# multi drug resistance for farms, body sites
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = farm, y = logcount, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_col() +
  labs(x = "Farm", y = "log(Number of genes)", title = "Multi-drug resistance") +
  coord_flip()
# body site
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = site, y = logcount, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_col() +
  labs(x = "Body Site", y = "log(Number of genes)", title = "Multi-drug resistance") +
  coord_flip()

## BOXPLOT
# multi drug resistance for farms, body sites
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = farm, y = logcount, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_boxplot() +
  labs(x = "Farm", y = "log(Number of genes)", title = "Multi-drug resistance")
# body site
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = site, y = logcount, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_boxplot() +
  labs(x = "Body Site", y = "log(Number of genes)", title = "Multi-drug resistance")
