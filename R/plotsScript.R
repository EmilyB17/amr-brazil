## Making Plots for USDA Grant - updated 3/20/20 with relative abundance data

require(dplyr)
require(tidyr)
require(stringr)
require(ggplot2)

# set ggplot theme
theme_set(theme_bw())

# read in data
counts <- read.table("https://github.com/EmilyB17/amr-brazil/blob/master/data/parsedRelativeAbundance.txt?raw=true", sep = "\t", header = TRUE)

# get list of controls
controls <- paste("G527_94_NoTemplate-DNAextraction2_S94",
                  "G527_96_NoTemplate-LibraryPrep_S96",
                  "G527_95_MockCommunity_S95",
                  "G527_93_NoTemplate-DNAextraction1_S93",
                  "G527_49_2019_DNA_S49", sep = "|")

# remove controls
nocont <- counts %>% 
  mutate(name = as.character(name)) %>% 
  filter(!str_detect(name, controls))

# summarize data for each gene and make wider by farm
sumfarm <- nocont %>% 
  select(-nhits) %>% 
  group_by(class, type, protein, gene, farm, site) %>% 
  summarize(count = sum(count))  %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = count, names_prefix = "farm_", values_fill = list(count = 5))

# add nicer column names
sumv <- nocont %>% 
  
  # fix names of everything for plots
  mutate(`Resistance Type` = gsub("_", " ", protein),
         `Resistance Class` = gsub("_", " ", type),
         `General Type` = gsub("_", " ", class)) %>% 
  rename(Gene = gene) 

## ---- Plots ----

## BARPLOT
# multi drug resistance for farms, body sites
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = farm, y = relabun, color = `Resistance Type`, fill = `Resistance Type`
           )) +
  geom_col() +
  labs(x = "Farm", y = "Relative Abundance", title = "Multi-drug resistance") +
  coord_flip()
# save
ggsave("./data/plots/relabun-barplot-farm.png", plot = last_plot(),
       width = 6, height = 5, units = "in")
# body site
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = site, y = relabun, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_col() +
  labs(x = "Body Site", y = "Relative Abundance", title = "Multi-drug resistance") +
  coord_flip()
# save
ggsave("./data/plots/relabun-barplot-body-site.png", plot = last_plot(),
       width = 6, height = 5, units = "in")

## BOXPLOT
# multi drug resistance for farms, body sites
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = farm, y = relabun, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_boxplot() +
  labs(x = "Farm", y = "Relative Abundance", title = "Multi-drug resistance")

# save
ggsave("./data/plots/relabun-boxplot-farm.png", plot = last_plot(),
       width = 6, height = 5, units = "in")

# body site
ggplot(data = filter(sumv, type == "Multi-drug_resistance"), 
       aes(x = site, y = relabun, color = `Resistance Type`, fill = `Resistance Type`)) +
  geom_boxplot() +
  labs(x = "Body Site", y = "Relative Abundance", title = "Multi-drug resistance")

# save
ggsave("./data/plots/relabun-boxplot-body-site.png", plot = last_plot(),
       width = 6, height = 5, units = "in")
