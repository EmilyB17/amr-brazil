## Summary statistics for AMR data

# goals:
# presence/absence surveillance of genes, proteins, types, classes for each farm and site
# also need to look at controls

# read in data
require(dplyr)
require(tibble)
require(stringr)
require(ggplot2)
dat <- read.table("https://github.com/EmilyB17/amr-brazil/blob/master/data/parsedRelativeAbundance.txt?raw=true", sep = "\t", header = TRUE)

# set theme
theme_set(theme_bw())

# 4 kinds of resistance
unique(dat$broadclass)

# for all of those 4 classes, 54 types of resistance (mostly drugs)
unique(dat$type)

# 186 unique proteins
unique(dat$protein)

# 798 unique genes
length(unique(dat$gene))

### filter only control samples
# pick out control samples
controls <- paste("G527_94_NoTemplate-DNAextraction2_S94",
                  "G527_96_NoTemplate-LibraryPrep_S96",
                  "G527_95_MockCommunity_S95",
                  "G527_93_NoTemplate-DNAextraction1_S93",
                  "G527_49_2019_DNA_S49", sep = "|")

# get dataframe of controls only
controlsDF <- relv %>% 
  filter(str_detect(name, controls)) %>% 
  # remove Alien sample
  filter(!str_detect(name, "G527_49_2019_DNA_S49"))

# get positive control
poscontrol <- controlsDF %>% filter(name == "G527_95_MockCommunity_S95") %>% 
  # remove zeros
  filter(relabun > 0)

# boxplot
ggplot(data = poscontrol, aes(x = broadclass, y = relabun)) +
  geom_boxplot() +
  labs(x = "Resistance Class", y = "Relative Abundance", title = "Mock Community")
# save to file
ggsave(filename = "C:/Users/emily/OneDrive - The Pennsylvania State University/Research/git/amr-brazil/data/plots/MockCommunityBoxplot.png",
       plot = last_plot())


# get DNA extraction negative controls
dna <- paste("G527_93_NoTemplate-DNAextraction1_S93",
             "G527_94_NoTemplate-DNAextraction2_S94", sep = "|")
dnacontrol <- controlsDF %>% 
  filter(str_detect(dna, name)) %>% 
  # get only present genes
  filter(relabun > 0) %>% 
  mutate(Control = factor(case_when(
    name %in% "G527_93_NoTemplate-DNAextraction1_S93" ~ 1,
    name %in% "G527_94_NoTemplate-DNAextraction2_S94" ~ 2
  )))

# boxplot
ggplot(data = dnacontrol, aes(x = broadclass, y = relabun)) +
  geom_boxplot() +
  labs(x = "Resistance Class", y = "Relative Abundance", title = "DNA Extraction Controls 1 & 2") +
  facet_wrap(~Control)
# save to file
ggsave(filename = "C:/Users/emily/OneDrive - The Pennsylvania State University/Research/git/amr-brazil/data/plots/DNAExtractionControlsBoxplot.png",
       plot = last_plot())


# get library prep control
libcontrol <- controlsDF %>% 
  filter(name == "G527_96_NoTemplate-LibraryPrep_S96") %>% 
  # filter out zeros
  filter(relabun > 0)

# boxplot
ggplot(data = libcontrol, aes(x = broadclass, y = relabun)) +
  geom_boxplot() +
  labs(x = "Resistance Class", y = "Relative Abundance", title = "Library Prep Control")
# save to file
ggsave(filename = "C:/Users/emily/OneDrive - The Pennsylvania State University/Research/git/amr-brazil/data/plots/LibPrepControlBoxplot.png",
       plot = last_plot())
