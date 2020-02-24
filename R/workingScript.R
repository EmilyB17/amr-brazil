
## WORKING SCRIPT: CUSTOM BLAST DATABASE PARSING

require(dplyr)
require(tidyr)
require(ggplot2)

# set wd
setwd("C:/Users/emily/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/emily/megares")

# read results dataframe
dat <- read.table("out.txt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)

# set column names
colnames(dat) <- c("readID", "fullID", "fullID2", "qcov", "qcov2")

# somehow we have extra columns
dat$fullID2 <- NULL
dat$qcov2 <- NULL

# how many unique reads?
length(unique(dat$readID)) # 2259

## string parsing
fulldat <- dat %>% 
  mutate(megID = sapply(strsplit(dat$fullID, split = "\\|"), `[`, 1),
         class = sapply(strsplit(dat$fullID, split = "\\|"), `[`, 2),
         drug = sapply(strsplit(dat$fullID, split = "\\|"), `[`, 3),
         protein = sapply(strsplit(dat$fullID, split = "\\|"), `[`, 4),
         gene = sapply(strsplit(dat$fullID, split = "\\|"), `[`, 5))

# since we don't really care about the specific reads, get summary for the whole sample
sumdat <- fulldat %>% 
  group_by(megID, class, drug, protein, gene) %>% 
  summarize(count = length(gene)) %>% 
  ungroup()

# split into two demo samples
sampdat <- sumdat %>% 
  mutate(sample = rep(c(1, 2), each = (254/2)))

## exploratory data analysis

# this looks better if we subset the top contributing proteins
subdat <- sampdat[order(sampdat$count, decreasing = TRUE), ]
subdat <- subdat[1:20, ]

# make a bar plot
ggplot(data = subdat, aes(x = drug, y = count, fill = drug)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~sample)

# make a boxplot
# this is not well representative of the data; some drug classes have more genes than others
ggplot(data = subdat, aes(x = drug, y = count, fill = drug)) +
  geom_boxplot() +
  facet_wrap(~sample) 

# what about the genes?
ggplot(data = subdat, aes(x = gene, y = count, fill = drug)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~sample)

# Hypothesis testing: ANOVA
mod <- aov(count ~ sample + class + drug + protein + gene, data = sampdat)
summary(mod)

## ---- metadata ----

# This is just annotations - I already have it
met <- read.csv("megares_modified_annotations_v2.00.csv")

# this shows what genes are coming from which database
# this likely explains the double-counting of some genes
met1 <- read.csv("megares_to_external_header_mappings_v2.00.csv",
                 stringsAsFactors = FALSE)

## test the theory - do the double-counted genes come from different databases?

# get unique megIDS
megids <- unique(fulldat$megID)

# parse unique megids from met1 
dbs <- met1 %>% 
  mutate(megID = sapply(strsplit(met1$MEGARes_v2_header, split = "\\|"), `[`, 1))

# get matching megIDS
matched <- dbs[dbs$megID %in% megids, ]

# CONFIRMED - each unique mutation described in CARD is collapsed into one category in MEGARES but assigned a different megID
