
## UPSET PLOTS

## ---- Data ----

# install & require packages

require(tidyverse)

if(!require(ComplexHeatmap)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ComplexHeatmap")
  require(ComplexHeatmap)
}

if(!require(UpSetR)) {
  install.packages("UpSetR")
  require(UpSetR)
}

# read data - presence/absence
dat <- read.table("https://github.com/EmilyB17/amr-brazil/blob/master/data/parsedPresenceAbsence.txt?raw=true",
                  sep = "\t", header = TRUE)

## ---- Colors ----

# get the pre-determined colors
source("./R/RColorBrewer.R")

## ---- Wrangle ----

# get data into binary matrix
datm <- dat %>% 
  select(pattern, name, presence) %>% 
  pivot_wider(names_from = pattern, values_from = presence) %>% 
  column_to_rownames(var = "name")


## ---- Drugs ----

# first group data by drug types and sum presence/absence
drug <- dat %>% 
  filter(broadclass == "Drugs") %>%
  mutate(name = paste(name, farm, sep = "_farm")) %>% 
  group_by(name, type) %>% 
  summarize(sumpres = sum(presence)) %>% 
  # determine which have presence/absence
  mutate(presence = case_when(
    sumpres != 0 ~ "1",
    sumpres == 0 ~ "0")) %>%
  mutate(presence = as.numeric(presence))

# string parsing to make type look better
drug$type1 <- str_replace_all(drug$type, "_", " ")        
drug$type1 <- str_replace_all(drug$type1, 
                              "Mycobacterium tuberculosis-specific Drug",
                              "Mycobacterium TB-Specific")

drug1 <- drug %>% 
  select(-type) %>% 
  # pivot wider
  select(-sumpres) %>% 
  pivot_wider(names_from = type1, values_from = presence) %>% 
  ungroup() %>% 
  column_to_rownames(var = "name") 

## plot each farm separately
f1 <- drug1 %>% 
  filter(grepl("_farm1", rownames(drug1)))

# make the plot using UpsetR package
upset(f1, 
      # number of sets to plot - by default, only shows a few
      nsets = 26,
      # orders by intersection frequency
      order.by = "freq", group.by = "sets",
      # name labels
      sets.x.label = "Total Samples", mainbar.y.label = "Intersection Size",
      # ratio of barplot to matrix
      mb.ratio = c(0.3, 0.5),
      # don't show numbers above the intersections
      show.numbers = FALSE,
      # show the total length of each set
      set_size.show = TRUE,
      
      ## COLORS
      # color of matrix shade and transparency (alpha)
      shade.color = accentcols[1], shade.alpha = 0.5,
      # barplot color
      main.bar.color = unname(farmcols[1]),
      # set bar color
      sets.bar.color = "grey"
)

## plot each farm separately
f2 <- drug1 %>% 
  filter(grepl("_farm2", rownames(drug1)))

# make the plot using UpsetR package
upset(f2, 
      # number of sets to plot - by default, only shows a few
      nsets = ncol(f2),
      # orders by intersection frequency
      order.by = "freq", #group.by = "sets",
      # name labels
      sets.x.label = "Total Samples", mainbar.y.label = "Intersection Size",
      # ratio of barplot to matrix
      mb.ratio = c(0.3, 0.5),
      # don't show numbers above the intersections
      show.numbers = FALSE,
      # show the total length of each set
      set_size.show = TRUE,
      
      ## COLORS
      # color of matrix shade and transparency (alpha)
      shade.color = accentcols[1], shade.alpha = 0.5,
      # barplot color
      main.bar.color = unname(farmcols[2]),
      # set bar color
      sets.bar.color = "grey"
)


## ---- all data with pattern in columns ----

# need: pattern as rowname, spread by sites

dat1 <- dat %>% 
  # make column for farm & site
  mutate(farm1 = case_when(
    farm %in% 1 ~ "A",
    farm %in% 2 ~ "B"
  )) %>% 
  mutate(farmsite = paste0("farm", farm1, "-", site)) %>% 
  # summarize by farmsite
  group_by(farmsite, pattern) %>% 
  summarize(sumpres = sum(presence))  %>% 
  ungroup() %>% 
  # make binary
  mutate(pres = case_when(
    sumpres == 0 ~ "0",
    sumpres != 0 ~ "1"
  )) %>% 
  mutate(pres = as.numeric(pres)) %>% 
  # spread by farmsite
  select(-sumpres) %>% 
  pivot_wider(names_from = farmsite, values_from = pres) %>% 
  column_to_rownames(var = "pattern")

# make the upset plot
# make the plot using UpsetR package
upset(dat1, sets = c("farmA-SNP", "farmA-rumen", "farmA-feces",
                     "farmB-SNP", "farmB-rumen", "farmB-feces"), keep.order = TRUE,
      ## COLORS
      # color of matrix shade and transparency (alpha)
      shade.color = accentcols[1], shade.alpha = 0.5,
      # barplot color
      main.bar.color = accentcols[2],
      # set bar color
      sets.bar.color = "grey")


