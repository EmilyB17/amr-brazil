
# Making plots for the manuscript

## ---- data and packages ----

require(tidyverse)
require(RColorBrewer)
require(ggpubr)

# presence/absence
pres <- read.table("https://github.com/EmilyB17/amr-brazil/blob/master/data/parsedPresenceAbsence.txt?raw=true",
                   sep = "\t", header = TRUE) %>% 
  # change farm column; should be Farm A and B
  mutate(farm = case_when(
    farm %in% 1 ~ "A",
    farm %in% 2 ~ "B"
  ))

# source the RColorBrewer code
source("./R/RColorBrewer.R")

# set theme
theme_set(theme_bw())

# define function for standard error
se <- function(x) {sqrt(stats::var(x)/length(x))}

## ---- barplot for sample distribution ----

## plot of heterogeneity of the 546 genes

# get range of samples each gene is in
wide <- pres %>% 
  select(name, pattern, presence) %>% 
  pivot_wider(names_from = pattern, values_from = presence) %>% 
  column_to_rownames(var = "name")

# get column sums (range of prevalence)
sum <- colSums(wide)
sumdf <- enframe(sum, name = "pattern", value = "sum") %>% 
  # add bins of ranges to make easier to visualize
  mutate(range = factor(case_when(
    sum > 0 & sum < 9 ~ "fewer than 10%",
    sum >= 9 & sum < 46 ~ "10-50%",
    sum >= 46 & sum < 73 ~ "51-80%",
    sum >= 73 & sum < 91 ~ "81-98%",
    sum == 91 ~ "100%"
  ), ordered = TRUE, levels = c("fewer than 10%", "10-50%", "51-80%", 
                                "81-98%", "100%"))) %>% 
  group_by(range) %>% 
  summarize(sumrange = length(range))

# plot showing distribution of gene prevalence
plot <- ggplot(data = sumdf) +
  geom_col(aes(x = range, y = sumrange, fill = range)) +
  labs(x = "Sample Percentage", y = "Number of genes", 
       # change the legend title
       fill = "Sample Percentage") +
  scale_fill_manual(values = accentcols) +
  geom_text(aes(x = 1, y = 325, label = "317")) +
  geom_text(aes(x = 2, y = 282, label = "274")) +
  geom_text(aes(x = 3, y = 67, label = "59")) +
  geom_text(aes(x = 4, y = 23, label = "15")) +
  geom_text(aes(x = 5, y = 9, label = "1"))

# save
#ggsave(filename = "./data/plots/sample-distribution-barplot.png", plot = plot, dpi = 300)


## ---- barplot for gene prevalence by farm ----

## we have to be careful here because the farms have unequal sample sizes
# therefore, we need to compare "by sample", or "average", instead of counts

# calculate per sample average count
tot <- pres %>% 
  group_by(farm, broadclass) %>% 
  summarize(sum = sum(presence)) %>% 
            #ster = se(presence),
            #sd = sd(presence)) 
            
  mutate(avg = case_when(
    # farms have unequal sample sizes, so divide by sample size
    farm %in% "A" ~ sum / 48,
    farm %in% "B" ~ sum / 43
  )) %>% 
  ungroup() %>% 
  mutate(Farm = factor(farm))

# plot
p <- ggplot(data = tot, aes(x = broadclass, y = avg, fill = Farm)) +
  geom_col(position = position_dodge(1)) +
  # customize barplot colors
  scale_fill_manual(values = farmcolscorrect) +
  labs(x = "Resistance Class", y = "Average Gene Presence") 
# save
#ggsave(filename = "./data/plots/class-prevalence-byfarm-barplot.png", plot = p, dpi = 300)

