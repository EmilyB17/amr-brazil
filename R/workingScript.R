# stat analysis working

require(ggpubr)


# non-parametric for one-way ANOVA

kruskal.test(relabun ~ site, data = counts)


# post-hoc for kruskal wallis
pairwise.wilcox.test(counts$relabun, counts$farm,
                     p.adjust.method = "bonferroni")

# this is the t test
wilcox.test(relabun ~ farm, data = counts)


# ---- WilcoxRank ----

## WILCOX RANK SUM FOR COMPARISONS BY FARM
counts <- counts %>% 
  mutate(farm = factor(farm, ordered = TRUE, levels = c("1", "2")),
         animal = factor(animal))

# difference in resistance classes 
wilcox.test(relabun ~ farm, data = filter(counts, broadclass == "Biocides"))
boxplot(relabun ~ broadclass, data = counts)

w <- wilcox.test(relabun ~ farm, data = filter(counts, broadclass == "Biocides"))

classes <- unique(counts$broadclass)
proteins <- unique(counts$protein)
genes <- unique(counts$gene)
patterns <- unique(counts$pattern)

outsig <- data.frame()
for(i in 1:length(patterns)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, pattern == patterns[i]),
                   exact = FALSE)
  
  outsig[i, "pattern"] <- patterns[i]
  outsig[i, "pval"] <- round(w$p.value, 3)
  
  if(w$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

s <- outsig %>% filter(sig == TRUE)
# get only significant variables and add median + IQR
sigs <- outsig %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "pattern") %>% 
  group_by(farm, pattern) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))



  
# for genes that are present in both farms, see which farm has more abundance
sum <- outsig %>% 
  filter(sig == TRUE) %>% 
  left_join(counts, by = "pattern")

# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = pattern)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different types by farm (n = 18)")

## ---- Kruskal-Wallis ----

## BOTH FARMS
outsig <- data.frame()
for(i in 1:length(types)) {
  
  k <- kruskal.test(relabun ~ site, data = filter(counts, type == types[i]))
  
  outsig[i, "type"] <- types[i]
  outsig[i, "pval"] <- round(k$p.value, 3)
  
  if(k$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

# get significant variables
sigs <- outsig %>% filter(sig == "TRUE") %>%  
  left_join(counts, by = "type") %>% 
  group_by(site, broadclass, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site))

# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by site (n = 46)")
