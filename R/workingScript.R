# stat analysis working

require(ggpubr)


# non-parametric for one-way ANOVA

kruskal.test(relabun ~ site, data = counts)


# post-hoc for kruskal wallis
pairwise.wilcox.test(counts$relabun, counts$farm,
                     p.adjust.method = "bonferroni")

# this is the t test
wilcox.test(relabun ~ farm, data = counts)


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

# get only significant variables and add median + IQR
sigs <- outsig %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "pattern") %>% 
  group_by(farm, pattern, broadclass, type, protein, gene) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))


# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = gene)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ggtitle("Significantly different genes by farm (n = 85)")
  