

## run normalizationDecontam RMD
counts = relfiltv

patterns <- unique(relfiltv$pattern)
# Wilcox rank-sum for pairwise differences between farms

# create empty df to fill in the loop
outsig <- data.frame()
# loop through all 666 genes
for(i in 1:length(patterns)) {
  
  # perform wilcox test
  w <- wilcox.test(relabun ~ farm, data = filter(counts, pattern == patterns[i]),
                   exact = FALSE)
  
  # get output p value
  outsig[i, "pattern"] <- patterns[i]
  outsig[i, "pval"] <- round(w$p.value, 3)
  
}

# adjust p values for false discovery rate
adj <- outsig %>% 
  mutate(
    # adjust p values with FDR
    padj = p.adjust(pval, method = "fdr", n = length(outsig$pval)),
    # Boolean for significance
    sig = case_when(
      padj < 0.05 ~ TRUE,
      padj >= 0.05 ~ FALSE))

# get only significant variables and add median + IQR
sigs <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "pattern") %>% 
  # add presence/absence
  mutate(presence = case_when(
    relabun == 0 ~ 0,
    relabun != 0 ~ 1
  )) %>% 
  group_by(farm, pattern, broadclass, type, protein, gene) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            prestotal = sum(presence)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))

## write summary statistics to table
sum <- sigs %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR, prestotal)) %>% 
  # add p value
  left_join(adj, by = "pattern") %>% 
  select(-c(pattern, pval, sig))
