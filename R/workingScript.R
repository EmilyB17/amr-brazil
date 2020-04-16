

## Stratify for each farm
counts <- counts %>% mutate(farm = factor(farm))
farms <- levels(counts$farm)
# create empty dataframe to fill within the loop
outsig <- data.frame()

for(j in 1:length(farms)) {
  
  # iterate through all 666 genes
  for(i in 1:length(patterns)) {
    
    # perform Kruskal-Wallis test
    k <- kruskal.test(relabun ~ site, data = filter(counts, pattern == patterns[i] &
                                                      farm == farms[j]))
    
   
    # some farms do not have the gene present; mark as 'absent'
    if(is.na(k$p.value)) {
      
      pvaldat <- NA
      
    } else {
      
      pvaldat <- round(k$p.value, 3)
      
    } 
        
    # collect variables to output DF
    outsig <- rbind(outsig,
                    data.frame(farm = farms[j],
                               pattern = patterns[i],
                               pval = pvaldat))
  }
}

# adjust p vals for multiple comparisons with false discovery rate
adj <- outsig %>% 
  mutate(
    # adjust p values with FDR
    padj = p.adjust(pval, method = "fdr", n = length(outsig$pval)),
    # Boolean for significance
    sig = case_when(
      padj < 0.05 ~ TRUE,
      padj >= 0.05 ~ FALSE))
 
# get genes that had a significant difference
sigs <- adj %>% filter(sig == "TRUE") %>%  
  left_join(counts, by = c("pattern", "farm")) %>% 
  group_by(site, broadclass, farm, pattern) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site),
         farm = factor(farm))

# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by site & farm (total n = 45)") +
  facet_wrap(~ farm)


