AMR Analysis
================
Emily Bean
March 12, 2020

## Analysis Structure

**Comparisons**

Pairwise comparisons between Farm 1 and Farm 2.  
Comparisons between body sites (rumen, feces, nasal swab) at each farm.

MEGARES database categories include: type of resistance, resistance
class within type (i.e. type of drug within Drug resistance), protein
class, and gene. Pairwise comparisons can be made for all of these but
it seems to make the most biological sense to compare genes, protein
classes, and resistance classes.

-----

Therefore, here is the proposed comparison structure:  
Genes: Farm 1 vs Farm 2  
Proteins: Farm 1 vs Farm 2  
Resistance class: Farm 1 vs Farm 2  
**Wilcox Rank Sum test**

Genes: SNP vs Rumen vs Feces (Farm 1 and Farm 2)  
Proteins: SNP vs Rumen vs Feces (Farm 1 and Farm 2)  
Resistance class: SNP vs Rumen vs Feces (Farm 1 and Farm 2)  
**Kruskal Wallis test**

-----

### Data Exploration

``` r
# in the entire dataset, there are 174 unique proteins in 3 classes
exp <- counts %>% 
  group_by(broadclass, type, protein) %>% 
  summarize()

# show unique proteins
#knitr::kable(exp)
```

``` r
# data is not normally distributed and is very zero-skewed
hist(counts$relabun, main = "Histogram of Relative Abundance",
     xlab = "Relative abundance", ylab = "Frequency")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

### Pairwise comparisons between farms

Wilcox Rank-Sum test for non-parametric pairwise test of relative
abundnace between farms

#### Genes

``` r
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
  group_by(farm, pattern, broadclass, type, protein, gene) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))

## write summary statistics to table
sum <- sigs %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR)) %>% 
  # add p value
  left_join(adj, by = "pattern") %>% 
  select(-c(pattern, pval, sig))
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/wilcox-genes.txt", sep = "\t", row.names = FALSE)
```

There are 5 genes that are significantly different between farms, which
is difficult to visualize.

``` r
# summary statistics on significant genes
sum <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "pattern") %>% 
  group_by(farm, pattern) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR, mean, sd), names_prefix = "farm_") %>% 
  left_join(adj, by = "pattern") %>% 
  select(-c(pval, sig)) %>% 
  mutate(pattern = str_replace_all(pattern, "\\|", "-"))

# print table
knitr::kable(sum, format = "markdown", escape = FALSE)
```

| pattern                                                             | median\_farm\_1 | median\_farm\_2 | IQR\_farm\_1 | IQR\_farm\_2 | mean\_farm\_1 | mean\_farm\_2 | sd\_farm\_1 | sd\_farm\_2 | padj |
| :------------------------------------------------------------------ | --------------: | --------------: | -----------: | -----------: | ------------: | ------------: | ----------: | ----------: | ---: |
| Drugs-MLS-Lincosamide\_nucleotidyltransferases-LNUB                 |       0.0000000 |               0 |    0.0147288 |            0 |     0.0608790 |     0.0000000 |   0.1654859 |   0.0000000 |    0 |
| Drugs-Multi-drug\_resistance-Multi-drug\_ABC\_efflux\_pumps-LSAE    |       0.0000000 |               0 |    0.0147288 |            0 |     0.0608790 |     0.0000000 |   0.1654859 |   0.0000000 |    0 |
| Drugs-Multi-drug\_resistance-Multi-drug\_RND\_efflux\_pumps-MUXB    |       0.0681038 |               0 |    0.1757803 |            0 |     0.1387449 |     0.0360132 |   0.2030432 |   0.1179626 |    0 |
| Drugs-Nucleosides-Streptothricin\_acetyltransferase-SAT             |       0.0000000 |               0 |    0.0340065 |            0 |     0.0530855 |     0.0043774 |   0.1653840 |   0.0287047 |    0 |
| Metals-Multi-metal\_resistance-Multi-metal\_RND\_efflux\_pumps-CZCP |       0.0000000 |               0 |    0.1358076 |            0 |     0.0983806 |     0.0167831 |   0.1826834 |   0.0796361 |    0 |

This plot shows all genes between farms; most decrease from farm 1 to
farm 2.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = gene)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different genes by farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")"))
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

The data is so zero-skewed that it can be misleading to plot just the
median and IQR, so below is a plot of whole dataset in the significant
genes.

``` r
## plot all data (not just  median)

# get significant data
allsig <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "pattern") %>% 
  mutate(farm = factor(farm))

# plot
ggplot(data = allsig, aes(x = farm, y = relabun,fill = farm)) +
  geom_violin() +
  geom_point(position = position_dodge2(width = 0.5), size = 0.3) +
  facet_grid(pattern~broadclass, scales = "free") +
  labs(x = "Farm", y = "Relative Abundance", title = "Relative abundance of significant genes")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### Proteins

``` r
outsig <- data.frame()
for(i in 1:length(proteins)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, protein == proteins[i]),
                   exact = FALSE)
  
  outsig[i, "protein"] <- proteins[i]
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
  left_join(counts, by = "protein") %>% 
  group_by(farm, broadclass, type, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))

## write summary statistics to table
sum <- sigs %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR)) %>% 
  # add p value
  left_join(adj, by = "protein") %>% 
  select(-c(pval, sig))
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/wilcox-proteins.txt", sep = "\t", row.names = FALSE)
```

There are 16 significantly different resistance proteins; the plot shows
that most are only present in one of the two farms.

``` r
# summary statistics on significant genes
sum <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "protein") %>% 
  group_by(farm, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR, mean, sd), names_prefix = "farm_") %>% 
  left_join(adj, by = "protein") %>% 
  select(-c(pval, sig))

# print table
knitr::kable(sum)
```

| protein                                         | median\_farm\_1 | median\_farm\_2 | IQR\_farm\_1 | IQR\_farm\_2 | mean\_farm\_1 | mean\_farm\_2 | sd\_farm\_1 | sd\_farm\_2 |      padj |
| :---------------------------------------------- | --------------: | --------------: | -----------: | -----------: | ------------: | ------------: | ----------: | ----------: | --------: |
| 16S\_rRNA\_methyltransferases                   |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0309216 |     0.0212610 |   0.1323115 |   0.1279766 | 0.0425000 |
| 23S\_rRNA\_methyltransferases                   |       0.0000000 |               0 |    0.0509730 |    0.0000000 |     0.0796005 |     0.0787008 |   0.1827057 |   0.2004140 | 0.0000000 |
| Aminoglycoside\_O-phosphotransferases           |       0.0097070 |               0 |    0.0703389 |    0.0230695 |     0.0827141 |     0.0363395 |   0.1830912 |   0.0919702 | 0.0000000 |
| Class\_A\_betalactamases                        |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0313945 |     0.0175315 |   0.1525016 |   0.1099037 | 0.0130769 |
| Class\_B\_betalactamases                        |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0341897 |     0.0183991 |   0.1416332 |   0.1075469 | 0.0000000 |
| Drug\_and\_biocide\_RND\_efflux\_pumps          |       0.0000000 |               0 |    0.0132188 |    0.0000000 |     0.0558194 |     0.0239299 |   0.1649132 |   0.1047371 | 0.0000000 |
| MLS\_resistance\_ABC\_efflux\_pumps             |       0.0000000 |               0 |    0.1280914 |    0.0521572 |     0.0976780 |     0.0802149 |   0.1839002 |   0.1805745 | 0.0000000 |
| Multi-biocide\_ABC\_efflux\_pump                |       0.0000000 |               0 |    0.0096316 |    0.0000000 |     0.0570699 |     0.0363381 |   0.1673507 |   0.1255296 | 0.0000000 |
| Multi-biocide\_RND\_efflux\_pump                |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0467279 |     0.0176145 |   0.1611198 |   0.1048535 | 0.0000000 |
| Multi-drug\_ABC\_efflux\_pumps                  |       0.0000000 |               0 |    0.0534582 |    0.0000000 |     0.0655053 |     0.0749633 |   0.1466797 |   0.1954580 | 0.0340000 |
| Multi-drug\_RND\_efflux\_pumps                  |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0590809 |     0.0183988 |   0.1722372 |   0.1051113 | 0.0000000 |
| Mupirocin-resistant\_isoleucyl-tRNA\_synthetase |       0.0024421 |               0 |    0.0577595 |    0.0160325 |     0.0772348 |     0.0640532 |   0.1666730 |   0.1700721 | 0.0130769 |
| Para-aminosalicylic\_acid\_resistant\_mutant    |       0.0000000 |               0 |    0.1087074 |    0.0000000 |     0.0853541 |     0.0677343 |   0.1578868 |   0.2025623 | 0.0000000 |
| Streptothricin\_acetyltransferase               |       0.0000000 |               0 |    0.0340065 |    0.0000000 |     0.0530855 |     0.0043774 |   0.1653840 |   0.0287047 | 0.0000000 |
| Tellurium\_resistance\_protein                  |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0372383 |     0.0001202 |   0.1762108 |   0.0013655 | 0.0340000 |
| Tetracycline\_resistance\_MFS\_efflux\_pumps    |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0371465 |     0.0273027 |   0.1429913 |   0.1218001 | 0.0000000 |

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = protein)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")"))
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

There are 16 proteins with significant differences, so that’s difficult
to visualize on one plot.

``` r
## plot all data (not just  median)

# get significant data
allsig <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "protein") %>% 
  mutate(farm = factor(farm))

# plot
ggplot(data = allsig, aes(x = farm, y = relabun, fill = farm)) +
  geom_violin() +
  geom_point(position = position_dodge2(width = 0.5), size = 0.3) +
  facet_wrap(~protein) +
  labs(x = "Farm", y = "Relative Abundance", title = "Relative abundance of significant proteins")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

Visualize with density plots; this shows that most of the proteins have
marginally higher relative abundance in farm 1.

``` r
# get significant data
sigrel <- counts %>% 
  semi_join(filter(adj, padj < 0.05), by = "protein")

# plot with ggpubr
ggdensity(sigrel, x = "relabun", y = "..count..", add = "median", 
          color = "farm", fill = "farm", rug = TRUE) +
  facet_wrap(~protein, scales = "free")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

#### Resistance Type

``` r
outsig <- data.frame()
for(i in 1:length(types)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, type == types[i]),
                   exact = FALSE)
  
  outsig[i, "type"] <- types[i]
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
  left_join(counts, by = "type") %>% 
  group_by(farm, broadclass, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))

## write summary statistics to table
sum <- sigs %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR)) %>% 
  # add p value
  left_join(adj, by = "type") %>% 
  select(-c(pval, sig))
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/wilcox-types.txt", sep = "\t", row.names = FALSE)
```

There are 14 significantly different resistance types; again, the plot
shows that most are only present in one of the two farms.

``` r
# summary statistics on significant genes
sum <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "type") %>% 
  group_by(farm, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR, mean, sd), names_prefix = "farm_") %>% 
  left_join(adj, by = "type") %>% 
  select(-c(pval, sig))

# print table

knitr::kable(sum)
```

| type                                       | median\_farm\_1 | median\_farm\_2 | IQR\_farm\_1 | IQR\_farm\_2 | mean\_farm\_1 | mean\_farm\_2 | sd\_farm\_1 | sd\_farm\_2 |      padj |
| :----------------------------------------- | --------------: | --------------: | -----------: | -----------: | ------------: | ------------: | ----------: | ----------: | --------: |
| Aminoglycosides                            |       0.0000336 |               0 |    0.1141168 |    0.0680577 |     0.1139073 |     0.0860095 |   0.2286553 |   0.1787539 | 0.0000000 |
| betalactams                                |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0314342 |     0.0121345 |   0.1508087 |   0.0896989 | 0.0000000 |
| Drug\_and\_biocide\_resistance             |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0421138 |     0.0181552 |   0.1549526 |   0.1032013 | 0.0000000 |
| Glycopeptides                              |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0411861 |     0.0318189 |   0.1466645 |   0.1353008 | 0.0000000 |
| MLS                                        |       0.0000000 |               0 |    0.0345729 |    0.0000000 |     0.0794274 |     0.0657229 |   0.1913855 |   0.1783858 | 0.0000000 |
| Multi-biocide\_resistance                  |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0454174 |     0.0158222 |   0.1639115 |   0.0930139 | 0.0000000 |
| Multi-drug\_resistance                     |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0560412 |     0.0488867 |   0.1613654 |   0.1654866 | 0.0000000 |
| Multi-metal\_resistance                    |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0386289 |     0.0278699 |   0.1478089 |   0.1279490 | 0.0274615 |
| Mupirocin                                  |       0.0024421 |               0 |    0.0577595 |    0.0160325 |     0.0772348 |     0.0640532 |   0.1666730 |   0.1700721 | 0.0046364 |
| Mycobacterium\_tuberculosis-specific\_Drug |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0381180 |     0.0266663 |   0.1400829 |   0.1318297 | 0.0000000 |
| Nucleosides                                |       0.0000000 |               0 |    0.0340065 |    0.0000000 |     0.0530855 |     0.0043774 |   0.1653840 |   0.0287047 | 0.0000000 |
| Tellurium\_resistance                      |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0372383 |     0.0001202 |   0.1762108 |   0.0013655 | 0.0127500 |
| Tetracyclines                              |       0.0000000 |               0 |    0.0021114 |    0.0000000 |     0.0699916 |     0.0549477 |   0.1853600 |   0.1515760 | 0.0000000 |
| Trimethoprim                               |       0.0000000 |               0 |    0.0000000 |    0.0000000 |     0.0228799 |     0.0058738 |   0.1435111 |   0.0698217 | 0.0291429 |

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different types by farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")"))
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

``` r
## plot all data (not just  median)

# get significant data
allsig <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "type") %>% 
  mutate(farm = factor(farm))

# plot
ggplot(data = allsig, aes(x = farm, y = relabun, fill = farm)) +
  geom_violin() +
  geom_point(position = position_dodge2(width = 0.5), size = 0.3) +
  facet_wrap(~type) +
  labs(x = "Farm", y = "Relative Abundance", title = "Relative abundance of significant resistance type")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

#### Resistance Class

``` r
outsig <- data.frame()
for(i in 1:length(classes)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, broadclass == classes[i]),
                   exact = FALSE)
  
  outsig[i, "class"] <- classes[i]
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
  left_join(counts, by = c("class" = "broadclass")) %>% 
  group_by(farm, class) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))

## write summary statistics to table
sum <- sigs %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR)) %>% 
  # add p value
  left_join(adj, by = "class") %>% 
  select(-c(pval, sig))
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/wilcox-class.txt", sep = "\t", row.names = FALSE)
```

There are 4 significantly different resistance types.

``` r
# summary statistics on significant genes
sum <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = c("class" = "broadclass")) %>% 
  group_by(farm, class) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = farm, values_from = c(median, IQR, mean, sd), names_prefix = "farm_") %>% 
  left_join(adj, by = "class") %>% 
  select(-c(pval, sig))


# print table
knitr::kable(sum)
```

| class          | median\_farm\_1 | median\_farm\_2 | IQR\_farm\_1 | IQR\_farm\_2 | mean\_farm\_1 | mean\_farm\_2 | sd\_farm\_1 | sd\_farm\_2 |  padj |
| :------------- | --------------: | --------------: | -----------: | -----------: | ------------: | ------------: | ----------: | ----------: | ----: |
| Biocides       |               0 |               0 |            0 |            0 |     0.0419351 |     0.0120948 |   0.1627921 |   0.0782484 | 0.000 |
| Drugs          |               0 |               0 |            0 |            0 |     0.0577865 |     0.0437536 |   0.1728195 |   0.1477930 | 0.000 |
| Metals         |               0 |               0 |            0 |            0 |     0.0370398 |     0.0255965 |   0.1462375 |   0.1261878 | 0.001 |
| Multi-compound |               0 |               0 |            0 |            0 |     0.0412593 |     0.0171020 |   0.1558914 |   0.0971461 | 0.000 |

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = class)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different classes by farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")"))
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

``` r
## plot all data (not just  median)

# get significant data
allsig <- adj %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = c("class" = "broadclass")) %>% 
  mutate(farm = factor(farm))

# plot
ggplot(data = allsig, aes(x = farm, y = relabun, fill = farm)) +
  geom_violin() +
  geom_point(position = position_dodge2(width = 0.5), size = 0.3) +
  facet_wrap(~type) +
  labs(x = "Farm", y = "Relative Abundance", title = "Relative abundance of significant resistance class")
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

### Comparisons between body sites

Kruskal-Wallis test for non-parametric comparisons between the 3 body
sites, stratified by farm.

#### Genes

``` r
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

## write summary statistics to table
sum <- sigs %>% 
  left_join(adj, by = c("pattern", "farm")) %>% 
  select(-c(pval, sig)) %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR))
  
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/kruskal-genes.txt", sep = "\t", row.names = FALSE)
```

There are 345 significant genes between body sites. The p values and
summary statistics are not printed since there are so many.

``` r
# summary statistics on significant genes
sumfarm1 <- adj %>% 
  filter(sig == "TRUE" & farm == 1) %>% 
  left_join(counts, by = c("pattern", "farm")) %>% 
  group_by(site, pattern) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "pattern") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm1)

sumfarm2 <- adj %>% 
  filter(sig == "TRUE" & farm == 2) %>% 
  left_join(counts, by = c("pattern", "farm")) %>% 
  group_by(site, pattern) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "pattern") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm2)
```

This plot visualizes the differences; most are higher in feces, then
rumen, then SNP. However, one gene is higher in SNP.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = pattern)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different genes by site & farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")")) +
  facet_wrap(~ farm)
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

#### Proteins

``` r
## Stratify for each farm
counts <- counts %>% mutate(farm = factor(farm))
farms <- levels(counts$farm)
# create empty dataframe to fill within the loop
outsig <- data.frame()

for(j in 1:length(farms)) {
  
  # iterate through all 666 genes
  for(i in 1:length(proteins)) {
    
    # perform Kruskal-Wallis test
    k <- kruskal.test(relabun ~ site, data = filter(counts, protein == proteins[i] &
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
                               protein = proteins[i],
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
  left_join(counts, by = c("protein", "farm")) %>% 
  group_by(site, broadclass, type, farm, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site),
         farm = factor(farm))

## write summary statistics to table
## write summary statistics to table
sum <- sigs %>% 
  left_join(adj, by = c("protein", "farm")) %>% 
  select(-c(pval, sig)) %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR)) %>% 
  arrange(farm)
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/kruskal-proteins.txt", sep = "\t", row.names = FALSE)
```

There are 117 significant proteins in the dataset.

``` r
# summary statistics on significant genes
sumfarm1 <- adj %>% 
  filter(sig == "TRUE" & farm == 1) %>% 
  left_join(counts, by = c("protein", "farm")) %>% 
  group_by(site, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "protein") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm1)

sumfarm2 <- adj %>% 
  filter(sig == "TRUE" & farm == 2) %>% 
  left_join(counts, by = c("protein", "farm")) %>% 
  group_by(site, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "protein") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm2)
```

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = protein)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by site & farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")")) +
  facet_wrap(~ farm)
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

#### Resistance Type

``` r
## Stratify for each farm
counts <- counts %>% mutate(farm = factor(farm))
farms <- levels(counts$farm)
# create empty dataframe to fill within the loop
outsig <- data.frame()

for(j in 1:length(farms)) {
  
  # iterate through all 666 genes
  for(i in 1:length(types)) {
    
    # perform Kruskal-Wallis test
    k <- kruskal.test(relabun ~ site, data = filter(counts, type == types[i] &
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
                               type = types[i],
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
  left_join(counts, by = c("type", "farm")) %>% 
  group_by(site, broadclass, type, farm) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site),
         farm = factor(farm))

## write summary statistics to table
## write summary statistics to table
sum <- sigs %>% 
  left_join(adj, by = c("type", "farm")) %>% 
  select(-c(pval, sig)) %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR)) %>% 
  arrange(farm)
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/kruskal-types.txt", sep = "\t", row.names = FALSE)
```

There are 44 significant proteins in the dataset.

``` r
# summary statistics on significant genes
sumfarm1 <- adj %>% 
  filter(sig == "TRUE" & farm == 1) %>% 
  left_join(counts, by = c("type", "farm")) %>% 
  group_by(site, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "type") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm1)

sumfarm2 <- adj %>% 
  filter(sig == "TRUE" & farm == 2) %>% 
  left_join(counts, by = c("type", "farm")) %>% 
  group_by(site, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "type") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm2)
```

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different types by site & farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")")) +
  facet_wrap(~ farm)
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

#### Class

``` r
## Stratify for each farm
counts <- counts %>% mutate(farm = factor(farm))
farms <- levels(counts$farm)
# create empty dataframe to fill within the loop
outsig <- data.frame()

for(j in 1:length(farms)) {
  
  # iterate through all 666 genes
  for(i in 1:length(classes)) {
    
    # perform Kruskal-Wallis test
    k <- kruskal.test(relabun ~ site, data = filter(counts, broadclass == classes[i] &
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
                               class = classes[i],
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
  left_join(counts, by = c("class" = "broadclass", "farm")) %>% 
  group_by(site, class, farm) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site),
         farm = factor(farm))

## write summary statistics to table
## write summary statistics to table
sum <- sigs %>% 
  left_join(adj, by = c("class", "farm")) %>% 
  select(-c(pval, sig)) %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR)) %>% 
  arrange(farm)
#write.table(sum, "/Users/epb5360/git/amr-brazil/data/results-tables/kruskal-class.txt", sep = "\t", row.names = FALSE)
```

There are 4 significant classes in the dataset.

``` r
# summary statistics on significant genes
sumfarm1 <- adj %>% 
  filter(sig == "TRUE" & farm == 1) %>% 
  left_join(counts, by = c("class" = "broadclass", "farm")) %>% 
  group_by(site, class) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "class") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm1)

sumfarm2 <- adj %>% 
  filter(sig == "TRUE" & farm == 2) %>% 
  left_join(counts, by = c("class" = "broadclass", "farm")) %>% 
  group_by(site, class) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun),
            mean = mean(relabun),
            sd = sd(relabun)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = site, values_from = c(median, IQR, mean, sd)) %>% 
  left_join(adj, by = "class") %>% 
  select(-c(pval, sig))

# print table
#knitr::kable(sumfarm2)
```

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = class)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different classes by site & farm",
          subtitle = paste0("n = (", length(which(adj$sig == TRUE)), ")")) +
  facet_wrap(~ farm)
```

<img src="statAnalysis_files/figure-gfm/unnamed-chunk-31-1.png" style="display: block; margin: auto;" />
