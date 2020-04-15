AMR Analysis
================
Emily Bean
March 12, 2020

Analysis Structure
------------------

**Comparisons**

Pairwise comparisons between Farm 1 and Farm 2.
Comparisons between body sites (rumen, feces, nasal swab) at each farm.

MEGARES database categories include: type of resistance, resistance class within type (i.e. type of drug within Drug resistance), protein class, and gene. Pairwise comparisons can be made for all of these but it seems to make the most biological sense to compare genes, protein classes, and resistance classes.

------------------------------------------------------------------------

Therefore, here is the proposed comparison structure:
Genes: Farm 1 vs Farm 2
Proteins: Farm 1 vs Farm 2
Resistance class: Farm 1 vs Farm 2
**Wilcox Rank Sum test**

Genes: SNP vs Rumen vs Feces (Farm 1 and Farm 2)
Proteins: SNP vs Rumen vs Feces (Farm 1 and Farm 2)
Resistance class: SNP vs Rumen vs Feces (Farm 1 and Farm 2)
**Kruskal Wallis test**

------------------------------------------------------------------------

### Data Exploration

``` r
# in the entire dataset, there are 174 unique proteins in 3 classes
exp <- counts %>% 
  group_by(broadclass, type, protein) %>% 
  summarize()

# show unique proteins
knitr::kable(exp)
```

| broadclass     | type                                                | protein                                                       |
|:---------------|:----------------------------------------------------|:--------------------------------------------------------------|
| Biocides       | Acetate\_resistance                                 | Acetate\_resistance\_protein                                  |
| Biocides       | Acid\_resistance                                    | Acid\_resistance\_protein                                     |
| Biocides       | Acid\_resistance                                    | Acid\_resistance\_regulator                                   |
| Biocides       | Biguanide\_resistance                               | Biguanide\_cation\_efflux                                     |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_ABC\_efflux\_pump                              |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_MATE\_efflux\_pump                             |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_MFS\_efflux\_pump                              |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_resistance\_protein                            |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_resistance\_regulator                          |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_RND\_efflux\_pump                              |
| Biocides       | Multi-biocide\_resistance                           | Multi-biocide\_RND\_efflux\_regulator                         |
| Biocides       | Paraquat\_resistance                                | Paraquat\_resistance\_protein                                 |
| Biocides       | Peroxide\_resistance                                | peroxide\_resistance\_protein                                 |
| Biocides       | Peroxide\_resistance                                | Peroxide\_resistance\_stress\_protein                         |
| Biocides       | Phenolic\_compound\_resistance                      | Phenolic\_RND\_efflux\_pump                                   |
| Biocides       | Phenolic\_compound\_resistance                      | Triclosan-resistant\_mutation                                 |
| Biocides       | Quaternary\_Ammonium\_Compounds\_(QACs)\_resistance | QAC-resistant\_protein\_UDP\_glucose\_4\_epimerase            |
| Biocides       | Quaternary\_Ammonium\_Compounds\_(QACs)\_resistance | QAC\_efflux\_pump                                             |
| Drugs          | Aminocoumarins                                      | Aminocoumarin-resistant\_DNA\_topoisomerases                  |
| Drugs          | Aminocoumarins                                      | Aminocoumarin\_efflux\_pump                                   |
| Drugs          | Aminoglycosides                                     | 16S\_rRNA\_methyltransferases                                 |
| Drugs          | Aminoglycosides                                     | Aminoglycoside-resistant\_16S\_ribosomal\_subunit\_protein    |
| Drugs          | Aminoglycosides                                     | Aminoglycoside\_efflux\_pumps                                 |
| Drugs          | Aminoglycosides                                     | Aminoglycoside\_N-acetyltransferases                          |
| Drugs          | Aminoglycosides                                     | Aminoglycoside\_O-nucleotidyltransferases                     |
| Drugs          | Aminoglycosides                                     | Aminoglycoside\_O-phosphotransferases                         |
| Drugs          | Bacitracin                                          | Bacitracin\_ABC\_transporter                                  |
| Drugs          | Bacitracin                                          | Undecaprenyl\_pyrophosphate\_phosphatase                      |
| Drugs          | betalactams                                         | Class\_A\_betalactamases                                      |
| Drugs          | betalactams                                         | Class\_B\_betalactamases                                      |
| Drugs          | betalactams                                         | Class\_C\_betalactamases                                      |
| Drugs          | betalactams                                         | Class\_D\_betalactamases                                      |
| Drugs          | betalactams                                         | Mutant\_porin\_proteins                                       |
| Drugs          | betalactams                                         | Penicillin\_binding\_protein                                  |
| Drugs          | betalactams                                         | Penicillin\_binding\_protein\_regulator                       |
| Drugs          | Cationic\_antimicrobial\_peptides                   | Cationic\_peptide-resistant\_16S\_ribosomal\_subunit\_protein |
| Drugs          | Cationic\_antimicrobial\_peptides                   | Defensin-resistant\_mprF                                      |
| Drugs          | Cationic\_antimicrobial\_peptides                   | Lipid\_A\_modification                                        |
| Drugs          | Cationic\_antimicrobial\_peptides                   | Polymyxin\_B\_resistance\_regulator                           |
| Drugs          | Elfamycins                                          | EF-Tu\_inhibition                                             |
| Drugs          | Elfamycins                                          | Elfamycin\_efflux\_pumps                                      |
| Drugs          | Fluoroquinolones                                    | Fluoroquinolone-resistant\_DNA\_topoisomerases                |
| Drugs          | Fluoroquinolones                                    | Fluoroquinolone\_ABC\_efflux\_pump                            |
| Drugs          | Fluoroquinolones                                    | Quinolone\_active\_efflux                                     |
| Drugs          | Fluoroquinolones                                    | Quinolone\_resistance\_protein\_Qnr                           |
| Drugs          | Fosfomycin                                          | Fosfomycin\_MFS\_efflux\_pump                                 |
| Drugs          | Fosfomycin                                          | Fosfomycin\_phosphorylation                                   |
| Drugs          | Fosfomycin                                          | Fosfomycin\_target\_mutation                                  |
| Drugs          | Fosfomycin                                          | Fosfomycin\_thiol\_transferases                               |
| Drugs          | Fusidic\_acid                                       | Fusidic\_acid-resistant\_mutation                             |
| Drugs          | Glycopeptides                                       | Bleomycin\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanA-type\_accessory\_protein                                 |
| Drugs          | Glycopeptides                                       | VanA-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanA-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanB-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanB-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanC-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanC-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | Vancomycin-resistant\_mutation                                |
| Drugs          | Glycopeptides                                       | VanD-type\_accessory\_protein                                 |
| Drugs          | Glycopeptides                                       | VanD-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanD-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanE-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanE-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanF-type\_accessory\_protein                                 |
| Drugs          | Glycopeptides                                       | VanG-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanG-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanI-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanL-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanM-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanM-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanN-type\_resistance\_protein                                |
| Drugs          | Glycopeptides                                       | VanO-type\_regulator                                          |
| Drugs          | Glycopeptides                                       | VanO-type\_resistance\_protein                                |
| Drugs          | Lipopeptides                                        | Colistin-resistant\_mutant                                    |
| Drugs          | Lipopeptides                                        | Colistin\_phosphoethanolamine\_transferase                    |
| Drugs          | Lipopeptides                                        | Daptomycin-resistant\_beta-subunit\_of\_RNA\_polymerase\_RpoB |
| Drugs          | Lipopeptides                                        | Daptomycin-resistant\_beta-subunit\_of\_RNA\_polymerase\_RpoC |
| Drugs          | Lipopeptides                                        | Daptomycin-resistant\_mutant                                  |
| Drugs          | Lipopeptides                                        | Lysocin-resistant\_mutant                                     |
| Drugs          | Metronidazole                                       | nim\_nitroimidazole\_reductase                                |
| Drugs          | MLS                                                 | 23S\_rRNA\_methyltransferases                                 |
| Drugs          | MLS                                                 | Lincosamide\_nucleotidyltransferases                          |
| Drugs          | MLS                                                 | Macrolide-resistant\_23S\_rRNA\_mutation                      |
| Drugs          | MLS                                                 | Macrolide\_esterases                                          |
| Drugs          | MLS                                                 | Macrolide\_glycosyltransferases                               |
| Drugs          | MLS                                                 | Macrolide\_phosphotransferases                                |
| Drugs          | MLS                                                 | MLS\_resistance\_ABC\_efflux\_pumps                           |
| Drugs          | MLS                                                 | MLS\_resistance\_MFS\_efflux\_pumps                           |
| Drugs          | MLS                                                 | Streptogramin\_A\_O-acetyltransferase                         |
| Drugs          | MLS                                                 | Streptogramin\_B\_ester\_bond\_cleavage                       |
| Drugs          | Multi-drug\_resistance                              | MDR\_23S\_ribosomal\_RNA\_methyltransferase                   |
| Drugs          | Multi-drug\_resistance                              | MDR\_23S\_rRNA\_mutation                                      |
| Drugs          | Multi-drug\_resistance                              | MDR\_mutant\_porin\_proteins                                  |
| Drugs          | Multi-drug\_resistance                              | MDR\_regulator                                                |
| Drugs          | Multi-drug\_resistance                              | Multi-drug\_ABC\_efflux\_pumps                                |
| Drugs          | Multi-drug\_resistance                              | Multi-drug\_MFS\_efflux\_pumps                                |
| Drugs          | Multi-drug\_resistance                              | Multi-drug\_RND\_efflux\_pumps                                |
| Drugs          | Multi-drug\_resistance                              | Multi-drug\_RND\_efflux\_regulator                            |
| Drugs          | Mupirocin                                           | Mupirocin-resistant\_isoleucyl-tRNA\_synthetase               |
| Drugs          | Mycobacterium\_tuberculosis-specific\_Drug          | Ethambutol-resistant\_mutant                                  |
| Drugs          | Mycobacterium\_tuberculosis-specific\_Drug          | Ethionamide-resistant\_mutant                                 |
| Drugs          | Mycobacterium\_tuberculosis-specific\_Drug          | Isoniazid-resistant\_mutant                                   |
| Drugs          | Mycobacterium\_tuberculosis-specific\_Drug          | Para-aminosalicylic\_acid\_resistant\_mutant                  |
| Drugs          | Mycobacterium\_tuberculosis-specific\_Drug          | Pyrazinamide-resistant\_mutant                                |
| Drugs          | Nucleosides                                         | Streptothricin\_acetyltransferase                             |
| Drugs          | Oxazolidinone                                       | Oxazolidinone-resistant\_23S\_rRNA\_mutation                  |
| Drugs          | Pactamycin                                          | Pactamycin-resistant\_16S\_ribosomal\_subunit\_protein        |
| Drugs          | Phenicol                                            | Chloramphenicol\_acetyltransferases                           |
| Drugs          | Phenicol                                            | Chloramphenicol\_hydrolase                                    |
| Drugs          | Phenicol                                            | Chloramphenicol\_phosphotransferase                           |
| Drugs          | Phenicol                                            | Phenicol-resistant\_23S\_rRNA\_mutation                       |
| Drugs          | Phenicol                                            | Phenicol\_resistance\_MFS\_efflux\_pumps                      |
| Drugs          | Pleuromutilin                                       | Pleuromutilin-resistant\_23S\_rRNA\_mutation                  |
| Drugs          | Rifampin                                            | Monooxygenase                                                 |
| Drugs          | Rifampin                                            | Rifampin-resistant\_beta-subunit\_of\_RNA\_polymerase\_RpoB   |
| Drugs          | Rifampin                                            | Rifampin\_ADP-ribosyltransferase\_Arr                         |
| Drugs          | Rifampin                                            | Rifampin\_phosphotransferase                                  |
| Drugs          | Rifampin                                            | RNA-polymerase\_binding\_protein                              |
| Drugs          | Sulfonamides                                        | Sulfonamide-resistant\_dihydropteroate\_synthases             |
| Drugs          | Tetracenomycin                                      | Tetracenomycin\_MFS\_efflux\_pump                             |
| Drugs          | Tetracyclines                                       | Tetracycline-resistant\_16S\_ribosomal\_subunit\_protein      |
| Drugs          | Tetracyclines                                       | Tetracycline\_inactivation\_enzymes                           |
| Drugs          | Tetracyclines                                       | Tetracycline\_resistance\_ABC\_efflux\_pumps                  |
| Drugs          | Tetracyclines                                       | Tetracycline\_resistance\_MFS\_efflux\_pumps                  |
| Drugs          | Tetracyclines                                       | Tetracycline\_resistance\_MFS\_efflux\_regulator              |
| Drugs          | Tetracyclines                                       | Tetracycline\_resistance\_ribosomal\_protection\_proteins     |
| Drugs          | Tetracyclines                                       | Tetracycline\_transcriptional\_repressor                      |
| Drugs          | Trimethoprim                                        | Dihydrofolate\_reductase                                      |
| Metals         | Aluminum\_resistance                                | Aluminum\_ATPase\_                                            |
| Metals         | Arsenic\_resistance                                 | Arsenic\_resistance\_protein                                  |
| Metals         | Arsenic\_resistance                                 | Arsenic\_resistance\_regulator                                |
| Metals         | Arsenic\_resistance                                 | Arsenite\_oxidase\_regulator                                  |
| Metals         | Cadmium\_resistance                                 | Cadmium\_resistance\_regulator                                |
| Metals         | Chromium\_resistance                                | Chromium\_resistance\_protein                                 |
| Metals         | Cobalt\_resistance                                  | Cobalt\_transporting\_ATPase                                  |
| Metals         | Copper\_resistance                                  | Copper\_resistance\_protein                                   |
| Metals         | Copper\_resistance                                  | Copper\_resistance\_regulator                                 |
| Metals         | Iron\_resistance                                    | Iron\_resistance\_protein                                     |
| Metals         | Lead\_resistance                                    | Lead\_P-type\_ATPase\_transporter                             |
| Metals         | Mercury\_resistance                                 | Mercury\_resistance\_protein                                  |
| Metals         | Mercury\_resistance                                 | Mercury\_resistance\_regulator                                |
| Metals         | Multi-metal\_resistance                             | Multi-metal\_ABC\_efflux\_pumps                               |
| Metals         | Multi-metal\_resistance                             | Multi-metal\_resistance\_protein                              |
| Metals         | Multi-metal\_resistance                             | Multi-metal\_resistance\_regulator                            |
| Metals         | Multi-metal\_resistance                             | Multi-metal\_RND\_efflux\_pumps                               |
| Metals         | Multi-metal\_resistance                             | Multi-metal\_RND\_efflux\_regulator                           |
| Metals         | Nickel\_resistance                                  | Nickel\_ABC\_efflux\_pumps                                    |
| Metals         | Nickel\_resistance                                  | Nickel\_ABC\_efflux\_regulator                                |
| Metals         | Nickel\_resistance                                  | Nickel\_MFS\_efflux\_pumps                                    |
| Metals         | Nickel\_resistance                                  | Nickel\_resistance\_regulator                                 |
| Metals         | Sodium\_resistance                                  | Sodium\_resistance\_protein                                   |
| Metals         | Tellurium\_resistance                               | Tellurium\_resistance\_protein                                |
| Metals         | Zinc\_resistance                                    | Zinc\_resistance\_protein                                     |
| Metals         | Zinc\_resistance                                    | Zinc\_resistance\_regulator                                   |
| Multi-compound | Biocide\_and\_metal\_resistance                     | Biocide\_and\_metal\_resistance\_protein                      |
| Multi-compound | Biocide\_and\_metal\_resistance                     | Biocide\_and\_metal\_resistance\_regulator                    |
| Multi-compound | Drug\_and\_biocide\_and\_metal\_resistance          | Drug\_and\_biocide\_and\_metal\_MFS\_efflux\_pumps            |
| Multi-compound | Drug\_and\_biocide\_and\_metal\_resistance          | Drug\_and\_biocide\_and\_metal\_resistance\_regulator         |
| Multi-compound | Drug\_and\_biocide\_and\_metal\_resistance          | Drug\_and\_biocide\_and\_metal\_RND\_efflux\_pumps            |
| Multi-compound | Drug\_and\_biocide\_and\_metal\_resistance          | Drug\_and\_biocide\_and\_metal\_RND\_efflux\_regulator        |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_ABC\_efflux\_pumps                        |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_ABC\_efflux\_regulator                    |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_MATE\_efflux\_pumps                       |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_MFS\_efflux\_pumps                        |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_MFS\_efflux\_regulator                    |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_RND\_efflux\_pumps                        |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_RND\_efflux\_regulator                    |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_SMR\_efflux\_pumps                        |
| Multi-compound | Drug\_and\_biocide\_resistance                      | Drug\_and\_biocide\_SMR\_efflux\_regulator                    |

``` r
# data is not normally distributed and is very zero-skewed
hist(counts$relabun, main = "Histogram of Relative Abundance",
     xlab = "Relative abundance", ylab = "Frequency")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

### Pairwise comparisons between farms

Wilcox Rank-Sum test for non-parametric pairwise test of relative abundnace between farms

**Genes**

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
```

There are 85 genes that are significantly different between farms, which is difficult to visualize.

This plot shows all genes between farms; most decrease from farm 1 to farm 2.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = gene)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different genes by farm (n = 85)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

**Proteins**

``` r
outsig <- data.frame()
for(i in 1:length(proteins)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, protein == proteins[i]),
                   exact = FALSE)
  
  outsig[i, "protein"] <- proteins[i]
  outsig[i, "pval"] <- round(w$p.value, 3)
  
  if(w$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

# get only significant variables and add median + IQR
sigs <- outsig %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "protein") %>% 
  group_by(farm, broadclass, type, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))
```

There are 41 significantly different resistance proteins; the plot shows that most are only present in one of the two farms.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = protein)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by farm (n = 41)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

**Resistance Type**

``` r
outsig <- data.frame()
for(i in 1:length(types)) {
  
  w <- wilcox.test(relabun ~ farm, data = filter(counts, type == types[i]),
                   exact = FALSE)
  
  outsig[i, "type"] <- types[i]
  outsig[i, "pval"] <- round(w$p.value, 3)
  
  if(w$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

# get only significant variables and add median + IQR
sigs <- outsig %>% 
  filter(sig == "TRUE") %>% 
  left_join(counts, by = "type") %>% 
  group_by(farm, broadclass, type) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(farm = factor(farm))
```

There are 18 significantly different resistance types; again, the plot shows that most are only present in one of the two farms.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = farm, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Farm", y = "Median Relative Abundance") +
  ggtitle("Significantly different types by farm (n = 18)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

### Comparisons between body sites

Kruskal-Wallis test for non-parametric comparisons between the 3 body sites on both farms.

**Genes**

``` r
# create empty dataframe to fill within the loop
outsig <- data.frame()
# iterate through all 666 genes
for(i in 1:length(patterns)) {
  
  # perform Kruskal-Wallis test
  k <- kruskal.test(relabun ~ site, data = filter(counts, pattern == patterns[i]))
  
  # collect output variables
  outsig[i, "pattern"] <- patterns[i]
  outsig[i, "pval"] <- round(k$p.value, 3)
  
  if(k$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

# get genes that had a significant difference
sigs <- outsig %>% filter(sig == "TRUE") %>%  # 360 of 666 are difference between site
  left_join(counts, by = "pattern") %>% 
  group_by(site, pattern, broadclass, type, protein, gene) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site))
```

There are 360 significant genes between body sites.

This plot visualizes the differences; most are higher in feces, then rumen, then SNP. However, one gene is higher in SNP.

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = gene)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different genes by site (n = 360)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

**Proteins**

``` r
outsig <- data.frame()
for(i in 1:length(proteins)) {
  
  k <- kruskal.test(relabun ~ site, data = filter(counts, protein == proteins[i]))
  
  outsig[i, "protein"] <- proteins[i]
  outsig[i, "pval"] <- round(k$p.value, 3)
  
  if(k$p.value < 0.05) {
    outsig[i, "sig"] <- TRUE
  } else {outsig[i, "sig"] <- FALSE}
  
}

# get significant variables
sigs <- outsig %>% filter(sig == "TRUE") %>%  
  left_join(counts, by = "protein") %>% 
  group_by(site, broadclass, type, protein) %>% 
  summarize(median = median(relabun),
            IQR = IQR(relabun)) %>% 
  ungroup() %>% 
  mutate(site = factor(site))
```

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = protein)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by site (n = 118)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

**Resistance Type**

``` r
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
```

``` r
# visualize the major trends
ggplot(data = sigs, aes(x = site, y = median, group = type)) +
  geom_point() +
  geom_jitter() +
  geom_line() +
  theme_bw() +
  labs(x = "Body Site", y = "Median Relative Abundance") +
  ggtitle("Significantly different proteins by site (n = 46)")
```

<img src="statAnalysis_files/figure-markdown_github/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />
