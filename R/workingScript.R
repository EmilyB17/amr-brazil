## post-decontam normalization

## get the gene lengths for each gene
len <- countsDF %>% select(pattern, genelen) %>%
  mutate(pattern = as.character(pattern)) %>% 
  group_by(pattern) %>% 
  summarize(genelen = min(genelen))

# read table with sequencing depth info
seqs <- read.table("./data/sampleDepth.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>% 
  rename(names = name,
         nreads = length)

# make dataframe from OTU table
counts <- as.data.frame(nocontrol@otu_table) %>% 
  rownames_to_column(var = "pattern") %>% 
  pivot_longer(cols = starts_with("G527_"), names_to = "names", values_to = "count") %>% 
  # add gene length
  inner_join(len, by = "pattern")  %>% 
  # add sequencing depth
  inner_join(seqs, by = "names")

## NORMALIZE: to gene length and number of reads
norm <- counts %>% 
  mutate(
    # normalize to gene length
    normgenelen = 
      case_when(
        # if there are no genes, keep as 0
        count == 0 ~ 0,
        # otherwise, divide by gene length
        count != 0 ~ count / genelen
      ),
          
    # then normalize by number of reads
    normseqdepth = 
      case_when(
        # if there are no genes, keep as 0
        count == 0 ~ 0,
        # otherwise, divide by number of reads
        count != 0 ~ normgenelen / nreads
      )
  ) %>% 
  # keep only final normalization
  select(pattern, names, normseqdepth)


# calculate relative abundance

# make the dataframe horizontal to calculate relative abundance
normh <- norm %>% 
  pivot_wider(names_from = pattern, values_from = normseqdepth, values_fill = list(normseqdepth = 0)) %>% 
  # decostand() requires a matrix, so make the names into rownames
  column_to_rownames(var = "names")

# calculate relative abundance
relabun <- decostand(normh, method = "range", MARGIN = 2) %>% 
  # bring back the names column
  rownames_to_column(var = "name")

