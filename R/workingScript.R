## Further parsing to get data ready for statistical analysis

require(dplyr)
require(tidyr)
require(stringr)

# figure out what is happening with parsing
filespath <- "C:/Users/emily/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/emily/blasted"

files <- list.files(filespath, full.names = TRUE)

# create empty dataframe to fill with gene counts
countsDF <- data.frame()
# create empty dataframe to fill with gene length
lengthDF <- data.frame()

# read in each file and parse
for(i in 1:length(files)) {
  
  # read the file
  dat <- read.table(files[i], header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  # assign column names
  colnames(dat) <- c("readID", "megID", "pident", "slen", "qcovs", "qcovhsp", "qcovus")
  
  # parse subject length
  slens <- dat %>% 
    select(megID, slen) %>% 
    # get the minimum gene length for each megID
    group_by(megID) %>% 
    summarize(minlen = min(slen))
  
  # get subject length (gene length) for each megID
  megdat <- dat %>% 
    group_by(megID) %>% 
    summarize(idcount = n()) %>% 
    right_join(slens, by = "megID") %>% 
    # add sample name
    mutate(sample = files[i])
  
  
  # group by megID and then by pattern
  sumdat <- megdat %>% 
    # get number of megID hits with each gene length
    group_by(megID, minlen) %>% 
    summarize(count = sum(idcount)) %>% 
    ungroup() %>% 
    # get the pattern without megID
    # this is an important grouping variable since multiple megIDs map to a single gene
    mutate(pattern = factor(sapply(strsplit(megID, "MEG_[0-9]+\\|"), `[`, 2))) %>% 
    # group by pattern
    group_by(pattern) %>% 
    summarize(genecount = sum(count),
              genelen = min(minlen)) %>% 
    # add sample name and read length
    mutate(nhits = nrow(dat),
           sample = files[i])
  
  # append to output file
  countsDF <- rbind(sumdat, countsDF)
  lengthDF <- rbind(megdat, lengthDF)
  
  # print progress report
  cat("\n finished parsing", files[i], "... \n")
  
}

# get only number of hits per sample
hits <- countsDF %>% select(sample, nhits) %>% 
  mutate(name = sapply(strsplit(sapply(strsplit(countsDF$sample, "/"), `[`, 13), ".txt"), `[`, 1)) %>% 
  select(-sample) %>% 
  group_by(name) %>% 
  distinct() %>% 
  ungroup()

# make dataframe horizontal
countsDFh <- countsDF %>% 
  # get the sample name to make it shorter
  mutate(name = sapply(strsplit(sapply(strsplit(countsDF$sample, "/"), `[`, 13), ".txt"), `[`, 1)) %>% 
  select(-c(sample, nhits)) %>% 
  group_by(pattern, genelen) %>% 
  # spread horizontally 
  pivot_wider(names_from = name, values_from = genecount, values_fill = list(genecount = 0)) %>% 
  ungroup()

# normalize to sequencing depth
seq <- read.table("./data/sampleDepth.txt", header = TRUE) %>% 
  mutate(name = as.character(name))

# make horizontal version and add different normalizations
countsv <- countsDFh %>% 
  pivot_longer(cols = starts_with("G527_"), names_to = "name", values_to = "genecount") %>% 
  mutate(pscount = genecount + 1,
         logpscount = log1p(pscount),
         # normalize to gene length
         normgenelen = case_when(
           genecount == 0 ~ 0, # if there are no hits, keep 0
           genecount != 0 ~ genecount / genelen # else, divide by gene length
         ))  %>% 
  left_join(seq, by = "name") %>% 
  mutate(
    # normalize raw data
    normseqdepth = case_when(
      genecount == 0 ~ 0,
      genecount != 0 ~ genecount / length
    ),
    # normalize the gene length normalized
    normall = case_when(
      normgenelen == 0 ~ 0,
      normgenelen != 0 ~ normgenelen / length
    )
  )


counts <- countsv
# pick out control samples
controls <- paste("G527_94_NoTemplate-DNAextraction2_S94",
                  "G527_96_NoTemplate-LibraryPrep_S96",
                  "G527_95_MockCommunity_S95",
                  "G527_93_NoTemplate-DNAextraction1_S93",
                  "G527_49_2019_DNA_S49", sep = "|")

# get dataframe without controls or Alien sample
nocontrols <- countsDF %>% 
  mutate(name = sapply(strsplit(sapply(strsplit(countsDF$sample, "/"), `[`, 13), ".txt"), `[`, 1)) %>% 
  filter(!str_detect(name, controls))

# get dataframe of controls only
controlsDF <- countsDF %>% 
  mutate(name = sapply(strsplit(sapply(strsplit(countsDF$sample, "/"), `[`, 13), ".txt"), `[`, 1)) %>% 
  filter(str_detect(name, controls)) %>% 
  # remove Alien sample
  filter(!str_detect(name, "G527_49_2019_DNA_S49"))

wrangled <- nocontrols %>% 
  mutate(
    ## MEGARES 
    # resistance class
    class = sapply(strsplit(nocontrols$pattern, "\\|"), `[`, 1),
    # type of resistance
    type = sapply(strsplit(nocontrols$pattern, "\\|"), `[`, 2),
    # protein
    protein = sapply(strsplit(nocontrols$pattern, "\\|"), `[`, 3),
    # gene
    gene = sapply(strsplit(nocontrols$pattern, "\\|"), `[`, 4),
    ## SAMPLE METADATA
    # farm number
    farm = str_remove(str_extract(nocontrols$name, "P(\\d)"), "P"),
    # animal number
    animal = str_remove(str_extract(nocontrols$name, "A(\\d{1,2})"), "A"),
    # sample number
    sample = sapply(str_split(nocontrols$name, "_"), `[`, 2),
    # body site
    site = sapply(str_split(nocontrols$name, "_"), `[`, 4)
    
  )

# get a summary of all the normalizations
sum <- wrangled %>% 
  select(-c(genelen, genecount, length)) %>% 
  pivot_longer(cols = c(pscount, logpscount, normgenelen, normseqdepth,
               normall), names_to = "normtype", 
               values_to = "value")

# get all normalizations
norms <- unique(sum$normtype)

for(i in 1:length(norms)) {
  
  # make histogram
  hist(sum$value[sum$normtype == norms[i]], main = paste0("Histogram: ", norms[i]))
  
  # make boxplot
  boxplot(value ~ farm, data = filter(sum, normtype == norms[i]),
          main = paste0("Boxplot: ", norms[i]))
 
}

# does this make sense to use relative abundance?
require(vegan)

sum <- sum %>% mutate(ids = rownames(sum))
rel <- sum %>% 
  filter(normtype == "normall") %>% 
  select(-c(normtype, class, type, protein, gene, farm, animal, sample)) %>% 
  group_by(name, pattern, site) %>% 
  spread(key = name, value = value) 


all <- sum %>% 
  filter(normtype == "normall") %>%
  select(-c(normtype, class, type, protein, gene, farm, animal, sample)) %>% 
  mutate(pattern = factor(pattern),
         name = factor(name),
         site = factor(site)) %>%  
  group_by(name, pattern) %>% 
  pivot_wider(names_from = pattern, values_from = value)


all <- nocontrols %>% 
  select(pattern, name, normall) %>% 
  group_by(name) #%>% 
  pivot_wider(id_cols = name, names_from = pattern, values_from = normall)
