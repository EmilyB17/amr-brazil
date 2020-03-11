
# Script to parse and merge all blast output files

dat <- read.table("./parse/G527_10_P1A26_feces_S10.txt", header = FALSE, stringsAsFactors = FALSE)

# read in annotations CSV
ids <- read.csv("./megares/megares_to_external_header_mappings_v2.00.csv",
                 stringsAsFactors = FALSE)

# get unique genes (database has duplicates)
ids$pattern <- sapply(strsplit(ids$MEGARes_v2_header, "MEG_[0-9]+\\|"), `[`, 2)

# get path to all files
PATH <- paste0(getwd(), "/parse/")

# get files from the path 
files <- list.files(PATH, full.names = TRUE)

# error check: is each file unique (no duplicates)?
if(!length(unique(files)) == length(files)) {
  
  stop("There are duplicate files")
  
}

outDF <- data.frame()
# read in each file and parse
for(i in 1:length(files)) {
  
  # read the file
  dat <- read.table(files[i], header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  # assign column names
  colnames(dat) <- c("readID", "megID", "pident", "slen", "qcovs", "qcovhsp", "qcovus")
  
  # parse subject length
  slens <- dat %>% 
    select(megID, slen) %>% 
    distinct()
  
  # get subject length (gene length) for each megID
  megdat <- dat %>% 
    group_by(megID) %>% 
    summarize(idcount = n()) %>% 
    right_join(slens, by = "megID")
  
  # get the pattern without megID
  megdat$pattern <- sapply(strsplit(megdat$megID, "MEG_[0-9]+\\|"), `[`, 2)
  dat$pattern <- sapply(strsplit(dat$megID, "MEG_[0-9]+\\|"), `[`, 2)
  
  # group by pattern
  sumdat <- dat %>% 
    # get number of gene hits
    group_by(pattern) %>% 
    summarize(count = n()) %>% 
    ungroup() 
  
  # add gene length (slen)
  lendat <- sumdat %>% 
    right_join(megdat, by = "pattern") %>% 
    select(-c(megID, idcount)) %>% 
    # add sample name and read length
    mutate(sample = files[i],
           depth = nrow(dat))
  
  # append to output file
  outDF <- rbind(lendat, outDF)
  
}

# normalize by sequencing depth 
require(stringr)
seqs <- read.table("./parse/mergestats.txt", stringsAsFactors = FALSE, header = TRUE)  %>% 
  filter(!file == "file") %>% 
  mutate(name = sapply(strsplit(sapply(strsplit(seqs$file, "/"), `[`, 7), ".extendedFrags.fastq"), `[`, 1)) %>% 
  select(name, num_seqs)

seqs$length <- as.integer(str_remove_all(seqs$num_seqs, ","))

# merge
both <- outDF %>% 
  mutate(name = sapply(strsplit(sapply(strsplit(outDF$sample, "/"), `[`, 13)), ))
