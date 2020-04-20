
# calculate presence/absence
presab <- norm %>% 
  mutate(presence = case_when(
    normseqdepth == 0 ~ 0,
    normseqdepth != 0 ~ 1
  )) %>% 
  select(-normseqdepth)

# make horizontal by each gene
presabh <- presab %>% 
  pivot_wider(names_from = pattern, values_from = presence)

# calculate total presence for each gene
sums <- colSums(presabh %>% select(-names))


# add grouping variables
prestotal <- presab %>% 
  rename(name = names) %>% 
  mutate(
    ## MEGARES 
    # resistance class
    broadclass = sapply(str_split(pattern, "\\|"), `[`, 1),
    # type of resistance
    type = sapply(str_split(pattern, "\\|"), `[`, 2),
    # protein
    protein = sapply(str_split(pattern, "\\|"), `[`, 3),
    # gene
    gene = sapply(str_split(pattern, "\\|"), `[`, 4),
    ## SAMPLE METADATA
    # farm number
    farm = case_when(
      # if the sample is a control, name it "control"
      str_detect(name, controls) ~ "control",
      # otherwise, detect the farm number
      !str_detect(name, controls) ~ str_remove(str_extract(name, "P(\\d)"), "P")
    ),
    # animal number
    animal = case_when(
      # if the sample is a control, name it "control"
      str_detect(name, controls) ~ "control",
      # otherwise, detect the animal number
      !str_detect(name, controls) ~ str_remove(str_extract(name, "A(\\d{1,2})"), "A")
    ),
    # sample number - this does NOT matter if it's control or not
    samplenum = sapply(str_split(name, "_"), `[`, 2),
    # body site
    site = case_when(
      # if the sample is a control, name it "control"
      str_detect(name, controls) ~ "control",
      # otherwise, detect the animal number
      !str_detect(name, controls) ~ sapply(str_split(name, "_"), `[`, 4)
    )
    # close mutate() parenthesis
  )




f <- prestotal %>% 
  group_by(pattern, farm) %>% 
  summarize(count = sum(presence))
