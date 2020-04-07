## troubleshoot 

# run parsingScript loop on one file

## what is happening with multiple lengths
tshoot <- dat %>% group_by(megID) %>% 
  summarize(rangelo = range(slen)[1],
            rangehi = range(slen)[2]) %>% 
  mutate(matched = case_when(
    rangelo == rangehi ~ "match",
    rangelo != rangehi ~ "miss"
  ))# %>% 
  filter(matched == "miss")


# there are difference sequence lengths for some database entries matching one megID
head(tshoot)

# what has a difference of over length 20
tshoot <- tshoot %>% 
  mutate(diff = rangehi - rangelo)

tshoot[tshoot$diff < 20 , ]

big <- tshoot$pattern[tshoot$diff > 2400]
grep(big, ids$pattern, fixed = TRUE)
ids[6521:6529,]
