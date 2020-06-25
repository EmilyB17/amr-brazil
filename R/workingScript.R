
### PLOTS FROM STATISTICAL ANALYSIS
# sample sizes are uneven; 48 farm 1, 43 farm 2-- need to use per sample average

require(RColorBrewer)

# define function for standard error
se <- function(x) {sqrt(stats::var(x)/length(x))}

## ---- ALL DATA ----

# get presence/absence
pres <- counts %>% 
  mutate(pres = case_when(
    relabun == 0 ~ 0,
    relabun != 0 ~ 1
  ))

# calculate per sample average count
tot <- pres %>% 
  group_by(farm, broadclass) %>% 
  summarize(sum = sum(pres),
            ster = se(pres),
            sd = sd(pres)) %>% 
  mutate(avg = case_when(
    farm == 1 ~ sum / 48,
    farm == 2 ~ sum / 43
  )) %>% 
  ungroup() %>% 
  mutate(Farm = factor(farm))

# plot by resistance class
p <- ggbarplot(tot, x = "broadclass", y = "avg", fill = "Farm",  position = position_dodge()) + 
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Resistance Class", y = "Average gene prevalence")

# save plot
ggexport(p, filename = "./data/plots/avg-presence-barplot-byfarm.jpg")

## ---- pie chart ----

## pie chart of gene prevalence heterogeneity
df <- pres %>% 
  group_by(pattern, broadclass) %>% 
  summarize(tot = sum(pres)) %>% 
  mutate(num =
           case_when(
             tot < 5 ~ "< 5",
             tot %in% 5:20 ~ "5 - 20",
             tot %in% 21:50 ~ "21 - 50",
             tot > 51 ~ "> 51"
           )) %>% 
  mutate(num = factor(num, ordered = TRUE,
                      levels = c("< 5", "5 - 20", "21 - 50", "> 51")))

# make chart
ggplot(data = df, aes(x = '',  fill = num)) +
  geom_bar(width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "Gene Representation")


## ---- multi biocide resistance ----

# important proteins - ABC and RND efflux pumps
prots <- paste("Multi-biocide_ABC_efflux_pump", "Multi-biocide_RND_efflux_pump", sep = "|")

tot <- pres %>% 
  filter(str_detect(protein, prots)) %>% 
  group_by(farm, protein) %>% 
  summarize(sum = sum(pres),
            ster = se(pres),
            sd = sd(pres)) %>% 
  mutate(avg = case_when(
    farm == 1 ~ sum / 48,
    farm == 2 ~ sum / 43
  )) %>% 
  ungroup() %>% 
  mutate(Farm = factor(farm)) %>% 
  recode(protein,
         c(Multi-biocide_ABC_efflux_pump = "ABC efflux pump",
         Multi-biocide_RND_efflux_pump = "RND efflux pump"))


ggbarplot(tot, x = "protein", y = "avg", fill = "Farm",  position = position_dodge()) + 
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Resistance Class", y = "Average gene prevalence")
