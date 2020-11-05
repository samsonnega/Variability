library(tidyverse)
library(MESS)
cort <- read.csv("total_cort_2015-18.csv")
cort$bleed <- as.factor(cort$bleed)
cort$hare <- as.factor(cort$hare)
cort$year <- as.factor(cort$year)
cort$cort <- as.numeric(cort$cort)

cortAUC_full <- cort %>% ## cortImputed = the dataset with imputed values replacing the missing values
  group_by(year, season, hare) %>% 
  ## grouped by all these variables so they would stay in the dataset after summarizing
  summarize(AUC = MESS::auc(bleed,cort, from = 2, to = 5, type = "spline", na.rm = TRUE)) 
## auc calculates the area under the curve using 'sample' as the key 
##(base = 1, dex = 2, 30 = 3, 60 = 4) and 'Cort3' as the variable

write.csv(cortAUC_full, "cortAUC.csv")

### compute BASE - DEX
cortDEX <- cort %>% 
  group_by(year, season, hare) %>% 
  spread(key = bleed, value = cort) %>%
  mutate(diff = `1`-`2`)

write.csv(cortDEX, "cortDEX.csv")


