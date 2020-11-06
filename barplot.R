props <- read.csv("varprops.csv")
str(props)
barplot(props)
library(ggplot2)
library(tidyverse)
props <- props %>% 
  mutate(trait = factor(trait, levels = c("Dex","AUC","FCM", "Body condition","Foraging", "Hourly distance")))

  varcomps <- ggplot() + geom_bar(aes(y=prop, x = trait, fill =varcomp),
                    data = props,
                    stat = "identity") +
  scale_fill_manual(values = c("burlywood3", "burlywood4"))+
  ylab("Ratio of Variance") + xlab("") + labs(fill = "Effects")

  ggsave(varcomps, 
         file = "variance partitioning.jpg", 
         width = 10, height = 6, unit = "in", dpi = 300)
  