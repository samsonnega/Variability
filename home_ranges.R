###Using package amt -------
#load in the data
data<-read.csv("gps2015-18.csv")
#load packages
library(lubridate)
library(tidyverse)
library(ggplot2)
devtools::install_github("jmsigner/amt")
library(amt)
library(sf)

#Take what matters
data2 <- data %>%
  filter(!is.na('y'))  %>%
  select(x_ = 'x', y_ = 'y',
         t_ = 'date2', id = 'id')
#make sure the dataset is complete
all(complete.cases(data2))
#check for duplicated time stamps
any(duplicated(data2$t_))
#delete duplicates
data2 <- data2[!duplicated(data2$t_), ]

#make R understand dates and time
data2$t_<-as.POSIXct(strptime(as.character(data2$t_),"%d/%m/%Y %H:%M"))

#...build the track ----
trk <- data2 %>%
  make_track(x_, y_, t_, id = id, crs = CRS("+init=epsg:32607"))

has_crs(trk)
get_crs(trk)
#add step lengths and also month
trk <- trk %>% mutate(sl_ = step_lengths(.),
                      month = lubridate::month(t_),
                      year = lubridate::year(t_))

summary(trk$sl_)
# take a look at the mean sampling rate
summarize_sampling_rate(trk)

# resample within a range determined by the above summary
trk2 <- trk1 %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = minutes(45), tolerance = minutes(30)) %>% 
      steps_by_burst()))


#take a look at just one hare
hare.1 <- trk %>% filter(id == "29405")

#look at mcp and kde home ranges
mcp1 <- hr_mcp(hare.1, levels = c(0.5, 0.95))
kde1 <- hr_kde(hare.1, levels = c(0.5, 0.95))

#plot the home ranges overlaid
plot(kde1)
plot(mcp1, add.relocations = FALSE, add = TRUE, border = "red")

# now nest data by id to apply hr to all individuals
dat1 <- trk %>% 
  nest(data = -"id")

# now apply home range function for different estimators to all data
hr1 <- dat1 %>%
  mutate(
    hr_mcp = map(trk, hr_mcp),
    hr_kde = map(trk, hr_kde),
    hr_locoh = map(trk, ~ hr_locoh(., n = ceiling(sqrt(nrow(.))))),
    hr_akde_iid = map(trk, ~ hr_akde(., fit_ctmm(., "iid"))),
    hr_akde_ou = map(trk, ~ hr_akde(., fit_ctmm(., "ou")))
  )

