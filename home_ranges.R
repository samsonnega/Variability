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

#first fix date
data$date2<-as.POSIXct(strptime(as.character(data$date2),"%d/%m/%Y %H:%M"))

#Take what matters and add year
data2 <- data %>%
  filter(!is.na('y'))  %>%
  mutate(year = lubridate::year(date2),
         month = lubridate::month(date2),
         yday = lubridate::yday(date2)) %>% 
  select(x_ = 'x', y_ = 'y',
         t_ = 'date2', id = 'id',
         p = 'phase', yr = 'year',
         m = 'month', yd = 'yday')

#make sure the dataset is complete
ind <- complete.cases(data2[, c("x_", "y_", "t_")])
table(ind)
#check for duplicates of id, and coordinates
ind2 <- data2 %>% 
  select(t_, x_, y_, id) %>%
  duplicated
sum(ind2)

#delete duplicates
data2 <- data2 %>% filter(!ind2)

#...build the track ----
trk <- data2 %>%
  make_track(x_, y_, t_, id = id, yr = yr, p = p, crs = CRS("+init=epsg:32607"))

has_crs(trk)
get_crs(trk)

#add step lengths and also month
trk <- trk %>% mutate(sl_ = step_lengths(.),
                      month = lubridate::month(t_),
                      year = lubridate::year(t_))

#filter for a miniumum number of bursts

filter_min_n_burst(trk, min_n = 5)

summary(trk$sl_)
# take a look at the mean sampling rate
summarize_sampling_rate(trk)

# resample within a range determined by the above summary
trk2 <- trk %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = minutes(45), tolerance = minutes(30)) %>% 
      steps_by_burst()))


#take a look at just one hare
hare.1 <- trk %>% filter(id == "30291")

trk <- remove_capture_effect(trk$id == "30291", start = days(3), end = days(3))

#look at mcp and kde home ranges
mcp1 <- hr_mcp(hare.1, levels = c(0.5, 0.95))
kde1 <- hr_kde(hare.1, levels = c(0.5, 0.95))

#plot the home ranges overlaid
plot(kde1)
plot(mcp1, add.relocations = FALSE, add = TRUE, border = "red")

# now nest data by id to apply hr to all individuals------------

dat1 <- trk %>%
  nest(data = c(x_, y_, t_))

dat1 <- dat1 %>% 
  mutate(
    week = week(t_),
    month = month(t_, label=TRUE), 
    year = year(t_),
    hour = hour(t_)
  )
# now apply home range function for different estimators to all data and clean 
#relocations at the beginning and end of each track
hr1 <- dat1 %>%
  mutate(hr_kde = map(data, hr_kde))

saveRDS(hr1, "hr1.RDS")

hr2 <- hr1 %>% select(-data) %>%
  pivot_longer(hr_mcp:hr_akde_ou, names_to = "estimator",
               values_to = "hr")

hr2.area <- hr2 %>%
  mutate(hr_area = map(hr, hr_area)) %>%
  unnest(cols = hr_area)

saveRDS(hr2.area, "hr.long.RDS")

hr2.area <- readRDS("hr.long.RDS")
head(hr2.area, 2)

hr3.area <- hr2.area %>% 
  select(id = 'id',
         year = 'yr',
         phase = 'p',
         estimator = 'estimator',
         level = 'level',
         area = 'area') %>% 
  mutate(area = area / 1e4)

write.csv(hr3.area, "hr3.csv")

#plot the means
# arrange into more tidy format
hr2.area1 <- hr2.area %>%
  mutate(area = area / 1e4)
hr2.area1$est_lab <- factor(hr2.area1$estimator,
                            levels = c(
                              "hr_mcp", "hr_kde", "hr_locoh", "hr_akde_iid", "hr_akde_ou", "hr_akde_auto",
                              "hr_akde_iid_with_error", "hr_akde_ou_with_error", "hr_akde_auto_with_error"),
                            labels = c("MCP", "KDE", "LoCoH", "aKDE (iid)", "aKDE (ou)", "aKDE (select)",
                                       "aKDE (iid, e)", "aKDE (ou, e)", "aKDE (select, e)"))
ci <- hr2.area1 %>% group_by(est_lab, yr) %>%
  summarise(m = mean(area), se = sd(area) / sqrt(n()),
            me = qt(0.975, n() - 1) * se, lci = m - me, uci = m + me)

p1 <- hr2.area1 %>% ggplot(aes(est_lab, area, col = yr)) +
  geom_pointrange(aes(x = est_lab, y = m, ymin = lci, ymax = uci, col = yr), data = ci, inherit.aes = FALSE,
                  position = position_dodge2(width = 0.5)) +
  theme_light() +
  theme(axis.title.x = element_blank()) +
  labs(y = expression(paste("HRS [", ha, "]")))

p1
