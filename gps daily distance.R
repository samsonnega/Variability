data<-read.csv("gps2015-18.csv")
data$date2<-as.POSIXct(strptime(as.character(data$date2),"%d/%m/%Y %H:%M", tz="GMT")) 

#trajectories 
library(adehabitatLT)
ltraj<-as.ltraj(xy=data[,2:3], date=data$date2, id=as.character(data$id))
hare<-ld(ltraj) #gets you the data for each point along the trajectory (i.e., distances)
hist(ltraj[1])
is.regular(ltraj)
plot(ltraj)
plotltr(ltraj, which = "dist")
ltraj
for(i in 1:length(ltraj)){print(median(ltraj[[i]]$dt, na.rm = TRUE))}


refda <- min(data$date2)
hare_NA <- setNA(ltraj,refda,1,units="hour")
ltraj <- sett0(hare_NA,refda,1,units="hour")
is.regular(ltraj)

#need to regularize the ltraj object making all 1 hour interval
for(i in 1:length(ltraj)){
  ref <- min(data$date2), "hours")
  ltraj[i] %>% 
    setNA(ltraj = ., date.ref = ref, dt = 1, units = "hour") %>%
    sett0(ltraj = ., date.ref = ref, dt = 1, units = "hour") -> ltraj[i]
}

#turn angle correlation
TAC <- matrix(ncol=1, nrow=length(ltraj)) # create empty data frame to populate with for-loop

for (i in 1:length(ltraj)){
  SA <- adehabitatLT::acfang.ltraj(ltraj[i], which = "relative", plot = FALSE) 
  TAC[i,] <- 1/(SA[[1]][1,])
}

TAC


###Using package amt -------
library(lubridate)
library(tidyverse)
library(ggplot2)
devtools::install_github("jmsigner/amt")
library(amt)
data<-read.csv("gps2015-18.csv")

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
 make_track(x_, y_, t_, id = id)
#add step lengths and also month
trk <- trk %>% mutate(sl_ = step_lengths(.),
                      month = lubridate::month(t_))

summary(trk$sl_)
# take a look at the mean sampling rate
summarize_sampling_rate(trk)
# nest the data by individual id
trk1 <- trk %>% 
  nest(data = -"id")

trk1 %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
   select(id, sr) %>% unnest

# resample within a range determined by the above summary
trk2 <- trk1 %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = minutes(45), tolerance = minutes(30)) %>% steps_by_burst()))

trk3 <- trk2 %>% select(id, steps) %>% unnest(cols = "steps")


trk3 <- trk2 %>% unnest(trk2, keep_empty = TRUE)
unnest(trk3, cols = c(id,month))
hoist(.data = trk2, .col = "steps")

sample <- trk1 %>% unnest

unnest_auto(sample, col ="ïd")
  ggplot(aes(sl_, fill = factor(id))) + geom_density(alpha = 0.4)

