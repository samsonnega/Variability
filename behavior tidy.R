behav <- read.csv("behavior.csv")
behav$Date<-as.POSIXct(strptime(as.character(behav$Date),"%d/%m/%Y"))

behav2 <- behav %>% 
  select(id = 'Hare',
         date = 'Date',
         forage = 'Foraging',
         hour = 'Hour') %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date),
         day = lubridate::day(date),
         foragemin = (forage/60),
         phase = case_when(year==2015 ~ "peak",
                             year==2016 ~ "peak",
                             year==2017 ~ "decline",
                             year==2018 ~ "decline"))

write.csv(behav2, "foraging.csv")
