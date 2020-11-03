library(nlme);library(MCMCglmm);library(brms)
library(parallel); library(tidyverse); library(tidybayes);library(ggthemes);
library(lubridate)
options(scipen = 999)

cond <- read.csv("condition.csv")
str(cond)

##add phase and month, and scale body condition residuals
cond <- cond %>% mutate(
  "phase" = case_when(year==2015 ~ "peak",
                      year==2016 ~ "peak",
                      year==2017 ~ "decline",
                      year==2018 ~ "decline"),
  "month" = month(as.POSIXlt(cond$Date, format = "%d-%b-%y")),
  "scale_bc" = scale(bc))

#make a bunch of factors
cond$year <- as.numeric(cond$year)
cond$phase <- as.factor(cond$phase)
cond$hare <- as.factor(cond$hare)
cond$month <- as.factor(cond$month)
str(cond)

###fit basic bayesian model using brms


m1_brm <- brm(scale_bc ~ month +
                (1 | hare) + (1 | year/month),
              data = cond,
              warmup = 500,
              iter = 3000,
              thin=2,
              chains = 2,
              inits = "random",
              cores = parallel:::detectCores(),
              seed = 12345,
              control = list(adapt_delta = 0.85)) ##this helps convergence
m1_brm <- add_criterion(m1_brm, "waic")

save(m1_brm,file = "m1_brm_cond.rds")

summary(m1_brm)
plot(m1_brm)
pairs(m1_brm)

#extract posterior modes
colnames(posterior_samples(m1_brm))[1:8]

var.animal_id <- posterior_samples(m1_brm)$"sd_hare__Intercept"^2
var.year <- posterior_samples(m1_brm)$"sd_year__Intercept"^2
var.year.month <- posterior_samples(m1_brm)$"sd_year:month__Intercept"^2
var.res <- posterior_samples(m1_brm)$"sigma"^2

#repeatability
RDist <- var.animal_id / (var.animal_id + var.year.month + var.year + var.res)
mean(RDist);HPDinterval(as.mcmc(RDist),0.95)

RYear <- var.year / (var.animal_id + var.year.month + var.year + var.res)
mean(RYear)

RYearMonth <- var.year.month / (var.animal_id + var.year.month + var.year + var.res)
mean(RYearMonth)

RRes <- var.res / (var.animal_id + var.year.month + var.year + var.res)
mean(RRes)

# Similar, and as explained in the frequentist section, we can calculate $CV_i$ as:
CVi <- sqrt(var.animal_id) / mean(cond$scale_bc)
CV.cond <- mean(CVi);HPDinterval(as.mcmc(CVi),0.95)
CV.cond*100
