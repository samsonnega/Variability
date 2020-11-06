library(nlme);library(MCMCglmm);library(brms)
library(parallel); library(tidyverse); library(tidybayes);library(ggthemes)
library(broom)
options(scipen = 999)

cort <- read.csv("cortAUC.csv")
str(cort)

cort <- cort %>% mutate(
  "phase" = case_when(year==2015 ~ "peak",
                      year==2016 ~ "peak",
                      year==2017 ~ "decline",
                      year==2018 ~ "decline"))

cort$year <- as.numeric(cort$year)
cort$phase <- as.factor(cort$phase)
cort$hare <- as.factor(cort$hare)
cort$season <- as.factor(cort$season)

###fit basic bayesian model using brms
prior_summary(m1_brm)
prior.1=c(set_prior("normal(0,1)", class = "sd",),
          set_prior("normal(0,1)", class = "b", dpar="sigma"))
prior.2=c(set_prior("normal(0,10)", class = "sd"),
          set_prior("normal(0,10)", class = "b", dpar="sigma")) 

m1_brm <- (brm(log(AUC) ~ season +
                (1 | hare) + (1 | phase/year),
              data = cort,
              warmup = 1000,
              iter = 3000,
              thin=2,
              chains = 2,
              prior = prior.1,
              inits = "random",
              cores = parallel:::detectCores(),
              seed = 12345,
              control = list(adapt_delta = 0.999))) ##this helps convergence
           
m1_brm <- add_criterion(m1_brm, "waic")
##model checks
get_prior(AUC ~ season + (1 | hare) + (1 | phase/year), data = cort)
pp_check(m1_brm)
conditional_effects(m1_brm)
save(m1_brm,file = "m1_brm_cort.rds")
summary(m1_brm)
plot(m1_brm)
launch_shinystan(m1_brm)

# extract variance components
colnames(posterior_samples(m1_brm))[1:8]
var.animal_id <- posterior_samples(m1_brm)$"sd_hare__Intercept"^2
var.year <- posterior_samples(m1_brm)$"sd_phase__Intercept"^2
var.year.season <- posterior_samples(m1_brm)$"sd_phase:year__Intercept"^2
var.res <- posterior_samples(m1_brm)$"sigma"^2

#repeatability
RDist <- var.animal_id / (var.animal_id + var.year.season + var.year + var.res)
mean(RDist);HPDinterval(as.mcmc(RDist),0.95)

RYear <- var.year / (var.animal_id + var.year.season + var.year + var.res)
mean(RYear)

RYearSeason <- var.year.season / (var.animal_id + var.year.season + var.year + var.res)
mean(RYearSeason)

RRes <- var.res / (var.animal_id + var.year.season + var.year + var.res)
mean(RRes)

# Similar, and as explained in the frequentist section, we can calculate $CV_i$ as:
CVi <- sqrt(var.animal_id) / mean(cort$AUC)
CV.cort <- mean(CVi);HPDinterval(as.mcmc(CVi),0.95)
CV.cort*100
  