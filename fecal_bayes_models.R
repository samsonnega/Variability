# Load packages
library(nlme);library(MCMCglmm);library(brms)
library(parallel); library(tidyverse); library(tidybayes);library(ggthemes)
options(scipen = 999)

#load data
fecal <- read.csv("fecal.csv")
#add phase
fecal <- fecal %>% mutate(
  "phase" = case_when(Year==2015 ~ "peak",
                      Year==2016 ~ "peak",
                      Year==2017 ~ "decline",
                      Year==2018 ~ "decline"))
str(fecal)
fecal$Time <- as.numeric(fecal$Time)
fecal$Year <- as.numeric(fecal$Year)
fecal$phase <- as.factor(fecal$phase)
fecal$Hare <- as.factor(fecal$Hare)
fecal$Sex <- as.factor(fecal$Sex)
str(fecal)
#How many repeats across phase?
table(fecal$phase,fecal$Hare)
# id,day,and sex should be factor...all others numeric
#Model comparison
bf_m1=bf(log(ng.g) ~ phase + (1|Hare))
bf_m2=bf(log(ng.g) ~ phase + (0+phase||Hare))
bf_m3=bf(log(ng.g) ~ phase + (1|Hare), sigma ~ 0+phase)
bf_m4=bf(log(ng.g) ~ phase + (0+phase||Hare), sigma ~ 0+phase)


fit_m1 <- brm(bf_m1, data = fecal,
              cores = parallel:::detectCores(),
              refresh=0,
              save_pars = save_pars(group = T))
fit_m2 <- brm(bf_m2, data = fecal,
              cores = parallel:::detectCores(),
              refresh=0,
              save_pars = save_pars(group = T))
fit_m3 <- brm(bf_m3, data = fecal,
              cores = parallel:::detectCores(),
              refresh=0,
              save_pars = save_pars(group = T))
fit_m4 <- brm(bf_m4, data = fecal,
              cores = parallel:::detectCores(),
              refresh=0,
              save_pars = save_pars(group = T))

#Create WAIC table
WAIC_fits=WAIC(fit_m1,fit_m2,fit_m3,fit_m4)
WAIC_fits

#Or use LOO criterion in Stan
LOO_fits=LOO(fit_m1,fit_m2,fit_m3,fit_m4)
LOO_fits

### Examining levels of variability within and between ----

# Vh & Vw ??? by environments
model.brms=bf(log(ng.g) ~ phase + Sex + Time + (0+phase||Hare), sigma ~ 0+phase)
fit_model.brms <- brm(model.brms, data = fecal,
                      cores    = parallel:::detectCores(),
                      refresh = 0)
summary(fit_model.brms) 
plot(fit_model.brms)

#make sure priors work
prior_summary(fit_model.brms)

colnames(posterior_samples(fit_model.brms))[1:10]
var.Hare <- posterior_samples(fit_model.brms)$"sd_Hare__Intercept"^2


#extract variance components 
var.brms <- posterior_samples(fit_model.brms, 
                              c("sd_Hare__phasepeak", # among unit sd E1
                                "sd_Hare__phasedecline", # among unit sd E2
                                "sigma_phasepeak", # within unit sd E1, log scale
                                "sigma_phasedecline")) # within unit sd E2, log scale

colnames(var.brms)=c("Vh.P", "Vh.D", "Vw.P", "Vw.D")
var.brms=as.data.frame(var.brms)
head(var.brms)
#We need to square these values to express theme as variances and take the 
#exponent for the residual standard deviation:

var.brms$Vh.P=(var.brms$Vh.P)^2
var.brms$Vh.D=(var.brms$Vh.D)^2
var.brms$Vw.P=exp(var.brms$Vw.P)^2
var.brms$Vw.D=exp(var.brms$Vw.D)^2

var.brms$tau.P=var.brms$Vh.P/(var.brms$Vh.P+var.brms$Vw.P)
var.brms$tau.D=var.brms$Vh.D/(var.brms$Vh.D+var.brms$Vw.D)

var.brms$delta.Vh=var.brms$Vh.D-var.brms$Vh.P
var.brms$delta.Vw=var.brms$Vw.D-var.brms$Vw.P
var.brms$delta.tau=var.brms$tau.D-var.brms$tau.P

var.brms=stack(var.brms)
var.brms.table=var.brms %>% group_by(ind) %>%
  summarise(mode=posterior.mode(as.mcmc(values)),
            low.ci=HPDinterval(as.mcmc(values))[1],
            up.ci=HPDinterval(as.mcmc(values))[2]) %>% data.frame()
var.brms.table[,-1]=round(var.brms.table[,-1],2)
var.brms.table

#plot the posterior mode differences
fecal.bayes <- var.brms %>% filter(.,ind=="delta.Vh" | 
                                     ind=="delta.Vw" | 
                                     ind=="delta.tau") %>% 
  ggplot(., aes(x=values, y=ind, fill=ind)) +
  geom_halfeyeh() + geom_vline(xintercept=0,linetype="dashed") +
  ylab("") + xlab("standardized posterior mode") +
  scale_fill_wsj() +
  theme_base() + theme(legend.position = "none")

fecal.bayes
#plot just the phase difference in metrics
fecal.bayes.2 <- var.brms %>% filter(.,ind=="Vh.P" | 
                                     ind=="Vh.D") %>% 
  ggplot(., aes(x=values, y=ind, fill=ind)) +
  geom_halfeyeh() + geom_vline(xintercept=0,linetype="dashed") +
  ylab("") + xlab("standardized posterior mode") +
  scale_fill_wsj() +
  theme_base() + theme(legend.position = "none")

fecal.bayes.2

fecal.bayes.3 <- var.brms %>% filter(.,ind=="Vw.P" | 
                                       ind=="Vw.D") %>% 
  ggplot(., aes(x=values, y=ind, fill=ind)) +
  geom_halfeyeh() + geom_vline(xintercept=0,linetype="dashed") +
  ylab("") + xlab("standardized posterior mode") +
  scale_fill_wsj() +
  theme_base() + theme(legend.position = "none")

fecal.bayes.3

### Now fit a basic model to examine trait repeatability 
m1_brm <- brm(log(ng.g) ~ Time + I(Time^2) + Sex +
                (1 | Hare) + (1 | Year/Time),
              data = fecal,
              warmup = 500,
              iter = 3000,
              thin=2,
              chains = 2,
              inits = "random",
              cores = parallel:::detectCores(),
              seed = 12345)
m1_brm <- add_criterion(m1_brm, "waic")

save(m1_brm,file = "m1_brm_fecal.rds")
summary(m1_brm)
plot(m1_brm)

#Can calculate repeatability by using posterior samples

colnames(posterior_samples(m1_brm))[1:8]
var.animal_id <- posterior_samples(m1_brm)$"sd_Hare__Intercept"^2
var.year <- posterior_samples(m1_brm)$"sd_Year__Intercept"^2
var.year.month <- posterior_samples(m1_brm)$"sd_Year:Time__Intercept"^2
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
CVi <- sqrt(var.animal_id) / mean(fecal$ng.g)
mean(CVi);HPDinterval(as.mcmc(CVi),0.95)



