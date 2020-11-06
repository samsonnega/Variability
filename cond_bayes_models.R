library(nlme);library(MCMCglmm);library(brms)
library(parallel); library(tidyverse); library(tidybayes);library(ggthemes);
library(lubridate)
options(scipen = 999)

cond <- read.csv("condition_males.csv")
str(cond)

##add phase and month, and scale body condition residuals
cond <- cond %>% mutate(
  "phase" = case_when(year==2015 ~ "peak",
                      year==2016 ~ "peak",
                      year==2017 ~ "decline",
                      year==2018 ~ "decline"),
  "scale_bc" = scale(residuals))

# just keep what you need

cond2 <- cond %>%
  filter(!is.na('scale_bc'))  %>%
  select(id ='hare', 
         bc = 'scale_bc',
         month = 'month', 
         year = 'year',
         phase = 'phase')

#make a bunch of factors
cond2$year <- as.numeric(cond2$year)
cond2$phase <- as.factor(cond2$phase)
cond2$id <- as.factor(cond2$id)
cond2$month <- as.factor(cond2$month)
str(cond2)

###fit basic bayesian model using brms

m1_brm <- brm(bc ~ month +
                (1 | id) + (1 | year/month),
              data = cond2,
              warmup = 1000,
              iter = 3000,
              thin=2,
              chains = 2,
              inits = "random",
              cores = parallel:::detectCores(),
              seed = 12345,
              control = list(adapt_delta = 0.99)) ##this helps convergence
m1_brm <- add_criterion(m1_brm, "waic")

save(m1_brm,file = "m1_brm_cond.rds")

load("m1_brm_cond.rds")
summary(m1_brm)
plot(m1_brm)
pp_check(m1_brm)
pairs(m1_brm)

#extract posterior modes
colnames(posterior_samples(m1_brm))[1:8]

var.animal_id <- posterior_samples(m1_brm)$"sd_id__Intercept"^2
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
CVi <- sqrt(var.animal_id) / mean(cond2$bc)
CV.cond <- mean(CVi);HPDinterval(as.mcmc(CVi),0.95)


mean(var.animal_id)


### examine levels of variability----
# Vh & Vw ??? by environments
model.brms=bf(bc ~ phase + month + (0+phase||id), sigma ~ 0+phase)
fit_model.brms <- brm(model.brms, data = cond2,
                      cores    = parallel:::detectCores(),
                      refresh = 0)
summary(fit_model.brms) 
plot(fit_model.brms)
VarCorr(fit_model.brms)

#make sure priors work
prior_summary(fit_model.brms)

colnames(posterior_samples(fit_model.brms))[1:10]

#extract variance components 
var.brms <- posterior_samples(fit_model.brms, 
                              c("sd_hare__phasepeak", # among unit sd E1
                                "sd_hare__phasedecline", # among unit sd E2
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
cond.bayes <- var.brms %>% filter(.,ind=="delta.Vh" | 
                                     ind=="delta.Vw" | 
                                     ind=="delta.tau") %>% 
  ggplot(., aes(x=values, y=ind, fill=ind)) +
  geom_halfeyeh() + geom_vline(xintercept=0,linetype="dashed") +
  ylab("") + xlab("Standardized Posterior Mode") +
  scale_fill_wsj() +
  theme_base() + theme(legend.position = "none")

cond.bayes
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


###stacked bar graph code

#### Figure 1B: Variances per population ####
Vh.P <- var.brms.table$mode[1]
Vh.D <- var.brms.table$mode[2]
Vw.P <- var.brms.table$mode[3]
Vw.D <- var.brms.table$mode[4]

data.vars=data.frame(V=c(Vh.P,Vh.D,Vw.P,Vw.D),
                     level=factor(c("Vi","Vi","Vw","Vw")),
                     env=factor(c("Peak","Decline","Peak","Decline"), 
                                levels = c("Peak","Decline")))
Fig1B=data.vars %>% 
  ggplot(., aes(y=V, x=env, fill=level)) + 
  geom_bar(position="stack", stat="identity", alpha=.8, colour="black") +
  scale_fill_grey(start = 0.1, end = 0.8) + 
  theme_base() +
  theme(legend.title=element_blank(),legend.position=c(0.1, .9),
        plot.background=element_blank()) +
  ylab("Total variation") + xlab("")
Fig1B
