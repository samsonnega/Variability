#####################Variability metrics########################################
#####################Sam Sonnega################################################
####################OCT 2020####################################################

###Text from Hertel et al. 2020:
# As we can see the credible interval around our repeatability estimate is relatively 
# wide. A reason for this could be that the number of individuals tested is still 
# relatively limited (n = 35). 

# We can see that about 22% of the remaining variance (after controlling for month 
# and sex) in daily movement distance can be attributed to differences between 
# individuals. This means that some elephants always move over shorter distances 
# compared to other elephants and this difference is not caused by predictable
# monthly variation or sex differences. In a real study we would probably fit
# other relevant explanatory variables, like age of the animal, herd identity, 
# and habitat composition in the elephants monthly range. 
# On the other hand there does not seem to be a year effect on movement with 
# elephants moving more or less in a particular year (R.year = 0.004).

# Alternatively to repeatability we can also calculate the coefficient of 
# variation for between individual variance ($CV_i$). $CV_i$ is not confounded 
# by within-individual effects and is hence a better measure to compare the extent 
# of between-individual differences among populations and traits. It is calculated as

# $CV_i$ = $\frac{\sqrt{V_{ID}}}{\bar{x}}$

# i.e. dividing the square root of the among-individual variance by the intercept 

#packages
library(lme4);library(arm);library(MuMIn);library(tidyverse)
library(plyr);library(broom);library(coda);library(grid)
library(gridExtra);library(brms); library(broom.mixed); library(merTools);
library(tidybayes);library(parallel); library(MCMCglmm)

fecal <- read.csv("fecal.csv")
#add phase
fecal <- fecal %>% mutate(
  "phase" = case_when(Year==2015 ~ "peak",
                      Year==2016 ~ "peak",
                      Year==2017 ~ "decline",
                      Year==2018 ~ "decline"))
####BAYESIAN MODELS----


model.brms=bf(log(ng.g) ~ phase + (0+phase||Hare), sigma ~ 0+phase)
fit_model.brms <- brm(model.brms, data = fecal,
                      cores    = parallel:::detectCores(),
                      refresh = 0)
summary(fit_model.brms) 
plot(fit_model.brms)

# We should make it a habit to inspect a) the effective sample size and b) the 
# Rhat parameter to make sure that our model converged. In a nutshell, we want 
# "Eff.Sample" close to the expected sample size for all parameters. The expected 
# sample size is the number of iterations minus the number of iterations which we 
# discard as warm-up and divided by the thinning interval and multiplied by the 
# number of chains defined in the model (2 chains in our case), i.e. ((3000 - 500) / 2) * 2 = 2500
# Rhat on the other hand should be 1. 

# Similar to the frequentist approach, we can calculate repeatability by dividing 
# the posterior distribution for the variance explained by the random intercept
# animal_id by the total variance. 
# *Importantly* - in brms variance parameters are given in standard deviations 
# and need to be squared to calculate the variance! 

# We can either calculate repeatability by hand:

(2.32)^2 / ((2.32)^2 + (0.59)^2 + (1.07)^2 + (4.29)^2)

# Or (better) take the mean and credible interval of the posterior distribution.

colnames(posterior_samples(m1_brm))[1:8]

var.animal_id <- posterior_samples(m1_brm)$"sd_animal_id__Intercept"^2
var.year <- posterior_samples(m1_brm)$"sd_year__Intercept"^2
var.year.month <- posterior_samples(m1_brm)$"sd_year:month__Intercept"^2
var.res <- posterior_samples(m1_brm)$"sigma"^2

RDist <- var.animal_id / (var.animal_id + var.year.month + var.year + var.res)
mean(RDist);HPDinterval(as.mcmc(RDist),0.95)


# The mean and credible interval for this posterior distribution is 0.21 [0.14, 0.3], 
# so similar but slightly higher than our result in the frequentist approach. It means 
# that after controlling for the fixed effects of month and Sex - 21% of the 
# remaining variance can be explained by individual differences in behavioral 
# expression among elephants.

RYear <- var.year / (var.animal_id + var.year.month + var.year + var.res)
mean(RYear)

RYearMonth <- var.year.month / (var.animal_id + var.year.month + var.year + var.res)
mean(RYearMonth)

RRes <- var.res / (var.animal_id + var.year.month + var.year + var.res)
mean(RRes)

# Similar, and as explained in the frequentist section, we can calculate $CV_i$ as:

CVi <- sqrt(var.animal_id) / mean(data$meanDailyDisplacement)
mean(CVi);HPDinterval(as.mcmc(CVi),0.95)


####BRMS model comparison with WAIC----

library(nlme);library(MCMCglmm);library(brms)
library(parallel); library(tidyverse); library(tidybayes);library(ggthemes)
options(scipen = 999)


# Vh & Vw ??? by environments
model.brms=bf(log(ng.g) ~ phase + (0+phase||Hare), sigma ~ 0+phase)
fit_model.brms <- brm(model.brms, data = fecal,
                      cores    = parallel:::detectCores(),
                      refresh = 0)
summary(fit_model.brms) 
plot(fit_model.brms)

#make sure priors work
prior_summary(fit_model.brms)
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







####  MCMC method ----
library(MCMCglmm)
prior.univ.hierarch<- list(R = list(V = diag(2), n = 1.002),
                           G = list(G1 = list(V = diag(2), n = 1.002))) # heterogeneous Vh & Vw

model.mcmc<-MCMCglmm(log(ng.g)~phase,random=~idh(phase):Hare, 
                     rcov=~idh(phase):units,
                     data=fecal, verbose=F, 
                     prior=prior.univ.hierarch) 
summary(model.mcmc)

Vh.E1.mcmc<-as.numeric(model.mcmc$VCV[,1])
Vh.E2.mcmc<-as.numeric(model.mcmc$VCV[,2])
Vw.E1.mcmc<-as.numeric(model.mcmc$VCV[,3])
Vw.E2.mcmc<-as.numeric(model.mcmc$VCV[,4])
tau.E1.mcmc<-Vh.E1.mcmc/(Vh.E1.mcmc+Vw.E1.mcmc)
tau.E2.mcmc<-Vh.E2.mcmc/(Vh.E2.mcmc+Vw.E2.mcmc)
Vh.diff.mcmc<-as.numeric(posterior.mode(as.mcmc(Vh.E2.mcmc-Vh.E1.mcmc)))
Vw.diff.mcmc<-as.numeric(posterior.mode(as.mcmc(Vw.E2.mcmc-Vw.E1.mcmc)))
tau.diff.mcmc<-as.numeric(posterior.mode(as.mcmc(tau.E2.mcmc-tau.E1.mcmc)))
Vh.diff.mcmc.ci<-as.numeric(HPDinterval(as.mcmc(Vh.E2.mcmc-Vh.E1.mcmc)))
Vw.diff.mcmc.ci<-as.numeric(HPDinterval(as.mcmc(Vw.E2.mcmc-Vw.E1.mcmc)))
tau.diff.mcmc.ci<-as.numeric(HPDinterval(as.mcmc(tau.E2.mcmc-tau.E1.mcmc)))




V.diff.mcmc.table=data.frame(V.component.diff=c("Vh","Vw","tau"),delta.simulated=c(0.75,0,.3),delta.estimated=c(Vh.diff.mcmc,Vw.diff.mcmc,tau.diff.mcmc),low.ci=c(Vh.diff.mcmc.ci[1],Vh.diff.mcmc.ci[1],tau.diff.mcmc.ci[1]),up.ci=c(Vh.diff.mcmc.ci[2],Vh.diff.mcmc.ci[2],tau.diff.mcmc.ci[2]))
V.diff.mcmc.table[,-1]=round(V.diff.mcmc.table[,-1],2)
V.diff.mcmc.table


plot(model.mcmc)
#looks good, priors ok




###################GRAPHS----
library(ggthemes)
library(plotMCMC)
plotTrace(allChains)


fecal.bayes <- var.brms %>% filter(.,ind=="delta.Vh" | 
                      ind=="delta.Vw" | 
                      ind=="delta.tau") %>% 
  ggplot(., aes(x=values, y=ind, fill=ind)) +
  stat_eye() + geom_vline(xintercept=0,linetype="dashed") +
  ylab("") + xlab("standardized posterior mode") +
  theme_base() + theme(legend.position = "none")

fecal %>%
  modelr::data_grid(phase) %>%
  add_fitted_draws(fit_model.brms, dpar = c("mu", "sigma")) %>%
  sample_draws(30) %>%
  ggplot(aes(y = phase)) +
  stat_dist_slab(aes(dist = "norm", arg1 = mu, arg2 = sigma),
                 slab_color = "gray65", alpha = 1/10, fill = NA
  ) +
  geom_point(aes(x = ng.g), data = fecal, shape = 21, fill = "#9ECAE1", size = 2)

ggsave(fecal.bayes, 
       file = "fecal.bayes.jpg", 
       width = 6, height = 6, unit = "in", dpi = 300)

#Theme from Hertel et al.
posteriorBT <- posterior_samples(m1_brm)[,9:43] %>%
  gather(animal_id, value, 
         "r_animal_id[elephant1,Intercept]" : "r_animal_id[elephant9,Intercept]")%>%
  separate(animal_id, 
           c(NA,NA,NA,"animal_id",NA), 
           sep = "([\\_\\[\\,])", fill = "right") %>%
  left_join(select(data[!duplicated(data$animal_id),], animal_id, Sex))

posteriorBT[posteriorBT$Sex == "F",]$value <- 
  posteriorBT[posteriorBT$Sex == "F",]$value + fixef(m1_brm, pars = "Intercept")[1]

posteriorBT[posteriorBT$Sex == "M",]$value <- 
  posteriorBT[posteriorBT$Sex == "M",]$value + fixef(m1_brm, pars = "Intercept")[1] + 
  fixef(m1_brm, pars = "SexM")[1] 

posteriorBT$col <- ifelse(posteriorBT$animal_id %in% 
                            c("elephant17", "elephant4", "elephant8",
                              "elephant36","elephant20"),
                          posteriorBT$animal_id, "Other individuals")

posteriorBT <- posteriorBT %>%
  dplyr::group_by(animal_id) %>%
  dplyr::mutate(meanBT = mean(value))%>%
  dplyr::ungroup()

BT <- ggplot()+
  ggridges::geom_density_ridges(data = posteriorBT, 
                                aes(x = value,
                                    y = reorder(as.factor(animal_id), meanBT),
                                    height = ..density.., 
                                    fill = col,scale = 3), alpha = 0.6)+
  geom_point(data = posteriorBT[!duplicated(posteriorBT$animal_id),],
             aes(x = meanBT, 
                 y = as.factor(animal_id),
                 col = Sex),
             size = 1)+
  labs(y = "", 
       x = "BT mean daily distance (km)", 
       fill = "ID")+
  theme_classic()+
  scale_fill_manual(values = c("#F8766D","#C77CFF","#7CAE00",
                               "#FFCC00","#00BFC4","gray"))

