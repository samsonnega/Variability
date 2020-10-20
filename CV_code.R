### ASSESSING VARIABILITY IN FECAL CORTISOL
### Code written Jan 2020
# load packages
library(scales)
library(lme4)
library(nlme)
library(bbmle)
library(stargazer)
library(dplyr)
library(raster)
library(ggplot2)
library(rptR)
library(patchwork)

#read in data and clean---------------------------------------------------------
#set as factors
data <- read.csv("New_Fecal_2019.csv",header=TRUE)
data$Year<-as.factor(data$Year)
data$Hare<-as.factor(data$Hare)
#insert column for julian day
date <- as.POSIXlt(data$Date, format = "%d-%b-%y")
data$jday <- date$yday
#scale the cort levels
data$CORT <- rescale(data$ng.g, to = c(0,100))


### Run models comparing FCM variability ---------------------------------------
#MODEL 1: NO RANDOM SLOPES OR RANDOM INTERCEPTS
mod1<-lm(log(ng.g)~Year + scale(poly(jday,2)),
         data,na.action=na.omit)

#MODEL 2: NO RANDOM SLOPES BUT WITH RANDOM INTERCEPTS
mod2<-lmer(log(ng.g)~scale(poly(jday,2)) + Year
           +(1|Hare),
           data,na.action=na.omit,REML = F)

#MODEL 3: RANDOM SLOPES AND RANDOM INTERCEPTS BUT NO COVARIANCE BETWEEN SLOPES AND INTERCEPTS
mod3<-lmer(log(ng.g)~scale(poly(jday,2)) + Year
           +(Year|Hare),
           data,na.action=na.omit,REML = F)

#MODEL 4: RANDOM SLOPES AND RANDOM INTERCEPTS, WITH COVARIANCE BETWEEN SLOPES AND INTERCEPT
mod4<-lmer(log(ng.g)~assay_id+breed_status+scale(poly(juldate,2))+scale(loc_density)+year.f
           +(1+scale(loc_density)|squirrel_id),
           data,na.action=na.omit,REML = F)


#GENERATE AICc TABLE
bbmle::AICctab(mod1,mod2,base = T,weights = T)


#GENERATE MODEL OUTPUT AND CONFIDENCE INTERVALS VIA BOOTSTRAP FOR THE TWO TOP MODELS 
#(IE MODEL 3 AND MODEL 4)
summary(mod2)
confint.merMod(mod2,method=c("boot"))
summary(mod4)
confint.merMod(mod4,method=c("boot"))


# fit your 'null' model
# by default, variance is assumed equal between the 2 groups
fit <- lme(CORT~Year, data, random = ~1|Hare)

# update the model to include a new weights argument,
# where we will call the varIdent function.
# it lets us split the variance across reference levels.

fit1 <- update(fit, weights = varIdent(form=~1|Year))

# run the anova to test fit
mod <- anova(fit, fit1)

# fit1 has an extra degree of freedomfrom estimating a 2nd parameter!
# the model still contains significantly more information,
# because 'group 1' has sd=10 and 'group 2' has sd=2!

stargazer(mod, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          out = "table.doc")

stargazer(mod, type = "html",
          out = "table2.doc")

#p value adjustment
p.fdr <- p.adjust(p, method = "fdr")

### Same for Body Condition ----------------------------------------------------
BC <- read.csv("condition.csv")
BC$cond <- rescale(BC$residuals, to = c(0,100))
BC$cond <- as.factor(BC$cond)
fitBC <- lme(cond~Exp_Year, BC, random = ~1|Hare)

fit1BC <- update(fitBC, weights = varIdent(form=~1|Exp_Year))

modBC <- anova(fitBC,fit1BC)

stargazer(modBC, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          out = "tableBC.htm")

### Same for Accelerometer ----------------------------------------------------
behav <- read.csv("2015-18 Behaviors Cleaned.csv")
str(behav)
behav$Year <- format(as.Date(behav$Date, format="%d/%m/%Y"),"%Y")
behav$Year <- as.factor(behav$Year) 
behav$Foraging <- as.factor(behav$Foraging)
behav$Forag_std <- rescale(behav$Foraging, to = c(0,100))

fitForage <- lme(Forag_std~Year, behav, random = ~1|Hare)

fit1Forage <- update(fitForage, weights = varIdent(form=~1|Year))

modForage <- anova(fitForage,fit1Forage)

modForage

stargazer(modForage, type = "html",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "",
          out = "tableForage.doc")



### Same for hormone challenge data
total <- read.csv("total_cort_2015-18.csv")
total <- na.omit(total)
total$year <- as.factor(total$year)
total$cort <- rescale(total$cort, to = c(0,100))

outTotal <- 
  total%>%
  group_by(year,bleed)%>%
  summarise(cv= cv(cort))

outDEX <- 
  total%>%
  filter(bleed==2)%>%
  group_by(year)%>%
  summarise(cv= cv(cort))

outACTH <- 
  total%>%
  filter(bleed==3)%>%
  group_by(year)%>%
  summarise(cv= cv(cort))

outBASE <- 
  total%>%
  filter(bleed==1)%>%
  group_by(year)%>%
  summarise(cv= cv(cort))

### Summarize CVs

outCORT <- 
  data %>%
  group_by(Year)%>%
  summarise(cv= cv(CORT))

outForage <- 
  behav %>%
  group_by(Year)%>%
  summarise(cv= cv(Forag_std))

outBC <- 
BC %>%
  group_by(Exp_Year)%>%
  summarise(cv= cv(cond))

outBC$Exp_Year <- as.factor(outBC$Exp_Year)

### Examine by metric

CVfecal <- data %>%
  summarise(cv= cv(CORT))

CVForage <- 
  behav %>%
  summarise(cv= cv(Forag_std))

CVBC <- 
  BC %>%
  summarise(cv= cv(cond))

CVDex2 <- 
  total%>%
  filter(bleed == 2) %>%
  summarise(cv= cv(cort))
  
CVACTH <- 
  total%>%
  filter(bleed == 3) %>%
  summarise(cv= cv(cort))

CVBase <- 
  total%>%
  filter(bleed == 1) %>%
  summarise(cv= cv(cort))
            
#combine
AllCV <- rbind(CVcort, CVForage, CVBC, CVDex2, CVACTH, CVBase)
write.csv(AllCV, "Metric variances.csv")
#read back in after changing names in EXCEL
AllCV <- read.csv("Metric variances.csv")


###REpeatability----------------------------------------------------------------
#using ICC
ICC.cort <- ICCest(Hare, ng.g, data = data)
ICC.behav <- ICCest(Hare, forage.prop, data=behav)
ICC.cond <- ICCest(Hare, cond, data =BC)

#using rptR
behavRpT <- rpt(forage.prop ~ (1 | Hare), grname = "Hare", data = behav, 
            datatype = "Gaussian", nboot = 1000, npermut = 0)

CortRpT <- rpt(CORT ~ (1 | Hare), grname = "Hare", data = data, 
                datatype = "Gaussian", nboot = 1000, npermut = 0)

BCRpT <- rpt(cond ~ (1 | Hare), grname = "Hare", data = BC, 
            datatype = "Gaussian", nboot = 1000, npermut = 0)

Repeatabilities <- (behavRpTR[[Hare]])

### graph-----------------------------------------------------------------------

(CORT.cv.plot <- ggplot(outCORT, aes(x=Year, y=cv, fill=Year))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    ggtitle("Fecal") +
    scale_fill_grey() +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
  # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

ggsave(CORT.cv.plot, 
       file = "Variability of Cort.jpg", 
       width = 6, height = 6, unit = "in", dpi = 300)
 
(BC.cv.plot <- ggplot(outBC, aes(x=Exp_Year, y=cv, fill=Exp_Year))+  
  # Add a boxplot
  geom_bar(stat = "identity") +
  ggtitle("Condition") +
  scale_fill_grey() +
  labs(y = "Standardized CV (%)") +
  coord_flip() +
  # The colors to assign each time period
  # Black and white theme
  theme_classic() +
  # Remove the x-axis title
  theme(legend.position = "none", axis.title.y=element_blank()))

(Forage.cv.plot <- ggplot(outForage, aes(x=Year, y=cv, fill=Year))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    ggtitle("Behavior") +
    scale_fill_grey() +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
    # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

(DEX.cv.plot <- ggplot(outDEX, aes(x=year, y=cv, fill=year))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    ggtitle("DEX") +
    scale_fill_grey() +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
    # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

(ACTH.cv.plot <- ggplot(outACTH, aes(x=year, y=cv, fill=year))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    ggtitle("ACTH") +
    scale_fill_grey() +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
    # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

(BASE.cv.plot <- ggplot(outBASE, aes(x=year, y=cv, fill=year))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    ggtitle("Baseline") +
    scale_fill_grey() +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
    # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

AllCV$metric <- as.factor(AllCV$metric)

(All.cv.plot <- ggplot(AllCV, aes(x=reorder(metric, -cv), y=cv, fill=metric))+  
    # Add a boxplot
    geom_bar(stat = "identity") +
    scale_fill_grey() +
    ggtitle("Variability across metrics") +
    labs(y = "Standardized CV (%)") +
    coord_flip() +
    # The colors to assign each time period
    # Black and white theme
    theme_classic() +
    # Remove the x-axis title
    theme(legend.position = "none",axis.title.y=element_blank()))

CORT.cv.plot + BC.cv.plot + Forage.cv.plot + DEX.cv.plot + ACTH.cv.plot + BASE.cv.plot

#save it
ggsave(All.cv.plot, 
       file = "Variability of Metrics.jpg", 
       width = 6, height = 6, unit = "in", dpi = 300)
