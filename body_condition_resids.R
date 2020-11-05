library(lubridate)
# Bring in the data
hare<-read.csv("body_condition_males.csv",header=TRUE)

#Assumption 1: The data is linear 
#Plot the data to make sure it follows a linear trend 
plot(weight~RHF, data=hare)

#Split by Sex
males<-hare[hare$sex=="M",]

# Model 1 Regression
maleBC = lm(weight ~ RHF, data=hare) 

#Put the residuals in a table so you can plot them
males.res<-resid(maleBC)


#All hares
BC <- lm(Weight ~ RHF, data = hare)
resids <- resid(BC)
# Assumpition 2: Homogenaeity of variance
#Plot the residuals to make sure they have a "shotgun" patturn
plot(hare$RHF, males.res) 
abline(0, 0)

#Assumpiton 3: There is normality
# Make a Q-Q plot of the residuals
qqnorm(males.res);
qqline(males.res)

#Create a column in the data set with the residuals from the model
hare$residuals <- males.res

#add month and year
hare$month <- lubridate::month(hare$date)
hare$year <- lubridate::year(hare$date)
#create a file with your data - residuals are the BC index!
write.csv(hare, "condition_males.csv", row.names = F)
