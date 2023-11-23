# title: "UA_2023_CPUE_STANDARDIZATION"
# author: "MARÍA SOTO"
# date: "2023-11-03"

####################################
## CPUE STANDARDIZATION EXERCISE  ##
####################################

# In this exercise we  use data from observers on board for the Spanish fresh
# trawling fleet in Mauritania targeting Black Hakes (Merluccius spp.)
# Data from observers are geo-referenced by fishing operation (hauls).
# Also information about the Sea Surface Temperature (SST) is available.
# The goal is to illustrate the main objectives and outputs for the 
# CPUE standardization process of commercial data.
# Depending on data availability, some models can deal better than others with 
# spatial structure in the observations. 
# We start with the simpler model based on GLM to illustrate the general 
# principles of the standardization.
# Later we will add complexity and show other models that can deal with spatial
# structure in more detail.
# We will explore example of models:

# 1. GLM Gamma model with fixed effects 
# 2. GLM Gamma model with vessel random effect
# 3. GAM Gamma model with area as lat-lon.
# 4. Bayesian Hierarchical Spatio-temporal Hurdle model Gamma lonlat.

####################################
## Install packages
####################################
library(ggplot2)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(carData) # para effect
library(effects) # para effect
library(Hmisc)
library(MASS) # PARA fitdistr
library(carData)
library(car) # para qqp
library(fitdistrplus) # PARA fitdistr
library(emmeans)
library(scales) #para scales del ggplot
library(ggpubr) # para ggarrange
library(Matrix)# RANDOM EFFECTS
library(lme4) # RANDOM EFFECTS, GAM
library(mgcv) # GAM
library(gtools)
library(plyr)
library(dplyr)
library(ggstatsplot)
library(psych)
library(stats)
library(tidyr)
library(reshape2)


####################################
# load data
####################################
load("C:/MARIA_SOTO/MIS_CURSOS/UniversidadAlicante/2023/R/D4.RData")

# Filter observer data to haul targeting black hakes (strategy=="hake").
# This will probably give an index for matures

data <- subset(D4,strategy=="hake")

# load("./input/data.Rdata")

####################################
### set up 
####################################

## Set up dependent variable to be standardized
#  cpue = catch hake (Kg) / fishing effort (trawling hours)
data$cpue<-data$hake/data$trawling 

## Set up covariates 
data$year<-as.factor(data$year) # year as variable ID
data$year <- droplevels(data$year)
data$year_id<-as.numeric(data$year) # year as variable ID

## Create factors for the spatial variables lat & lon & depth
# based on the quantiles

library(gtools)
with(data,plot(lon,cpue))
Qlon <- quantcut(data$lon) # create a factor for longitude
with(data, plot(lat,cpue))
Qlat <- quantcut(data$lat) # create a factor for latitude
with(data, plot(depth,cpue))
Qdepth <- quantcut(data$depth) # create a factor for depth

library(plyr)
library(dplyr)

levels(Qlat)[1] <- '1'
levels(Qlat)[2] <- '2'
levels(Qlat)[3] <- '3'
levels(Qlat)[4] <- '4'

data$area <- Qlat

Qdepth <- quantcut(data$depth)

data$Qdepth <- Qdepth

with(data,plot(month,cpue))
# Based on the knowledge of the fishery and data exploratory analysis
# we divide month in fishing seasons: 
# 1 After Christmas (1,2,3)
# 2 Rest of the year (4,5,6,7,8)
# 3 at the end of the year near Christmas and spawning season (9,10,11,12)

data$season <- revalue(data$month,c("1"="1", "2"="1","3"="1",
                                    "4"="2","5"="2","6"="2","7"="2","8"="2",
                                    "9"="3","10"="3","11"="3","12"="3"))

data$season <- ordered(data$season, levels = c("1", "2", "3"))

with(data,plot(experience,cpue))
data$skill <- revalue(data$experience, c("1"="1","2"="2","3"="3","4"="3","5"="4"))

with(data,plot(vessel_size,cpue))
# Vessel size has a range with few values, that better fit into a factor 
# with vessel categories classified between intervals of vessel lenghts.
Qvessel_size <- quantcut(data$vessel_size) # create a factor for vessel size
levels(Qvessel_size)[1] <- '1'
levels(Qvessel_size)[2] <- '2'
levels(Qvessel_size)[3] <- '3'
levels(Qvessel_size)[4] <- '4'

data$Qvessel_size <- Qvessel_size

save(data,file = "./input/data.RData")

####################################
#### IMPORTANT 
####################################
# in this case we have used all data from the observers surveys
# If we are working with commercial data from a whole fleet
# it is very important to "clean" observations of the 
# opportunistic behavior of vessels that only appear in the fishery
# in short periods of time. This only contributes to add noise to the model
# and doesn´t add information about the trend in abundance of the population

# for this pupose we have to calculate the number of vessels each year
# and the number of days of each vessel each year.
# Then set a minimum limit of days per year and a minimum limit of years in the period.
# Then, validate that the selected vessels account for the majority of total catch
# (more than 90%)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Explore dependent variable CPUE distribution
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################
# check outliers https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
####################################

# if Q1= 25th percentile, Q3= 75th percentile, then, IQR= Q3 – Q1
# An outlier would be a point below [Q1- (1.5)IQR] or above [Q3+(1.5)IQR].

library(ggstatsplot)

## look for outliers in the whole cpue distribution
outliers <- boxplot(data$cpue, plot=FALSE)$out
length(outliers)
x<-data 
x<- x[-which(x$cpue %in% outliers),]
boxplot(x$cpue)
newdata <- x

## if we want to analise the distribution of observations that has been 
# identified as outliers...
outliersdata <- data[which(data$cpue %in% outliers),]
hist(outliersdata$cpue)
boxplot(outliersdata$cpue)
with(outliersdata, plot(year, cpue))
length(outliersdata$cpue)

## look for outliers year by year
with(data, plot(year, cpue))

require(stats)
out.year <- by(data, data[,"year"], 
               function(x) boxplot(x$cpue, plot=FALSE)$out)
x <- data.frame(unlist(out.year))
outliers2 <- x$unlist.out.year.
length(outliers2) 
# the number of outliers by year is higher than outliers in global CPUE
# we remove outliers from the whole cpue

# once we have eliminated the outliers we explore the distribution for
# the resulting CPUE
dim(newdata) # we have 1333 observations
hist(newdata$cpue) # distribution is asimetric
boxplot(newdata$cpue)

summary(newdata$cpue) # probably values in the first quantile will give problems in model fitting
length(which(newdata$cpue ==0)) # No 0 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT OF DISTRIBUTIONS BY MAXIMUM LIKELIHOOD ESTIMATION
# https://rpubs.com/aafernandez1976/fitdistrplus
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# there will be several model distributions that could fit reasonably well the CPUE
# Gamma is a good candidate for x>0
# We first observe the empirical distrubution of the data

plotdist(newdata$cpue, histo=TRUE, demp=TRUE)

# Then, we explore the statistics of kurtosis and skewness.
# If skewness != 0, then the empirical distribution is asymmetric
# kurtosis quantifies the weight of the extremes and can be comparable to the
# value of kurtosis of a normal distribution that is 3.

descdist(newdata$cpue)

# Criteria: The density plot and the CDF plot may be considered as the basic 
# classical goodness-of-fit plots. 
# QQPlot  and PP plot are complementary and can be very informative in some cases. 
# The Q-Q plot emphasizes the lack-of-fit at the distribution tails 
# The P-P plot emphasizes the lack-of-fit at the distribution center.

library(MASS)
library("fitdistrplus")
a = fitdistr(newdata$cpue,"gamma") 
b = fitdistr(newdata$cpue,"lognormal")


# Creating x-values for density 
u <- seq(1:length(newdata$cpue))

# pgamma gives the Gamma cummulative distribution function with parameters 
# shape = α and scale = σ
vab=pgamma(u,a$estimate[1],a$estimate[2])

# plnorm gives the lognormal cummulative distribution function with parameters 
# mu and σ
vms=plnorm(u,b$estimate[1],b$estimate[2])

# empirical cummulative distribution function

plot(ecdf(newdata$cpue),main = "ecdf(CPUE)",xlab="CPUE",ylab="Fn(CPUE)")
lines(u,vab,col="red",lty=1) # Gamma best fit
lines(u,vms,col="blue",lty=2)
legend(6,0.4, legend=c("x","Gamma", "lognormal"),
       col=c("black","red", "blue"), lty=c(1,1,2), cex=0.5,box.lty=0)



save(newdata,file = "./input/newdata.RData")

#####################################
# We fit a Gamma model for the CPUE
#####################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GLM 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We select columns of explanatory variables 

# lets try with Depth and vessel_size as continous variables
# and the factorized version Qdepth and Qvessel_size

# continuous
glmvarscont <- c("cpue","year","season","area","depth","SST","vessel_size","skill","vessel")
glmdatacont <- data.frame(newdata[,names(newdata) %in% glmvarscont])
save(glmdatacont,file = "./input/glmdatacont.RData")


# Factors
glmvars <- c("cpue","year","season","area","Qdepth","SST","Qvessel_size","skill","vessel")
glmdata <- data.frame(newdata[,names(newdata) %in% glmvars])
save(glmdata,file = "./input/glmdata.RData")

#####################################
# Explore visual relationships/correlations in glmdata
#####################################

library(psych)
pairs.panels(glmdatacont, 
             method = "pearson", # Method for correlations
             hist.col = "aquamarine",
             density = TRUE,  # show densities
             ellipses = FALSE) #Ignore correlation ellipses 
# From the plot the relationship between CPUE and covariates is not obvious

pairs.panels(glmdata, 
             method = "pearson", # Method for correlations
             hist.col = "aquamarine",
             density = TRUE,  # show densities
             ellipses = FALSE) #Ignore correlation ellipses 

#####################################
## Explore TIME-SPACE interactions
#####################################

# 1. season:area
season_area.df <- glmdata[,c("season","area","cpue")]
m1<-lm(cpue ~ season + area, data = season_area.df)    
m2<-lm(cpue ~ season * area, data = season_area.df)    
s_a <- effect('season*area', m2,se=TRUE) #The Interaction
seasonarea.df<-as.data.frame(s_a) #Data Frame

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))

sa <- ggplot(data=seasonarea.df, aes(x=season, y=fit, group=area))+
  geom_line(size=1, aes(linetype=area,color=area))+
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se,fill= area),alpha=.2)+
  ylab("CPUE")+
  xlab("season")+
  ggtitle("CPUE interaction season-area")+
  theme_bw()+
  theme(text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")+
  theme(axis.text.x = element_text(face="bold", color="black", size=14))+
  theme(axis.text.y = element_text(face="bold", color="black", size=14))
sa
ggsave("SeasonArea.pdf")

# 2. season:depth

season_depth.df <- glmdatacont[,c("season","depth","cpue")]
m1<-lm(cpue ~ season + depth, data = season_depth.df)    
m2<-lm(cpue ~ season * depth, data = season_depth.df)    
s_a <- effect('season*depth', m2,se=TRUE) #The Interaction
seasondepth.df<-as.data.frame(s_a) #Data Frame

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))

sd <- ggplot(data=seasondepth.df, aes(x=depth, y=fit, group = season))+
  geom_line(size=1, aes(linetype=season,color=season))+
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se,fill= season),alpha=.2)+
  ylab("CPUE")+
  xlab("depth")+
  ggtitle("CPUE interaction season-depth")+
  theme_bw()+
  theme(text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top")+
  theme(axis.text.x = element_text(face="bold", color="black", size=14))+
  theme(axis.text.y = element_text(face="bold", color="black", size=14))
sd
ggsave("SeasonDepth.pdf")

# Graphs show evidences that we have to consider interactions in the models
# We can explore all the interactions that we consider that could influence
# commercial CPUE

library(lme4)

#####################################
# FIXED EFFECTS
#####################################

# Consider different GLMs without/with time-space interactions and
# continuous/factors Qdepth and Qvessel_size

# formulas
fglm1 <- cpue ~ year + season + Qdepth + area + SST + Qvessel_size + skill
fglm2 <- cpue ~ year + SST + Qvessel_size + skill+ season:area + season:Qdepth
fglm3 <- cpue ~ year + Qdepth + Qvessel_size + skill+ season:area + season:Qdepth # remove SST (we will see later that is not significant)
fglm4 <- cpue ~ year + season + depth + area + SST + vessel_size + skill
fglm5 <- cpue ~ year + SST + vessel_size + skill + season:area + season:depth
fglm6 <- cpue ~ year + depth + vessel_size + skill + season:area + season:depth

glm1 <- glm(fglm1,
            data = glmdata, 
            family  = Gamma(link = "log"))


glm2 <- glm(fglm2,
            data = glmdata, 
            family  = Gamma(link = "log"))

glm3 <- glm(fglm3,
            data = glmdata, 
            family  = Gamma(link = "log"))

glm4 <- glm(fglm4,
            data = glmdatacont, 
            family  = Gamma(link = "log"))

glm5 <- glm(fglm5,
            data = glmdatacont, 
            family  = Gamma(link = "log"))

glm6 <- glm(fglm6,
            data = glmdatacont, 
            family  = Gamma(link = "log"))

#####################################
# RANDOM EFFECTS for vessel and possible interactions
#####################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## GLMM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# formulas

fglmm1 <- cpue ~ year + Qvessel_size + season:Qdepth + season: area + (1|vessel) 
fglmm2 <- cpue ~ year + Qvessel_size + (1|season:Qdepth) + (1|season: area) + (1|vessel)
fglmm3 <- cpue ~ year + vessel_size + season:depth + season: area + (1|vessel)
fglmm4 <- cpue ~ year + vessel_size + (1|season:depth) + (1|season: area) + (1|vessel)

glmm1 <- glmer(fglmm1,
               data = glmdata, 
               family  = Gamma(link = "log"))

glmm2 <- glmer(fglmm2,
               data = glmdata, 
               family  = Gamma(link = "log"))

glmm3 <- glmer(fglmm3,
               data = glmdatacont, 
               family  = Gamma(link = "log"))

glmm4 <- glmer(fglmm4,
               data = glmdatacont, 
               family  = Gamma(link = "log"))

#####################################
# DIAGNOSTIC AND MODEL COMPARISON FOR FIXED EFFECTS GLMs
#####################################

## Partial effects
library(effects)
x11()
plot(allEffects(glm1))
plot(allEffects(glm2))
plot(allEffects(glm3))
plot(allEffects(glm4))
plot(allEffects(glm5))
plot(allEffects(glm6))


## Significance of covariates

Deviance_table_func <- function(my_glm) {
  
  x <- anova(my_glm,test="Chisq")
  x <- as.data.frame(x)
  x$PercDevExp <- 100*(x$Deviance/(max(x[,4])-min(x[,4])))
  row.names(x)[x$PercDevExp >= 1.5]
  DevTable <- x
  Dev_expl <- (max(x[,4])-min(x[,4]))/max(x[,4])
  
  return(DevTable)
  return(Dev_expl)
}


Deviance_explained_func <- function(my_glm) {
  
  x <- anova(my_glm,test="Chisq")
  x <- as.data.frame(x)
  x$PercDevExp <- 100*(x$Deviance/(max(x[,4])-min(x[,4])))
  row.names(x)[x$PercDevExp >= 1.5]
  DevTable <- x
  Dev_expl <- (max(x[,4])-min(x[,4]))/max(x[,4])
  return(Dev_expl)
}

Deviance_table_func(glm1)
Deviance_table_func(glm2)
Deviance_table_func(glm3)
Deviance_table_func(glm4)
Deviance_table_func(glm5)
Deviance_table_func(glm6)


DevGLM <- data.frame(cbind(c("glm1","glm2","glm3","glm4","glm5","glm6"),
round(rbind(
Deviance_explained_func(glm1),
Deviance_explained_func(glm2),
Deviance_explained_func(glm3),
Deviance_explained_func(glm4),
Deviance_explained_func(glm5),
Deviance_explained_func(glm6)),digits = 4)))
names(DevGLM) <- c("model","DevExpl")
DevGLM

# glm5 explain 41% of the variance in the data

#####################################
# MODEL COMPARISON FOR RANDOM EFFECTS GLMMs
#####################################

anova(glmm1,glmm2) # Better glmm1
anova(glmm3,glmm4) # Better glmm4

#####################################
# FIXED EFFECTS VS. RANDOM EFFECTS
#####################################
# COMPARE ALL MODELS PREDICTIVE AKAIKE INFO CRITERION. 
# Predictive capacity: LOWER AIC, THE BEST MODEL
# You cannot use likelihood-based statistics like AIC to compare across models 
# with different likelihood functions - the underlying formulas are different!!!

library(stats)
GLMs <- data.frame(AIC(glm1,glm2,glm3,glm4,glm5,glm6,glmm1,glmm2,glmm3,glmm4))
GLMs

#######################
# Better model glmm4 
#######################

#####################################
# DIAGNOSTIC MODEL FITTING FOR RANDOM EFFECTS GLMMs
#####################################

summary(glmm4)
# Residuals Diagnostic plots
x11()
par(mfrow=c(2,2))

hist(glmdatacont$cpue)
boxplot(glmdatacont$cpue)
hist(resid(glmm4),xlab="residuals",main="",col="green")
qqnorm(resid(glmm4),main="")
qqline(resid(glmm4),col="red",main="")

# The model is good in the distribution center and fails in the tails.
# The way to improve the model is to increase the number of observations
# and inspect other sources of information (environmental, fishing power, etc)


###########################
## SELECTED MODEL GLMM4  ##
###########################

# We will use glmm4 to construct the standardized CPUE
# from the glmm
#####################################
# ESTIMATE LS-MEANS FROM SELECTED GLMM
#####################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LS-MEANS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Transformación de los coeficientes del modelo que pondera por
## el número de observaciones de cada estrato (year-vessel_size) 
## y calcula el promedio anual cuando el muestreo no es uniforme. 

ls<-lsmeans(glmm4,"year",type="response") # type="response" estimates LSMeans
# in the scale of the response (cpue), i.e. exp(fitted values)
ls
plot(ls)
LSM.df <- as.data.frame(ls)
LSM.df$year <- as.numeric(as.character(LSM.df$year))

# Si queremos comparar el índice con la CPUE nominal, debemos tener cuidado con 
# la escala. La escala de la CPUE nominal debe estar en la misma escala que la 
# CPUE estandarizada para poder comparar.
# Al utilizar la opción type="response" en lsmeans estamos trabajando con la 
# escala original. Como hemos seguido un modelo Gamma con link = log, si 
# trabajamos en la escala de los fitted values entonces las LSMeans están 
# también en escala logarítmica. En ese caso, PRIMERO SE HARÍA EL LOGARITMO Y 
# LUEGO LA MEDIA ANUAL SOBRE CPUE NOMINAL

# Function to arrange columns for lsmeans and nominal values

CPUE_nominal_DF_function <- function(data,LSM.df){
  
  nominal <- aggregate.data.frame(data$cpue,list(year=data$year),FUN=mean)
  names(nominal) <- c("year","CPUE_nom")
  nominal$sd_nom <- as.vector(aggregate.data.frame(data$cpue,list(year=data$year),FUN=sd)[,2])
  estimate <- cbind(nominal,LSM.df[,-1])
  estimate$lwr <- with(estimate,lsmean-1.96*SE)
  estimate$upr <- with(estimate,lsmean+1.96*SE)
  return(estimate)
}

estimate <- CPUE_nominal_DF_function (glmdata,LSM.df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot standardized and nominal CPUE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# funtion to put years in numeric without decimals in the ggplot.
# Redondeamos year para que aparezca entero en el eje x con la función
# int_breaks

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 

estimate$year <- as.numeric(as.character(estimate$year))

I<-ggplot(estimate) + 
  geom_line(aes(year, lsmean),color = "blue",size=2) +
  geom_ribbon(aes(year, ymin = lwr, ymax = upr), alpha = .2) +
  geom_point(aes(year, y = CPUE_nom), color = "red",size=2)+
  labs(x = "year",y = "kg/fishing hours",title = "Black hakes CPUE")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold", color="black", size=14),plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks= int_breaks)
I
# I+geom_smooth(aes(year, lsmean), color="red",methods="loess")

estimate
write.table(estimate,file="./output/index.csv")


######################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## GAM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If we have continuous variables like lat and lon or any other, one alternative
# is to use GAMs.  
# GAM models spatial information as actual fishing locations rather than areas.
# The latitude and longitude points can be smoothed as tensor products, 
# te(lat, lon).
# We explore the time-space interaction is modeled as an approximation to 
# capture potential paterns in space and over time through a tensor 
# te(lat, lon, seasonday) term.
# Seasonday is the continuous date along the year.

#################################
### Select columns for GAM models
#################################
gamvars <- c("cpue","year","seasonday","lat","lon","depth","SST","vessel_size","skill","vessel")
gamdata <- data.frame(newdata[,names(newdata) %in% gamvars])
save(gamdata, file="./input/gamdata.RData")

#################################
## Explore correlations
#################################
pairs.panels(gamdata, 
             method = "pearson", # Método para el cálculo de las correlaciones
             hist.col = "aquamarine",
             density = TRUE,  # Mostrando las densidades
             ellipses = FALSE #Ignorando las elipses de correlación
)

#################################
## fit GAM models
#################################

# To fit GAM models we will use 80% of the data for training the model and
# 20% as test data.

set.seed(123)
train.index <- sample(1:nrow(gamdata), 0.8 * nrow(gamdata))
test.index <- setdiff(1:nrow(gamdata), train.index)
gamtrain <- gamdata[train.index,]
gamtest <- gamdata[test.index,]

library(mgcv)

# formulas

fgam1 <- cpue ~ year + s(seasonday) + s(depth) + te(lat, lon)+ vessel_size + skill
fgam2 <- cpue ~ year + s(seasonday) + s(depth) + te(lat, lon)+ s(vessel, bs = 're')
fgam3 <- cpue ~ year + s(depth) + te(lat, lon, seasonday)+ s(vessel, bs = 're')
fgam4 <- cpue ~ as.numeric(year) + s(seasonday) + s(depth) + 
                te(lat, lon, as.numeric(year)) + s(vessel, bs = 're')   
# models

gam1 <- gam(fgam1, 
            family = Gamma(link = log), data = gamtrain)

gam2 <- gam(fgam2, 
            family = Gamma(link = log), data = gamtrain)

gam3 <- gam(fgam3, 
            family = Gamma(link = log), data = gamtrain)

gam4 <- gam(fgam4, 
            family = Gamma(link = log), data = gamtrain)

summary(gam1)
summary(gam2)
summary(gam3)
summary(gam4)

#################################
# Compare GAM models
#################################

# statistically we can compare through ANOVA
anova(gam1, gam2, gam3, gam4, test ="Chisq")

# now we have additional statistical evidence to suggest that incorporating
# random effect vessel and interactions time-space improves the model.

summary(gam1)
summary(gam2)
summary(gam3)
summary(gam4)

# model gam 3 accounts for much of the variance in cpue, with an adjusted R-squared of 
# R-sq.(adj) =  0.417 and  Deviance explained =   47%

# GCV, similar to AIC 
c(gam1$gcv.ubre,gam2$gcv.ubre,gam3$gcv.ubre,gam4$gcv.ubre)

#################################
# Compare predition skill for 
# glm, glmm and gam
#################################

AIC(glm1,glm2,glm3,glm4,glm5,glm6,glmm1,glmm2,glmm3,glmm4,gam1,gam2,gam3,gam4)

#################################
# Visualize gam effects
#################################

plot(gam3)
plot(gam4)

# GAM improves model prediction with vessel as random effect

#################################
# Diagnostic GAM
#################################
# model fitting --> residuals
# model prediction --> AIC

x11()
par(mfrow = c(2,2))
hist(resid(gam3),breaks= 30, xlab="residuals",main="",col="green")
qqnorm(resid(gam3),main="gam3")
qqline(resid(gam3),col="red",main="")

hist(resid(gam4),breaks= 30, xlab="residuals",main="",col="green")
qqnorm(resid(gam4),main="gam4")
qqline(resid(gam4),col="red",main="")

# Residuals show some fitting problems in the tails of the distribution

#################################
# prediction
#################################
# validate prediction of gam models for test data

gam3pred <- predict.gam(gam3, newdata = gamtest)
header <- c("year","gam3.p", "gam3.o") ## Names of columns
datos.o.p <- data.frame(as.integer(as.character(gamtest[,"year"])),gam3pred,log(gamtest[,"cpue"])) 
colnames(datos.o.p) <- header ##Name the columns
library(tidyr)
library(reshape2)
datos.op.long <- melt(datos.o.p, id.vars = c("year"))

#Construcción del gráfico de predicción vs. observado

g3 <- ggplot(datos.op.long, aes(x=year, y=value, color=variable,shape=variable)) +
  geom_point() +
  scale_shape_manual(values = c(3, 5)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  scale_alpha_manual(values = c(0.1, 1)) +
  theme_minimal() +
  theme(legend.position = "top")
# +facet_wrap("variable")
x11()
g3

# For some years prediction works better than others

# constrution of the CPUE index
# first we replicate for all years in the data all combinations of the locations
# and year, including all those not fished in a particular year,
# the mean value of depth, mean value of seasonday and the
# vessel 1 as it is the most representative

summary(gamdata$vessel)

# create a combined variable of (lat,lon) in gamdata

library(dplyr)
dfpred1 <- gamtrain %>% 
  mutate(ID = group_indices_(gamtrain, .dots=c("lat", "lon"))) 
dim(dfpred1)
mean(gamdata$seasonday)
mean(gamdata$depth)

#################################
# prediction grid
#################################

# for prediction we have to set a reference to fix covariates for all
# possible locations each year. To fis a reference vessel we look for the one 
# that isset the most representative, i.e. present most of the years
table(gamdata$year,gamdata$vessel)
# in this case v1

dim(dfpred1)
dfpred2 <- data.frame(year = as.factor(as.character(
                      rep(c("2007","2009","2010","2011","2016","2017","2018","2019","2021"),
                                                              each=1066))),

                      depth = as.double(rep(mean(gamdata$depth),1066*9)),
                      lat = as.double(rep(dfpred1$lat,times=9)),
                      lon = as.double(rep(dfpred1$lon,times=)),
                      seasonday = as.double(rep(mean(gamdata$seasonday),1066*9)),
                      vessel = as.factor(as.character(rep("v1",1066*9))))

gam3pred <- predict.gam(gam3, newdata = dfpred2,type= "response",se=TRUE)
gam3pred[["fit"]]
gam3pred[["se.fit"]]
gam3estimate <- cbind(dfpred2,gam3pred[["fit"]],gam3pred[["se.fit"]])
names(gam3estimate) <- c( "year" ,"depth" ,"lat" ,"lon" ,"seasonday" ,"vessel", "fit" ,"se")

#################################
## CPUE index
#################################

# Calculate the CPUE index as the mean values across year
Igam <- data.frame(year = unique(gam3estimate$year))
Igam$Igammu <- as.vector(with(gam3estimate,aggregate.data.frame(fit,list(year=year),FUN=mean)[,2]))
Igam$Igamsd <- as.vector(with(gam3estimate,aggregate.data.frame(se,list(year=year),FUN=sd)[,2]))
Igam$lowgam <- with(Igam, Igammu - 1.96 * Igamsd)
Igam$upgam <- with(Igam, Igammu + 1.96 * Igamsd)

# Add to GLMM CPUE estimate GAM CPUE

Iglmmgam <- cbind(estimate,Igam[,-1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot standardized and nominal CPUE 
# GLMM and GAM indices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Iglmmgam$year <- as.numeric(as.character(Iglmmgam$year))

colors <- c("GAM" = "blue", "GLM" = "red", "nominal" = "black")
I <- ggplot(Iglmmgam) + 
  geom_line(aes(year, Igammu,color = "GAM"),size=2) +
  geom_ribbon(aes(year, ymin = Igammu-1.96*Igamsd, ymax = Igammu+1.96*Igamsd,color="GAM"), alpha = .2) +
  geom_line(aes(year, lsmean,color = "GLM"),size=2) +
  geom_ribbon(aes(year, ymin = asymp.LCL, ymax = asymp.UCL,color="GLM"), alpha = .2) +
  geom_point(aes(year,CPUE_nom,color = "nominal"),size=2) +
  labs(x = "year",y = "kg/fishing hours",title = "Black hakes CPUE")+
  scale_color_manual(values = colors, name = "") +
  theme_classic()+
  theme(axis.text.x = element_text(face="bold", color="black", size=14),
        plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks= int_breaks)
I

write.table(Iglmmgam,file="./output/indexglmmgam.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# APPENDIX INLA Model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Code to implement INLA but not prediction
# year term indicates the trend in CPUE

###################################################################
# Load support files and packages

library(sp)
library(INLA)
library(rgeos)
library(raster)
library(fields)
library(ggplot2)
library(ggspatial)
library("rnaturalearth")

# AREA DE ESTUDIO ###############################

Mau <- getData('GADM',country="MRT",level=0)
ext<-extent(-19.01852, -15.98148, 15.98154, 21.34012) # pruebo con este
Mau <- crop(Mau, ext) # recorta la parte de Mauritania dentro de ese rectángulo. 
# Por la parte de tierra limita el rectángulo y por la derecha la costa de Mauritania.
x11()
plot(Mau)
points(gamdata$lon, gamdata$lat)
#Este es el polígono de la tierra.

#  2) Defino el área de estudio, es decir una caja rectangular de la costa y los puntos. 
#  hay que jugar un poco haciendo el xym porque hay que definir una matriz de 4 
#  puntos que coincidan para cerrar la caja:
#  es decir, un rectángulo tendrá parte de Mauritania y parte de mar que incluya los lances
#xym<- as.matrix(data.frame(x =c(-19, -16, -16,-19),
#                          y = c(16,16,21.5,21.5)))

xym<- as.matrix(data.frame(x =c(-19.01852, -15.98148, -15.98148,-19.01852),
                           y = c(15.98154,15.98154,21.34012,21.34012)))
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

x11()
plot(sps)
plot(Mau, add=T, col="grey")
points(gamdata$lon, gamdata$lat)

Mau_rec<-crop(Mau, sps)
plot(Mau_rec)
proj4string(sps)<-proj4string(Mau_rec)

## --- Defino el cuadrado para pintar --- ###

p2 = Polygon(xym)
ps2 = Polygons(list(p2),1)
sps2 = SpatialPolygons(list(ps2))

plot(sps2)
plot(Mau_rec, add=TRUE)

#  3) se hace la diferencia entre estos dos poligonos para crear un tercer 
# poligono espacial que será el de mar

### --- Selección del polígono que contiene los datos de los lances --- ###
coast <- gDifference(sps2, Mau_rec) # superpone la TIERRA en sps2 y eso es mar
plot(coast, col="blue")


# mesh 11

# distancias en º
Loc <- cbind(gamdata$lon,gamdata$lat)
D <- dist(Loc)

bound=bound<-inla.sp2segment(coast)
mesh11 <- inla.mesh.2d(loc=as.matrix(Loc),boundary = bound,
                       offset = c(0.5, 0.5),
                       cutoff = 0.05, max.edge = c(0.17, 0.4))
save(mesh11,file="./input/MESH.RData")
save(Mau_rec,file="./input/Mau_rec.RData")

# MESH
x11()
plot(mesh11,main="")
points(gamdata[,6:5], col="red", pch=20) # 17:16 son lon lat de D5
plot(Mau_rec, add=TRUE,  col="gray")


# stack ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gamvars <- c("cpue","year","seasonday","lat","lon","depth","SST","vessel_size","skill","vessel")
gamdata <- data.frame(newdata[,names(newdata) %in% gamvars])
save(gamdata, file="./input/gamdata.RData")

data <- gamdata[,c("cpue","lon","lat")]

# SPDE ------------------------------------------------------------------------------------------------
spde  <- inla.spde2.pcmatern(mesh = mesh11, alpha = 1.5,
                             prior.range = c(2, 0.6),
                             prior.sigma = c(20, 0.01))

A.est <- inla.spde.make.A(mesh11, loc=cbind(gamdata$lon, gamdata$lat)) 

# stack
stk.est<-inla.stack(data=list(y=data$cpue),
                      A=list(A.est, 1),
                      effects=list(spatial=1:spde$n.spde,
                                   data.frame(beta0=1, gamdata)),
                      tag="est",compress = TRUE, remove.unused = TRUE)

# fórmula
hyper.prior = list(theta1 = list(prior="pc.prec", param=c(0.06, 0.008)),
                   theta2 = list(prior="pc.cor1", param=c(0.9, 0.9)) )

f <- y~-1 + beta0 + f(spatial, model=spde) +
  f(depth,model="rw1") + f(year,model = "rw1") + 
  seasonday + SST + 
  vessel_size + skill + 
  f(vessel, model ='iid') 

#modelo

m1 <- inla(f, 
               data=inla.stack.data(stk.est), family="gamma",
               control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
               control.predictor=list(A=inla.stack.A(stk.est), compute=TRUE, 
                                      quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
               num.threads = 3,
               verbose=T)

inla.spde.make.index(name = "year",n.spde = mesh11$n,n.group = n_groups)

################################################################################
# Mapas de media y sd espaciales
bbox(coast)
(dxy <- apply(bbox(coast),1, diff)) # longitud de los rangos de lat y lon
(r <- dxy[1]/dxy[2]) # ratio entre le rango de log y lat
m<-150
proj.grid.mat <- inla.mesh.projector(mesh11, # dentro de mesh9 establece un rango 
                                     xlim=bbox(coast)[1,], # donde vamos a predecir en puntos
                                     ylim=bbox(coast)[2,], # que no son los de gamdata, son nuevos
                                     dims=c(r, 1)*m)
x11()
plot(coast)
points(proj.grid.mat$lattice$loc, pch=20, cex=0.5)
# limpiar (fijar NA para los valores fuera de boundary)
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, coast@proj4string), coast)
# comprobar los puntos de la red dentro de map
i.map <- is.na(ov)
# puntos sobre los que vamos a predecir
par(mar=c(0,0,0,0))
plot(sps)
points(proj.grid.mat$lattice$loc[!i.map,], col="red", cex=0.2)
points(proj.grid.mat$lattice$loc[i.map,], col="blue", cex=0.2)

#consideramos solo los que están dentro de map
proj.grid.mat$lattice$loc[i.map, ]

# variables espaciales del modelo
mean.g <- inla.mesh.project(proj.grid.mat, m1$summary.random$spatial$mean)
sd.g <- inla.mesh.project(proj.grid.mat, m1$summary.random$spatial$sd)
quantile_0.025 <- inla.mesh.project(proj.grid.mat, m1$summary.random$spatial$`0.025quant`)
quantile_0.975 <- inla.mesh.project(proj.grid.mat, m1$summary.random$spatial$`0.975quant`)

sd.g[i.map] <- mean.g[i.map] <- quantile_0.025[i.map] <- quantile_0.975[i.map] <- NA

### --- Spatial effect --- ###
### --- Posterior mean --- ###
x11()
plot(sps2, col="gray", main="Spatial mean CPUE gamma")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           mean.g, add=TRUE)
plot(Mau_rec, add=TRUE)

sp.mean.raster.m1<-raster(list(x = proj.grid.mat$x, 
                                   y = proj.grid.mat$y,
                                   z = mean.g)) 

jpeg("CPUE_mean.jpeg", width = 950, height = 1500, res = 300)
plot(sp.mean.raster.m1, col=tim.colors(100)[1:100],main=" ", axes=T)
plot(Mau_rec,add=T, axes=TRUE,col='dark grey')

dev.off()

writeRaster(sp.mean.raster.m1, filename="sp_mean_m1.asc", 
            format="ascii", overwrite=TRUE) # Luego se puede leer en ArcGis

#spatial effect sd
sp.mean.raster.m1_sd<-raster(list(x = proj.grid.mat$x, 
                                      y = proj.grid.mat$y,
                                      z = sd.g)) 
### --- Posterior sd --- ###
x11()
plot(sps2, col="gray", main="Spatial sd CPUE")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           sd.g, add=TRUE)
plot(Mau_rec, add=TRUE)
dev.off()

jpeg("CPUE_sd.jpeg", width = 950, height = 1500, res = 300)
x11()
plot(sp.mean.raster.m1_sd, col=tim.colors(100)[1:100],main=" ", axes=T)
plot(Mau_rec,add=T, axes=TRUE,col='dark grey')

dev.off()

writeRaster(sp.mean.raster.m1_sd, filename="sp_mean_m1_Sd.asc", 
            format="ascii", overwrite=TRUE) # Luego se puede leer en ArcGis

# Figura de la media y la sd juntas
jpeg("CPUE.jpeg", width = 950, height = 1500, res = 300)
x11()
par(mfrow=c(1,2))
plot(sp.mean.raster.m1, col=tim.colors(100)[1:100],main="Spatial mean CPUE", axes=T)
plot(Mau_rec,add=T, axes=TRUE,col='dark grey', main="")
plot(sp.mean.raster.m1_sd, col=tim.colors(100)[1:100],main="Spatial sd CPUE", axes=T)
plot(Mau_rec,add=T, axes=TRUE,col='dark grey')
##########################################################################

# los estadísticos del modelo
summary(m1)

x11()
plot(m1)
dev.off()

# probabilidad de que los betas estén dentro del 0
1-inla.pmarginal(0, m1$marginals.fixed$beta0)
1-inla.pmarginal(0, m1$marginals.fixed$seasonday) 
1-inla.pmarginal(0, m1$marginals.fixed$SST)
1-inla.pmarginal(0, m1$marginals.fixed$vessel_size)
1-inla.pmarginal(0, m1$marginals.fixed$skill2)
1-inla.pmarginal(0, m1$marginals.fixed$skill3)
1-inla.pmarginal(0, m1$marginals.fixed$skill4)

# plot posteriors
x11()
par(mfrow=c(2,2))
plot(m1$marginals.fixed$beta0, type="l",main="CPUE Intercept")
abline(v=0, col="red", lwd=2)
plot(m1$marginals.fixed$seasonday, type="l",main="seasonday")
abline(v=0, col="red", lwd=2)
plot(m1$marginals.fixed$SST, type="l",main="SST")
abline(v=0, col="red", lwd=2)
plot(m1$marginals.fixed$vessel_size, type="l",main="vessel_size")
abline(v=0, col="red", lwd=2)

#batimetria
x11()
Smoother <- m1$summary.random$depth
plot(x = Smoother[,"ID"],
     y = Smoother[,"mean"],
     type = 'l', lwd = 3,
     xlab = 'Bathymetry', 
     ylab = 'Smoother',
     ylim = c(-1, 1),
     main = "depth")
abline(h = 0, lty =3)
lines(x = Smoother[, "ID"],
      y = Smoother[, "0.025quant"], lty = 2)
lines(x = Smoother[, "ID"],
      y = Smoother[, "0.975quant"], lty = 2)

dev.off()

################################################################################
# year as RW1
x11()
Smoother2 <- m1$summary.random$year
plot(x = Smoother2[,"ID"],
     y = Smoother2[,"mean"],
     type = 'l', lwd = 3,
     xlab = 'Year', 
     ylab = 'Smoother',
     ylim = c(-0.7, 1),
     main = "year")
abline(h = 0, lty =3)
lines(x = Smoother2[, "ID"],
      y = Smoother2[, "0.025quant"], lty = 2)
lines(x = Smoother2[, "ID"],
      y = Smoother2[, "0.975quant"], lty = 2)

# This will be used as a proxy for the CPUE (not the best way, only for illustrative purposes)
# Prediction with INLA is complex and takes time ... another course

A.est@Dim
