##############################################################
# Author: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# MPIO (Seewiesen) and ICL (Silwood Park) 
# Email: alfredo.tojar@gmail.com

# Script first created in August 2017
# inspiration: from Shinichi Nakagawa and Malgorzata Lagisz's code: https://github.com/mlagisz

##############################################################
# Description of script and Instructions
##############################################################

# This script is to run multilevel meta-analyses and meta-regressions
# to test whether bib size correlates with dominance rank in
# house sparrows across populations, research groups, etc.
# The script also generates plots for those analyses.

# Meta-analytic challenge to a textbook example of status
# signalling: publication trends and biases

# Sanchez-Tojar et al. In prep.

# more info at the Open Science Framework: osf.io/cwkxb


##############################################################
# Packages needed
##############################################################

# load pacakges
library(MCMCglmm)

# Clear memory
rm(list=ls())

##############################################################
# Functions needed
##############################################################

# function to convert Zr to r        
Zr.to.r<-function(Zr){
  r<-(exp(2*Zr)-1)/(exp(2*Zr)+1)
}


##############################################################
# Data
##############################################################

meta1 <- read.table("processed_data/metadatasets/Meta1.csv",
                    header=TRUE,sep=",")

meta2 <- read.table("processed_data/metadatasets/Meta2.csv",
                    header=TRUE,sep=",")

meta3 <- read.table("processed_data/metadatasets/Meta3.csv",
                    header=TRUE,sep=",")

meta4 <- read.table("processed_data/metadatasets/Meta4.csv",
                    header=TRUE,sep=",")

prior1 <- list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),
                                           G2=list(V=1,nu=0.002)))


##################################################################
# META 1: ORIGINAL PUBLISHED ESTIMATES: Supplementary Information
##################################################################

meta1.1 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta1$VZr,
                    data=meta1,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta1.2 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta1$VZr,
                    data=meta1,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta1.3 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta1$VZr,
                    data=meta1,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta1.1,meta1.2,meta1.3,
     file="processed_data/metadatasets/model_meta1_paper&pop.full.RData")

load("processed_data/metadatasets/model_meta1_paper&pop.full.RData")

m1 <- meta1.1
m2 <- meta1.2
m3 <- meta1.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))

#check chain autocorrelation
autocorr(m1$Sol[,1])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,1])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,1])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m3 #lowest DIC

plot(chosen)
summary(chosen)


###############################
# heterogeneity and variances #
###############################

# WI = weight, s2I = measurement error variance = sigma2m
WI <- na.omit(1/meta1$VZr)

s2I <- sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

total_var <- chosen$VCV[,"paperID"]+
  chosen$VCV[,"popID2"]+
  chosen$VCV[,"units"]+s2I # total variance (sigma2t)


total_het <- (total_var-s2I)/total_var # total heterogeneity I2
mean(total_het) 
HPDinterval(total_het)


#studyID heterogeneity
studyID_het <- chosen$VCV[,"paperID"]/total_var
mean(studyID_het)
HPDinterval(studyID_het)

#populationID heterogeneity
popID_het <- chosen$VCV[,"popID2"]/total_var
mean(popID_het)
HPDinterval(popID_het)


### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])


#############################################################
# PUBLICATION BIAS
#############################################################

# Egger's regressions performed on model residuals and measurement
# errors. Evidence for publication bias if intercept is different
# from 0.

chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta1$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta1$Precision<-sqrt(1/meta1$VZr)

meta1$MAR<-meta1$Zr-meta1$Prediction # meta-analytic residual!

meta1$zMAR<-meta1$MAR*meta1$Precision

# Egger's regression
Egger.meta1<-MCMCglmm(zMAR~Precision,family="gaussian",
                      verbose=FALSE,data=meta1,
                      nitt=2000000, thin=1800, burnin=200000)

save(Egger.meta1,
     file="processed_data/metadatasets/model_meta1_paper&pop_Egger.full.RData")

load("processed_data/metadatasets/model_meta1_paper&pop_Egger.full.RData")

plot(Egger.meta1)
autocorr(Egger.meta1$Sol)
summary(Egger.meta1)
HPDinterval(Egger.meta1$Sol[,1])



##############################################################
# META 2: FULL: Supplementary Information
##############################################################

meta2.1 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta2.2 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta2.3 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta2.1,meta2.2,meta2.3,
     file="processed_data/metadatasets/model_meta2_paper&pop.full.RData")


load("processed_data/metadatasets/model_meta2_paper&pop.full.RData")

m1 <- meta2.1
m2 <- meta2.2
m3 <- meta2.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,1])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,1])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,1])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m2 #lowest DIC

plot(chosen)
summary(chosen)


###############################
# heterogeneity and variances #
###############################

# WI = weight, s2I = measurement error variance = sigma2m
WI <- na.omit(1/meta2$VZr)

s2I <- sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

total_var <- chosen$VCV[,"paperID"]+
  chosen$VCV[,"popID2"]+
  chosen$VCV[,"units"]+s2I # total variance (sigma2t)


total_het <- (total_var-s2I)/total_var # total heterogeneity I2
mean(total_het)
HPDinterval(total_het)


#studyID heterogeneity
studyID_het <- chosen$VCV[,"paperID"]/total_var
mean(studyID_het)
HPDinterval(studyID_het)


#populationID heterogeneity
popID_het <- chosen$VCV[,"popID2"]/total_var
mean(popID_het)
HPDinterval(popID_het)


### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])


#############################################################
# PUBLICATION BIAS
#############################################################

# Egger's regressions performed on model residuals and measurement
# errors. Evidence for publication bias if intercept is different
# from 0.

chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta2$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta2$Precision<-sqrt(1/meta2$VZr)

meta2$MAR<-meta2$Zr-meta2$Prediction # meta-analytic residual!

meta2$zMAR<-meta2$MAR*meta2$Precision

# Egger's regression
Egger.meta2<-MCMCglmm(zMAR~Precision,family="gaussian",
                      verbose=FALSE,data=meta2,
                      nitt=2000000, thin=1800, burnin=200000)

save(Egger.meta2,
     file="processed_data/metadatasets/model_meta2_paper&pop_Egger.full.RData")

load("processed_data/metadatasets/model_meta2_paper&pop_Egger.full.RData")

plot(Egger.meta2)
autocorr(Egger.meta2$Sol)
summary(Egger.meta2)
HPDinterval(Egger.meta2$Sol[,1])



##############################################################
# META 2: FULL -> META-REGRESSION: TIME-LAG BIAS TEST
##############################################################

meta_reg2.1 <- MCMCglmm(Zr~scale(year),
                        random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta_reg2.2 <- MCMCglmm(Zr~scale(year),
                        random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta_reg2.3 <- MCMCglmm(Zr~scale(year),
                        random=~paperID+popID2,
                    mev=meta2$VZr,
                    data=meta2,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta_reg2.1,meta_reg2.2,meta_reg2.3,
     file="processed_data/metadatasets/model_meta_reg2_year_paper&pop.full.RData")


load("processed_data/metadatasets/model_meta_reg2_year_paper&pop.full.RData")

m1 <- meta_reg2.1
m2 <- meta_reg2.2
m3 <- meta_reg2.3


#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,c(1,2)])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,c(1,2)])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,c(1,2)])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m2 #lowest DIC

plot(chosen)
summary(chosen)


############################
# Estimates and R2marginal #
############################

### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])

mean(chosen$Sol[,2])
HPDinterval(chosen$Sol[,2])


# # Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(chosen$Sol[i,c(1:2)] %*% t(chosen$X)))
  vmVarF[i]<-Var}
R2m<-100*(vmVarF/(vmVarF+chosen$VCV[,"paperID"]+chosen$VCV[,"popID2"]+chosen$VCV[,"units"]))
mean(R2m)
HPDinterval(R2m)



##############################################################
# META 3: FULL: meta 1 in main text
##############################################################

meta3.1 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta3$VZr,
                    data=meta3,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta3.2 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta3$VZr,
                    data=meta3,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta3.3 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta3$VZr,
                    data=meta3,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta3.1,meta3.2,meta3.3,
     file="processed_data/metadatasets/model_meta3_paper&pop.full.RData")

load("processed_data/metadatasets/model_meta3_paper&pop.full.RData")

m1 <- meta3.1
m2 <- meta3.2
m3 <- meta3.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,1])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,1])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,1])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m1 #lowest DIC

plot(chosen)
summary(chosen)


###############################
# heterogeneity and variances #
###############################

# WI = weight, s2I = measurement error variance = sigma2m
WI <- na.omit(1/meta3$VZr)

s2I <- sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

total_var <- chosen$VCV[,"paperID"]+
  chosen$VCV[,"popID2"]+
  chosen$VCV[,"units"]+s2I # total variance (sigma2t)


total_het <- (total_var-s2I)/total_var # total heterogeneity I2
mean(total_het)
HPDinterval(total_het)


#studyID heterogeneity
studyID_het <- chosen$VCV[,"paperID"]/total_var
mean(studyID_het)
HPDinterval(studyID_het)


#populationID heterogeneity
popID_het <- chosen$VCV[,"popID2"]/total_var
mean(popID_het)
HPDinterval(popID_het)


### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])


#############################################################
# PUBLICATION BIAS
#############################################################

# Egger's regressions performed on model residuals and measurement
# errors. Evidence for publication bias if intercept is different
# from 0.

chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta3$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta3$Precision<-sqrt(1/meta3$VZr)

meta3$MAR<-meta3$Zr-meta3$Prediction # meta-analytic residual!

meta3$zMAR<-meta3$MAR*meta3$Precision

# Egger's regression
Egger.meta3<-MCMCglmm(zMAR~Precision,family="gaussian",
                      verbose=FALSE,data=meta3,
                      nitt=2000000, thin=1800, burnin=200000)

save(Egger.meta3,
     file="processed_data/metadatasets/model_meta3_paper&pop_Egger.full.RData")

load("processed_data/metadatasets/model_meta3_paper&pop_Egger.full.RData")

plot(Egger.meta3)
autocorr(Egger.meta3$Sol)
summary(Egger.meta3)
HPDinterval(Egger.meta3$Sol[,1])



###################################################################
# META 3: FULL -> META-REGRESSION: Biological moderators for meta 1
###################################################################

meta3$season.num <- ifelse(meta3$season=="nonbreeding",0,1)
meta3$int.num <- ifelse(meta3$interactions=="both",0,1)
meta3$inttype.num <- ifelse(meta3$interactiontype=="mix",0,1)


meta_reg3.1 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg3.2 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg3.3 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)

save(meta_reg3.1,meta_reg3.2,meta_reg3.3,
     file="processed_data/metadatasets/model_meta_reg3_paper&pop.full_inttype.RData")


load("processed_data/metadatasets/model_meta_reg3_paper&pop.full_inttype.RData")

m1 <- meta_reg3.1
m2 <- meta_reg3.2
m3 <- meta_reg3.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,c(1,2,3)])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,c(1,2,3)])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,c(1,2,3)])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m1 #lowest DIC

plot(chosen)
summary(chosen)


############################
# Estimates and R2marginal #
############################

### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])

mean(chosen$Sol[,2])
HPDinterval(chosen$Sol[,2])

mean(chosen$Sol[,3])
HPDinterval(chosen$Sol[,3])

mean(chosen$Sol[,4])
HPDinterval(chosen$Sol[,4])


# # Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(chosen$Sol[i,c(1:4)] %*% t(chosen$X)))
  vmVarF[i]<-Var}
R2m<-100*(vmVarF/(vmVarF+chosen$VCV[,"paperID"]+chosen$VCV[,"popID2"]+chosen$VCV[,"units"]))
mean(R2m)
HPDinterval(R2m)



##############################################################
# META 4: FULL: meta 2 in main text
##############################################################

meta4.1 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta4.2 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta4.3 <- MCMCglmm(Zr~1,random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta4.1,meta4.2,meta4.3,
     file="processed_data/metadatasets/model_meta4_paper&pop.full.RData")

load("processed_data/metadatasets/model_meta4_paper&pop.full.RData")

m1 <- meta4.1
m2 <- meta4.2
m3 <- meta4.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))

#check chain autocorrelation
autocorr(m1$Sol[,1])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,1])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,1])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m2 #lowest DIC

plot(chosen)
summary(chosen)


###############################
# heterogeneity and variances #
###############################

# WI = weight, s2I = measurement error variance = sigma2m
WI <- na.omit(1/meta4$VZr)

s2I <- sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

total_var <- chosen$VCV[,"paperID"]+
  chosen$VCV[,"popID2"]+
  chosen$VCV[,"units"]+s2I # total variance (sigma2t)


total_het <- (total_var-s2I)/total_var # total heterogeneity I2
mean(total_het)
HPDinterval(total_het)


#studyID heterogeneity
studyID_het <- chosen$VCV[,"paperID"]/total_var
mean(studyID_het)
HPDinterval(studyID_het)


#populationID heterogeneity
popID_het <- chosen$VCV[,"popID2"]/total_var
mean(popID_het)
HPDinterval(popID_het)


### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])


#############################################################
# PUBLICATION BIAS
#############################################################

# Egger's regressions performed on model residuals and measurement
# errors. Evidence for publication bias if intercept is different
# from 0.

chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta4$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta4$Precision<-sqrt(1/meta4$VZr)

meta4$MAR<-meta4$Zr-meta4$Prediction # meta-analytic residual!

meta4$zMAR<-meta4$MAR*meta4$Precision


# Egger's regression
Egger.meta4<-MCMCglmm(zMAR~Precision,family="gaussian",
                      verbose=FALSE,data=meta4,
                      nitt=2000000, thin=1800, burnin=200000)

save(Egger.meta4,
     file="processed_data/metadatasets/model_meta4_paper&pop_Egger.full.RData")

load("processed_data/metadatasets/model_meta4_paper&pop_Egger.full.RData")


plot(Egger.meta4)
autocorr(Egger.meta4$Sol)
summary(Egger.meta4)
HPDinterval(Egger.meta4$Sol[,1])



###################################################################
# META 4: FULL -> META-REGRESSION: Biological moderators for meta 2
###################################################################

meta4$season.num <- ifelse(meta4$season=="nonbreeding",0,1)
meta4$int.num <- ifelse(meta4$interactions=="both",0,1)
meta4$inttype.num <- ifelse(meta4$interactiontype=="mix",0,1)


meta_reg4.1 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                    random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta_reg4.2 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                    random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)
meta_reg4.3 <- MCMCglmm(Zr~scale(season.num)+scale(int.num)+scale(inttype.num),
                    random=~paperID+popID2,
                    mev=meta4$VZr,
                    data=meta4,
                    family="gaussian",
                    nitt=2000000,thin=1800,burnin=200000,
                    pr=TRUE,verbose=FALSE,
                    prior=prior1)

save(meta_reg4.1,meta_reg4.2,meta_reg4.3,
     file="processed_data/metadatasets/model_meta_reg4_paper&pop.full_inttype.RData")

load("processed_data/metadatasets/model_meta_reg4_paper&pop.full_inttype.RData")

m1 <- meta_reg4.1
m2 <- meta_reg4.2
m3 <- meta_reg4.3


#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))

#check chain autocorrelation
autocorr(m1$Sol[,c(1,2,3)])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,c(1,2,3)])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,c(1,2,3)])
autocorr(m3$VCV[,c(1,2,4)])

#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))

chosen <- m1 #lowest DIC


plot(chosen)
summary(chosen)


############################
# Estimates and R2marginal #
############################

### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])

mean(chosen$Sol[,2])
HPDinterval(chosen$Sol[,2])

mean(chosen$Sol[,3])
HPDinterval(chosen$Sol[,3])

mean(chosen$Sol[,4])
HPDinterval(chosen$Sol[,4])


# # Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(chosen$Sol[i,c(1:4)] %*% t(chosen$X)))
  vmVarF[i]<-Var}
R2m<-100*(vmVarF/(vmVarF+chosen$VCV[,"paperID"]+chosen$VCV[,"popID2"]+chosen$VCV[,"units"]))
mean(R2m)
HPDinterval(R2m)



###################################################################
# META 5: RAW DATA -> META-REGRESSION: sampling effort as moderator
###################################################################

meta5 <- meta4[!(is.na(meta4$ratio)),]
meta5$ratio2 <- meta5$ratio^2

meta_reg5.1 <- MCMCglmm(Zr~scale(ratio)+scale(ratio2),
                        random=~paperID+popID2,
                        mev=meta5$VZr,
                        data=meta5,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg5.2 <- MCMCglmm(Zr~scale(ratio)+scale(ratio2),
                        random=~paperID+popID2,
                        mev=meta5$VZr,
                        data=meta5,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg5.3 <- MCMCglmm(Zr~scale(ratio)+scale(ratio2),
                        random=~paperID+popID2,
                        mev=meta5$VZr,
                        data=meta5,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)

save(meta_reg5.1,meta_reg5.2,meta_reg5.3,
     file="processed_data/metadatasets/model_meta_reg5_ratio_paper&pop.full.RData")


load("processed_data/metadatasets/model_meta_reg5_ratio_paper&pop.full.RData")

m1 <- meta_reg5.1
m2 <- meta_reg5.2
m3 <- meta_reg5.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,c(1,2,3)])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,c(1,2,3)])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,c(1,2,3)])
autocorr(m3$VCV[,c(1,2,4)])


#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))


chosen <- m1 #lowest DIC

plot(chosen)
summary(chosen)


############################
# Estimates and R2marginal #
############################

### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])

mean(chosen$Sol[,2])
HPDinterval(chosen$Sol[,2])

mean(chosen$Sol[,3])
HPDinterval(chosen$Sol[,3])


# # Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(chosen$Sol[i,c(1:3)] %*% t(chosen$X)))
  vmVarF[i]<-Var}
R2m<-100*(vmVarF/(vmVarF+chosen$VCV[,"paperID"]+chosen$VCV[,"popID2"]+chosen$VCV[,"units"]))
mean(R2m)
HPDinterval(R2m)



###############################################################
# META 6: RAW DATA -> META-REGRESSION: published vs unpublished
###############################################################

pub <- as.character(meta1$paperID)

meta3$pubvsunpub <- as.factor(ifelse(as.character(meta3$paperID) %in% pub,1,0))


meta_reg6.1 <- MCMCglmm(Zr~pubvsunpub,
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg6.2 <- MCMCglmm(Zr~pubvsunpub,
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)
meta_reg6.3 <- MCMCglmm(Zr~pubvsunpub,
                        random=~paperID+popID2,
                        mev=meta3$VZr,
                        data=meta3,
                        family="gaussian",
                        nitt=2000000,thin=1800,burnin=200000,
                        pr=TRUE,verbose=FALSE,
                        prior=prior1)

save(meta_reg6.1,meta_reg6.2,meta_reg6.3,
     file="processed_data/metadatasets/model_meta_reg6_pubvsunpub_paper&pop.full_0unp_v2.RData")


load("processed_data/metadatasets/model_meta_reg6_pubvsunpub_paper&pop.full_0unp_v2.RData")


m1 <- meta_reg6.1
m2 <- meta_reg6.2
m3 <- meta_reg6.3

#check model convergence
#The shrink factors should be below 1.05
gelman.diag(list(m1$Sol[,1],m2$Sol[,1],m3$Sol[,1]))
gelman.diag(list(m1$VCV[,c(1)],m2$VCV[,c(1)],m3$VCV[,c(1)]))
gelman.diag(list(m1$Deviance,m2$Deviance,m3$Deviance))


#check chain autocorrelation
autocorr(m1$Sol[,c(1,2)])
autocorr(m1$VCV[,c(1,2,4)])
autocorr(m2$Sol[,c(1,2)])
autocorr(m2$VCV[,c(1,2,4)])
autocorr(m3$Sol[,c(1,2)])
autocorr(m3$VCV[,c(1,2,4)])


#extract and summarise DIC
print(c(m1$DIC,m2$DIC,m3$DIC))


chosen <- m2 #lowest DIC

plot(chosen)
summary(chosen)


############################
# Estimates and R2marginal #
############################

### fixed effects (overall intercept)
mean(chosen$Sol[,1])
HPDinterval(chosen$Sol[,1])

mean(chosen$Sol[,2])
HPDinterval(chosen$Sol[,2])


# # Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(chosen$Sol[i,c(1:2)] %*% t(chosen$X)))
  vmVarF[i]<-Var}
R2m<-100*(vmVarF/(vmVarF+chosen$VCV[,"paperID"]+chosen$VCV[,"popID2"]+chosen$VCV[,"units"]))
mean(R2m)
HPDinterval(R2m)




##############################################################
# MANUAL FOREST PLOT: Figure 1 in main text
##############################################################
#loading models
load("processed_data/metadatasets/model_meta3_paper&pop.full.RData")
meta1 <- meta3.1

load("processed_data/metadatasets/model_meta4_paper&pop.full.RData")
meta2 <- meta4.2


tiff("Plots/Figure1_forest_plot_meta-analysis.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)

yaxis <- rev(c(1:2))
# means <- c(0.229,0.196)
# lowerCI <- c(0.008,-0.008)
# upperCI <- c(0.462,0.414)
means <- c(mean(meta1$Sol[,1]),mean(meta2$Sol[,1]))
lowerCI <- c(HPDinterval(meta1$Sol[,1])[1],
             HPDinterval(meta2$Sol[,1])[1])
upperCI <- c(HPDinterval(meta1$Sol[,1])[2],
             HPDinterval(meta2$Sol[,1])[2])
analysis <- c(rep(c("black"),2))

labels <- c("meta 1",
            "meta 2")

k <- c(85,87)


op <- par(mar = c(4.25,4,1,1)) #bottom, left, top, and right. 

plot(means,yaxis,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(0.6,2.4),
     xlim=c(-0.5,1),
     frame.plot=FALSE)

abline(v=0, lwd=1.5, lty=1)


axis(1,at=seq(-0.5,1,0.5),
     labels=c("-0.5","0","0.5","1"),
     cex.axis=1.25,tck=-0.02)

axis(2,
     at=rev(c(1:2)),
     labels=labels,
     cex.axis=1.25,las=2,tck=0,lty=0,line=-1)

title(xlab = "effect size (Zr)", line = 2.75, cex.lab=1.75)

polygon(c(c(0.1003,0.310),rev(c(0.1003,0.310))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(0.310,0.549),rev(c(0.310,0.549))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.16))

polygon(c(c(0.549,1),rev(c(0.549,1))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.25))

polygon(c(c(-0.1003,-0.310),rev(c(-0.1003,-0.310))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(-0.310,-0.5),rev(c(-0.310,-0.5))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.16))


points(means,yaxis,col=analysis,pch=19,cex=2)

arrows(lowerCI,yaxis,
       upperCI,yaxis,
       angle=90,code=3,col=analysis,
       length = 0,lwd=2.75)

text(0.8,yaxis,k,adj=1,cex=1.15)
text(0.77,2.25,"k",adj=1,cex=1.15)

dev.off()



##############################################################
# MANUAL FOREST PLOT: SUPPLEMENTS
##############################################################

#loading models
load("processed_data/metadatasets/model_meta1_paper&pop.full.RData")
pub1 <- meta1.3

load("processed_data/metadatasets/model_meta2_paper&pop.full.RData")
pub2 <- meta2.2


tiff("Plots/FigureS1_forest_plot_meta-analysis.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)

yaxis <- rev(c(1:2))
means <- c(mean(pub1$Sol[,1]),mean(pub2$Sol[,1]))
lowerCI <- c(HPDinterval(pub1$Sol[,1])[1],
             HPDinterval(pub2$Sol[,1])[1])
upperCI <- c(HPDinterval(pub1$Sol[,1])[2],
             HPDinterval(pub2$Sol[,1])[2])
analysis <- c(rep(c("black"),2))

labels <- c("published 1",
            "published 2")

k <- c(20,53)


op <- par(mar = c(4.25,5.1,1,1)) #bottom, left, top, and right. 

plot(means,yaxis,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(0.6,2.4),
     xlim=c(-0.5,1),
     frame.plot=FALSE)

abline(v=0, lwd=1.5, lty=1)


axis(1,at=seq(-0.5,1,0.5),
     labels=c("-0.5","0","0.5","1"),
     cex.axis=1.25,tck=-0.02)

axis(2,
     at=rev(c(1:2)),
     labels=labels,
     cex.axis=1.15,las=2,tck=0,lty=0,line=-1)

title(xlab = "effect size (Zr)", line = 2.75, cex.lab=1.75)

polygon(c(c(0.1003,0.310),rev(c(0.1003,0.310))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(0.310,0.549),rev(c(0.310,0.549))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.16))

polygon(c(c(0.549,1),rev(c(0.549,1))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.25))

polygon(c(c(-0.1003,-0.310),rev(c(-0.1003,-0.310))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(-0.310,-0.5),rev(c(-0.310,-0.5))),
        c(c(0,0),rev(c(13,13))),
        border=NA,col=rgb(0,0,0, 0.16))


points(means,yaxis,col=analysis,pch=19,cex=2)

arrows(lowerCI,yaxis,
       upperCI,yaxis,
       angle=90,code=3,col=analysis,
       length = 0,lwd=2.75)

text(0.9,yaxis,k,adj=1,cex=1.15)
text(0.87,2.25,"k",adj=1,cex=1.15)

dev.off()



##############################################################
# TIME-LAG BIAS PLOT
##############################################################

#first obtaining estimates from the model

load("processed_data/metadatasets/model_meta_reg2_year_paper&pop.full.RData")

chosen <- meta_reg2.2

newdat<-expand.grid(year = seq(1985,2016,0.01))

newdat$year.z <- scale(newdat$year)

xmat<-model.matrix(~year.z,data=newdat)


fitmatboth <- matrix(NA, 
                     ncol = nrow(chosen$Sol), 
                     nrow = nrow(newdat))


for(i in 1:nrow(chosen$Sol)) {
  fitmatboth[,i] <- xmat%*%chosen$Sol[i,c(1,2)]
}

newdat$fit<-apply(fitmatboth, 1, mean)
newdat$lower<-apply(fitmatboth, 1, quantile, prob= 0.025)
newdat$upper<-apply(fitmatboth, 1, quantile, prob= 0.975)


#############
# Actual plot

tiff("Plots/main/time-lag_bias_plot_Meta2.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=800)

xaxis <- meta2$year
yaxis <- meta2$Zr
cex.study <- meta2$N/8


op <- par(mar = c(4.5,4.5,1,1)) #bottom, left, top, and right. 

plot(xaxis,yaxis,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-2.5,2.5),
     xlim=c(1985,2016))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(1985,2016,10),
     cex.axis=1.25,tck=-0.02)

axis(2,
     at=seq(-2,2,1),
     cex.axis=1.25,las=2,tck=-0.02)

title(xlab = "year of publication", 
      ylab = "effect size (Zr)",
      line = 2.75, cex.lab=1.75)

polygon(c(newdat$year,rev(newdat$year)),
        c(newdat$lower,rev(newdat$upper)),
        border=NA,col=rgb(58/255,95/255,205/255, 0.25))

lines(newdat$year, newdat$fit, lwd=3.5,col="royalblue3") 


points(jitter(xaxis,2),yaxis,
       bg=rgb(58/255,95/255,205/255, 0.5),
       pch=21,
       cex=cex.study)

cex.legend <- c(4/8,10/8,25/8)

legend(1985,-1,
       c("n = 4 birds",
         "n = 10 birds",
         "n = 25 birds"),
       pt.bg=rgb(58/255,95/255,205/255, 0.75),
       pt.cex=cex.legend,
       pch=21,
       inset=c(0,0),
       y.intersp=1,x.intersp=1.3)


dev.off()



##############################################################
# PUBLISHED VS. UNPUBLISHED
##############################################################

#first obtaining estimates from the model
load("processed_data/metadatasets/model_meta_reg6_pubvsunpub_paper&pop.full_0unp_v2.RData")

chosen <- meta_reg6.2

summary(chosen)


##values
mean.unp <- mean(chosen$Sol[,1])
lower.unp <- HPDinterval(chosen$Sol[,1])[1]
upper.unp <- HPDinterval(chosen$Sol[,1])[2]

mean.pub <- mean(chosen$Sol[,2])
lower.pub <- HPDinterval(chosen$Sol[,2])[1]
upper.pub <- HPDinterval(chosen$Sol[,2])[2]


#############
# Actual plot

tiff("Plots/main/published_vs_unpublished.tiff",
     height=10.5, width=10.5,
     units='cm', compression="lzw", res=600)


op <- par(mar = c(5,5,1,1)) #bottom, left, top, and right. 

plot(c(0,1),newdat$fit,
     type="n",
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n",
     ylim=c(-1,1),
     xlim=c(-0.4,1.4))


abline(a=0,b=0, lwd=1, lty=1)


axis(1,at=seq(0,1,1),
     labels=c("published","unpublished"),
     cex.axis=1.25,tck=-0.02)

axis(2,
     at=seq(-1,1,0.5),
     cex.axis=1.25,las=2,tck=-0.02)

title(ylab = "effect size (Zr)",
      xlab = "primary data",
      line = 3.25, cex.lab=1.75)

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(0.1003,0.310),rev(c(0.1003,0.310))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(0.310,0.549),rev(c(0.310,0.549))),
        border=NA,col=rgb(0,0,0, 0.16))

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(0.549,1.5),rev(c(0.549,1.5))),
        border=NA,col=rgb(0,0,0, 0.25))

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(-0.1003,-0.310),rev(c(-0.1003,-0.310))),
        border=NA,col=rgb(0,0,0, 0.07))

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(-0.310,-0.549),rev(c(-0.310,-0.549))),
        border=NA,col=rgb(0,0,0, 0.16))

polygon(c(c(-1,-1),rev(c(2,2))),
        c(c(-0.549,-1.5),rev(c(-0.549,-1.5))),
        border=NA,col=rgb(0,0,0, 0.25))

arrows(1,lower.unp,
       1,upper.unp,
       angle=90,code=3,
       col=rgb(1,165/255,0,0.85),
       length = 0,lwd=3.5)

arrows(0,lower.pub,
       0,upper.pub,
       angle=90,code=3,
       col=rgb(58/255,95/255,205/255, 0.85),
       length = 0,lwd=3.5)

points(1,mean.unp,
       col=rgb(1,165/255,0,0.95),
       pch=19,
       cex=2)

points(0,mean.pub,
       col=rgb(58/255,95/255,205/255, 0.95),
       pch=19,
       cex=2)


text(0.05,0.96,"k = 53",adj=0.5,cex=1.15)
text(1.05,0.38,"k = 32",adj=0.5,cex=1.15)


dev.off()
