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

# This script is to build up the datasets that will be used
# to test whether bib size correlates with dominance rank in
# house sparrows across populations, research groups, etc...
# Part of the data comes from previously published studies, the 
# other part correspond to unpublished data obtained from
# researchers that previously published on the topic.

# Meta-analytic challenge to a textbook example of status
# signalling: publication trends and biases

# Sanchez-Tojar et al. In prep.

# more info at the Open Science Framework: osf.io/cwkxb


##############################################################
# Packages needed
##############################################################

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# Calculating Pearson's r from Spearman's rho
# Following Lajeunesse 2013 (Chapter 13)
spearman_to_pearson <- function(rho){
  r <- 2*sin((pi*rho)/6)
}


# Calculating Pearson's r from Kendall's tau 
# Following Lajeunesse 2013 (Chapter 13)
kendall_to_pearson <- function(tau){
  r <- sin((pi*tau)/2)
}

# function to convert r to Zr
r.to.Zr<-function(r){
  Zr<-round(0.5*(log(1+r)-log(1-r)),3)
}


# function to obtain variance
VZr <- function(N){
  VZr<-1/(N-3)
}


##############################################################
# Data
##############################################################

columnstokeep <- c("study","authorID","groupID","year",
                   "Ttri","Ttri_pvalue","repeatability.sim",
                   "ratio","prop.unknowndyads","N","bib_type",
                   "originalestimate","published","rawdata",
                   "typeofestimate","r")

# Original estimates from published data
original <- read.table("processed_data/norawdata_estimates.csv",
                       header=TRUE,sep=",")

# The following code assigns each estimate its corresponding type of estimate
study <- c("Ritchison","Moller","Solberg",
           "Liker","Gonzalez","Hein",
           "Riters","Lindstrom","Bokony",
           "Buchanan","Rojas","Dolnik")
typeofestimate <- c("chi-squared","Spearman","Spearman",
                    "Spearman","ANCOVA","Kendall",
                    "Spearman","Spearman","partial.regression",
                    "GLM","LMM","Pearson")

type.estimate.db <- data.frame(authorID=study,
                               typeofestimate=typeofestimate,
                               stringsAsFactors=FALSE)

original.2 <- merge(original,type.estimate.db,by="authorID",all.x=TRUE)
original.2$r <- ifelse(original.2$typeofestimate=="Spearman",
                       spearman_to_pearson(original.2$r_spearman),
                       ifelse(original.2$typeofestimate=="Kendall",
                              kendall_to_pearson(original.2$r_spearman),
                              original.2$r_spearman))

# r = 1 and r = -1, are transformed to 0.975 and -0.975
# to be able to estimate Zr (due to log(0))
original.2$r <- ifelse(round(original.2$r,9)==1,
                       0.975,
                       ifelse(round(original.2$r,9)==-1,
                              -0.975,
                              original.2$r))

original.3 <- original.2[,columnstokeep]


# Published re-analyzed data and unpublished analyzed data
reanalyzed <- read.table("processed_data/re-analysis_rElo_and_prioritybib.csv",
                         header=TRUE,sep=",")

reanalyzed$typeofestimate <- "Spearman"
reanalyzed$r <- spearman_to_pearson(reanalyzed$r_spearman)

# r = 1 and r = -1, are transformed to 0.975 and -0.975
# to be able to estimate Zr (due to log(0))
reanalyzed$r <- ifelse(round(reanalyzed$r,9)==1,
                       0.975,
                       ifelse(round(reanalyzed$r,9)==-1,
                              -0.975,
                              reanalyzed$r))

reanalyzed.2 <- reanalyzed[,columnstokeep]


# Non-reported non-significant estimates assigned to r = 0
nonreported <- read.table("processed_data/non-reported_estimates.csv",
                          header=TRUE,sep=",")

nonreported$typeofestimate <- "Pearson"
nonreported$r <- nonreported$r_spearman

nonreported.2 <- nonreported[,columnstokeep]


# Reported and non-reported estimates reanalized from Buchanan's processed data
Buchanan <- read.table("processed_data/re-analysis_Buchanan2010.csv",
                       header=TRUE,sep=",")

Buchanan$typeofestimate <- "Spearman"
Buchanan$r <- spearman_to_pearson(Buchanan$r_spearman)

# r = 1 and r = -1, are transformed to 0.975 and -0.975
# to be able to estimate Zr (due to log(0))
Buchanan$r <- ifelse(round(Buchanan$r,9)==1,
                       0.975,
                       ifelse(round(Buchanan$r,9)==-1,
                              -0.975,
                              Buchanan$r))

Buchanan.2 <- Buchanan[,columnstokeep]


# adding: type of interactions (aggressive vs non-aggressive) as moderator
typeint <- read.table("processed_data/interactiontype.csv",header=TRUE,sep=",")


##############################################################
# DATASETS
##############################################################

# 1. meta-analysis including only the original estimates of published studies
# Supplementary information
Meta1 <- original.3
Meta1$Zr <- r.to.Zr(Meta1$r)
Meta1$VZr <- VZr(Meta1$N)
Meta1$id <- seq(1:nrow(Meta1)) #same as groupID2


# Adding random effects from different database
random.meta1 <- read.table("processed_data/meta1_random_effects_and_exp.csv",
                           header=TRUE,sep=",")

Meta1.random <- merge(Meta1,random.meta1,by="authorID",all.x=TRUE)

Meta1.random <- merge(Meta1.random,
                      typeint,
                      by="paperID",
                      all.x=TRUE)

write.csv(Meta1.random[Meta1.random$N>3,],
          "processed_data/metadatasets/Meta1.csv",row.names=FALSE)


# 2. meta-analysis including only published studies, priority given
# to the standardized new estimates, then the original ones
# Supplementary Information
Meta2 <- rbind(original.3[original.3$rawdata=="no",],
               reanalyzed.2[reanalyzed.2$published=="yes",],
               Buchanan.2)
Meta2$Zr <- r.to.Zr(Meta2$r)
Meta2$VZr <- VZr(Meta2$N)
Meta2$id <- seq(1:nrow(Meta2))


# for the re-analyzed ones
random.meta2 <- read.table("processed_data/meta2-4_random_and_fixed.csv",
                           header=TRUE,sep=",")

random.meta2_4 <- random.meta2[,c("study","popID","popID2","groupID2","paperID",
                                "experimental","location","season","interactions")]

Meta2.random_and_fixed <- merge(Meta2,random.meta2_4,
                                by="study",all.x=TRUE)

Meta2.random_and_fixed <- merge(Meta2.random_and_fixed,
                                typeint,
                                by="paperID",
                                all.x=TRUE)

write.csv(Meta2.random_and_fixed[Meta2.random_and_fixed$N>3,],
          "processed_data/metadatasets/Meta2.csv",row.names=FALSE)


# # double checking
# Meta2.random_and_fixed[,c("study","id","popID","popID2","groupID2",
#                           "paperID","experimental",
#                           "location","season","interactions")]


# 3. meta-analysis including published and unpublished studies, priority
# given to the standardized new estimates, then the original ones.
# meta1 in main text
Meta3 <- rbind(original.3[original.3$rawdata=="no",],
               reanalyzed.2,
               Buchanan.2)
Meta3$Zr <- r.to.Zr(Meta3$r)
Meta3$VZr <- VZr(Meta3$N)
Meta3$id <- seq(1:nrow(Meta3))


Meta3.random_and_fixed <- merge(Meta3,random.meta2_4,
                                by="study",all.x=TRUE)

# # double checking
# Meta3.random_and_fixed[,c("study","id","popID","popID2","groupID2",
#                           "paperID","experimental",
#                           "location","season","interactions")]

Meta3.random_and_fixed <- merge(Meta3.random_and_fixed,
                                typeint,
                                by="paperID",
                                all.x=TRUE)

write.csv(Meta3.random_and_fixed[Meta3.random_and_fixed$N>3,],
          "processed_data/metadatasets/Meta3.csv",row.names=FALSE)

# print the number of estimates discarded
nrow(Meta3.random_and_fixed)-nrow(Meta3.random_and_fixed[Meta3.random_and_fixed$N>3,])

# 4. meta-analysis including published and unpublished studies, priority
# given to the standardized new estimates, then the original ones.
# Additionally, this dataset includes an approximation for the 2 
# non-significant estimates not reported. The assumption is that those 
# 2 estimates are highly non-significant, and thus, the estimates are 
# considered to be 0.
# meta2 in main text

# Sample sizes were unclear, so we used the number that were
# presented in the materials and methods. Specifically:
#
# Andersson & Ahlund 1991: We assumed the number of males described in the 
#                          materials and methods, that is, 20 males (10
#                          adult males and 10 yearling males)
#
# Dolnik & Hoi 2010: We assumed that the number of males before the 
#                    experiment was the same as the one after the 
#                    experiment, i.e. 31 (one bird was excluded due to
#                    health issues, see section "Experiment scheme").

Meta4 <- rbind(original.3[original.3$rawdata=="no",],
               reanalyzed.2,
               Buchanan.2,
               nonreported.2)

Meta4$Zr <- r.to.Zr(Meta4$r)
Meta4$VZr <- VZr(Meta4$N)
Meta4$id <- seq(1:nrow(Meta4))


Meta4.random_and_fixed <- merge(Meta4,random.meta2_4,
                                by="study",all.x=TRUE)

# #double checking
# Meta4.random_and_fixed[,c("study","id","popID","popID2","groupID2",
#                           "paperID","experimental",
#                           "location","season","interactions")]

Meta4.random_and_fixed <- merge(Meta4.random_and_fixed,
                                typeint,
                                by="paperID",
                                all.x=TRUE)

write.csv(Meta4.random_and_fixed[Meta4.random_and_fixed$N>3,],
          "processed_data/metadatasets/Meta4.csv",row.names=FALSE)

