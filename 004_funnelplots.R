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

# This script is to create funnel plots based on meta-analytic
# residuals.

# Meta-analytic challenge to a textbook example of status
# signalling: publication trends and biases

# Sanchez-Tojar et al. In prep.

# more info at the Open Science Framework: osf.io/cwkxb


##############################################################
# Packages needed
##############################################################

# load pacakges
library(meta) 
library(metafor)
library(MCMCglmm)

# cleaning up
rm(list=ls())


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

tiff("plots/Figure2_Funnel_plots.tiff",
     #"plots/FigureS2_Funnel_plots.tiff",
     height=10.5, width=21,
     units='cm', compression="lzw", res=600)

par(mfrow=c(1,2))

op <- par(oma = c(3,3.5,0,1) + 0.1,
          mar = c(1,2.25,1,0.5) + 0.1)

# ##############################################################
# # META 1: FULL
# ##############################################################
# 
# load("processed_data/metadatasets/model_meta1_paper&pop.full.RData")
# 
# chosen <- meta1.2 #lowest DIC
# 
# chosen$Random$formula<-update(chosen$Random$formula,
#                               ~.+leg(mev, -1, FALSE):units)
# 
# # getting predictions
# # (raw data - predictions = meta-analytic residuals)
# meta1$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)
# 
# meta1$Precision<-sqrt(1/meta1$VZr)
# 
# meta1$MAR<-meta1$Zr-meta1$Prediction # meta-analytic residual!
# 
# 
# plot(meta1$MAR,1/meta1$VZr,
#      pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5,
#      yaxt="n",
#      xlim=c(-3,3),
#      ylim=c(0,60),
#      col=rgb(58/255,95/255,205/255, 0.5))
# 
# abline(v=0, lwd=1, lty=2)
# 
# text(0.8,56,"published 1",adj = 0 ,cex=1.25)
# 
# axis(2,at=seq(0,60,10),
#      cex.axis=1.25,tck=-0.02,las=2)
# 
# 
# ##############################################################
# # META 2: FULL
# ##############################################################
# 
# load("processed_data/metadatasets/model_meta2_paper&pop.full.RData")
# 
# chosen <- meta2.1 #lowest DIC
# 
# chosen$Random$formula<-update(chosen$Random$formula,
#                               ~.+leg(mev, -1, FALSE):units)
# 
# # getting predictions
# # (raw data - predictions = meta-analytic residuals)
# meta2$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)
# 
# meta2$Precision<-sqrt(1/meta2$VZr)
# 
# meta2$MAR<-meta2$Zr-meta2$Prediction # meta-analytic residual!
# 
# 
# plot(meta2$MAR,1/meta2$VZr,
#      pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5,
#      yaxt="n",
#      xlim=c(-3,3),
#      ylim=c(0,60),
#      col=rgb(58/255,95/255,205/255, 0.5))
# 
# abline(v=0, lwd=1, lty=2)
# 
# text(0.8,56,"published 2",adj = 0 ,cex=1.25)
# 
# axis(2,at=seq(0,60,10),
#      cex.axis=1.25,tck=-0.02,las=2)


##############################################################
# META 3: FULL
##############################################################

load("processed_data/metadatasets/model_meta3_paper&pop.full.RData")

chosen <- meta3.1 #lowest DIC

chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta3$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta3$Precision<-sqrt(1/meta3$VZr)

meta3$MAR<-meta3$Zr-meta3$Prediction # meta-analytic residual!


# for differential coloring
meta3$coloring <- ifelse(meta3$published=="yes",
                         rgb(58/255,95/255,205/255, 0.5),
                         rgb(1,165/255,0,0.5))

plot(meta3$MAR,1/meta3$VZr,
     pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5,
     yaxt="n",
     xlim=c(-3,3),
     ylim=c(0,45),
     col=meta3$coloring)


abline(v=0, lwd=1, lty=2)

text(1.3,42,"meta 1",adj = 0 ,cex=1.5)


axis(2,at=seq(0,40,10),
     cex.axis=1.25,tck=-0.02,las=2)



##############################################################
# META 4: FULL
##############################################################

load("processed_data/metadatasets/model_meta4_paper&pop.full.RData")

chosen <- meta4.2 #lowest DIC


chosen$Random$formula<-update(chosen$Random$formula,
                              ~.+leg(mev, -1, FALSE):units)

# getting predictions
# (raw data - predictions = meta-analytic residuals)
meta4$Prediction<-predict(chosen, marginal=~leg(mev, -1, FALSE):units)

meta4$Precision<-sqrt(1/meta4$VZr)

meta4$MAR<-meta4$Zr-meta4$Prediction # meta-analytic residual!


# for differential coloring

nonreported <- c("2010_xb_Dolnik","1991_x_Andersson")

meta4$coloring <- ifelse(meta4$published=="yes",
                         rgb(58/255,95/255,205/255, 0.5),
                         ifelse((meta4$study %in% nonreported),
                                rgb(0,0,0, 0.5),
                                rgb(1,165/255,0,0.5)))

plot(meta4$MAR,1/meta4$VZr,
     pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5,
     yaxt="n",
     xlim=c(-3,3),
     ylim=c(0,45),
     col=meta4$coloring)

abline(v=0, lwd=1, lty=2)

text(1.3,42,"meta 2",adj = 0 ,cex=1.5)

axis(2,at=seq(0,40,10),
     cex.axis=1.25,tck=-0.02,las=2)


legend(-3,45,
       c("published",
         "unpublished",
         "non-reported"),
       pt.bg="white",
       col=c(rgb(58/255,95/255,205/255, 0.5),
             rgb(1,165/255,0,0.5),
             rgb(0,0,0, 0.5)),
       pt.cex=1,
       pch=c(rep(19,3)),
       inset=c(0,0),
       y.intersp=1,x.intersp=1.3,
       cex=0.6)

title(ylab = "precision (1 / VZr)",
      xlab = "meta-analytic residuals (Zr)",
      outer = TRUE, line = 1.75, cex.lab=1.75)

par(mfrow=c(1,1))

dev.off()
