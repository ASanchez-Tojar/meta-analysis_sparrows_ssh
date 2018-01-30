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

# This script is to analyze raw interaction and bib size data
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

library(aniDom)
library(compete)


# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

#adaptation from the "aniDom" v.0.1.2 function, plot_hierarchy_shape()
plot_hierarchy_shape_2 <- function (identity, rank, winners, losers, fitted = FALSE) 
{
  winners.rank <- rank[match(winners, identity)]
  losers.rank <- rank[match(losers, identity)]
  xx <- winners.rank - losers.rank
  x <- 1:(max(abs(xx)))
  y <- rep(NA, length(x))
  CI.upper <- y
  CI.lower <- y
  for (i in 1:length(x)) {
    y[i] <- sum(xx == -x[i])/sum(abs(xx) == x[i])
    CI.upper[i] <- y[i] + sqrt(y[i] * (1 - y[i])/sum(abs(xx) ==
                                                       x[i])) + 0.5/sum(abs(xx) == x[i])
    CI.upper[i] <- min(CI.upper[i], 1)
    CI.lower[i] <- y[i] - sqrt(y[i] * (1 - y[i])/sum(abs(xx) ==
                                                       x[i])) - 0.5/sum(abs(xx) == x[i])
    CI.lower[i] <- max(CI.lower[i], 0)
    
  }
  x <- x[!is.na(y)]
  CI.upper <- CI.upper[!is.na(y)]
  CI.lower <- CI.lower[!is.na(y)]
  y <- y[!is.na(y)]
  plot(x, y, xlab = "Difference in rank", ylab = "Probability that higher rank wins", 
       ylim = c(0, 1), pch = 20, cex = 2)
  arrows(x, CI.lower, x, CI.upper, length = 0.1, angle = 90, 
         code = 3)
  if (fitted) {
    l <- loess(y ~ x)
    lines(x, l$fitted, col = "red", lwd = 2)
  }
  invisible(data.frame(Rank.diff = x, Prob.dom.win = y, CI.upper = CI.upper, 
                       CI.lower = CI.lower))
}


# Estimating bib area according to Moller 1987
# This is an equation estimated from analysing 29 museum skins, 
# for which hidden bib was measured. Hidden or true bib is: the
# area covered with feathers having black bases.
# This function has been widely used, also for estimating 
# visible bib.

bib_Moller1987 <- function(length,width){ #length and width in mm
  
  area <- 166.7 + 0.45 * length * width

}

# Estimating bib area according to Veiga 1993
# This is based on the equation to estimate the area of 
# a circular serctor knowing the radius and the chord of it.
# The specific equations were obtained from: 
# http://mathworld.wolfram.com/CircularSegment.html
# Though this equation appears to be a much more standardized
# way of estimating bib area, it hasn't been used as much as 
# Moller's.

bib_Veiga1993 <- function(length,width){ #length and width in mm
  
  radius <- length
  chord <- width
  angle <- 1/(2*sin(chord/(2*radius)))
  sectorarea <- 0.5*width^2*(angle-sin(angle))
  area<-2*sectorarea

}


# elo_scores() from "aniDom" v.0.1.2 improved (typo fixed)

elo_scores_2 <- function (winners, losers, identities = NULL, sigmoid.param = 1/100, 
                            K = 200, init.score = 0, randomise = FALSE, n.rands = 1000, 
                            return.as.ranks = FALSE, return.trajectories = FALSE) 
{
  if (is.null(identities)) {
    identities <- unique(c(winners, losers))
  }
  n.inds <- length(identities)
  if (randomise == FALSE) {
    n.rands <- 1
  }
  if (sum(c(winners, losers) %in% identities) < length(c(winners, 
                                                         losers))) {
    stop("Not all winners and/or losers are contained in identities")
  }
  T <- length(winners)
  if (return.trajectories) {
    if (randomise == TRUE) {
      all.scores <- array(init.score, c(n.inds, T + 1, 
                                        n.rands))
    }
    else {
      all.scores <- array(init.score, c(n.inds, T + 1))
    }
  }
  else {
    all.scores <- array(init.score, c(n.inds, n.rands))
  }
  rownames(all.scores) <- identities
  if (length(K) == 1) {
    K <- rep(K, T)
  }
  for (r in 1:n.rands) {
    if (randomise == FALSE) {
      ord <- 1:T
    }
    else {
      ord <- sample(1:T, T, replace = F)
    }
    winners.perm <- winners[ord]
    losers.perm <- losers[ord]
    scores <- array(NA, c(n.inds, T + 1))
    scores[, 1] <- init.score
    for (i in 1:T) {
      scores[, i + 1] <- scores[, i]
      winner <- which(identities == winners.perm[i])
      loser <- which(identities == losers.perm[i])
      p <- 1/(1 + exp(-sigmoid.param * abs((scores[winner, 
                                                   i] - scores[loser, i]))))
      if (scores[winner, i] >= scores[loser, i]) {
        scores[winner, i + 1] <- scores[winner, i] + 
          (1 - p) * K[i]
        scores[loser, i + 1] <- scores[loser, i] - (1 - 
                                                      p) * K[i]
      }
      else {
        scores[winner, i + 1] <- scores[winner, i] + 
          p * K[i]
        scores[loser, i + 1] <- scores[loser, i] - 
          p * K[i]
        }
    }
    if (return.trajectories) {
      if (randomise == TRUE) {
        all.scores[, , r] <- scores
      }
      else {
        all.scores <- scores
      }
    }
    else {
      all.scores[, r] <- scores[, T + 1]
    }
  }
  freq <- table(factor(c(winners, losers), levels = identities))
  all.scores[which(identities %in% names(freq)[which(freq == 
                                                       0)]), ] <- NA
  if (return.as.ranks == TRUE) {
    all.ranks <- apply(all.scores, 2, function(x) {
      rank(-x)
    })
    all.ranks[is.na(all.scores)] <- NA
    all.scores <- all.ranks
  }
  invisible(all.scores)
}


#############################################################
# Analysis of all raw data
#############################################################

# get metadata: the index contains a list of all accessible
# matrices and sequences
indexdb <- read.csv("001_data_index.csv",
                   stringsAsFactors=FALSE, header=TRUE)


cor.S.list <- c()
cor.S.p.list <- c()
cor.K.list <- c()
cor.K.p.list <- c()
identifier.list <- c()
Ttri.list <- c()
Ttri.p.list <- c()
ratio.list <- c()
unknowndyads.list <- c()
N.list <- c()
bib.type.list <- c()
rept.list <- c()
year.list <- c()
groupID.list <- c()
author.list <- c()
x.mean.length <- c()
x.sd.length <- c()
x.mean.width <- c()
x.sd.width <- c()

# loop through each dataset
for (z in 1:nrow(indexdb)) {
  
  year <- indexdb$year[z]
  group_ID <- indexdb$group_ID[z]
  author_ID <- indexdb$authorID[z]
  
  identifier <- paste(year,group_ID,author_ID,sep="_")
  
  if(indexdb$original_format[z]=="matrix"){
    
    mat.identifier <- paste0("matrix_",identifier)
    
    # reading each group
    mat <- read.csv(paste0(mat.identifier,".csv"),
                    header=TRUE, row.names=1)
    
    # Hein's matrices are coded unusually, with losers displayed by row
    # and winner by column. Therefore, these matrices need to be 
    # transpose using t(). The other matrices are correct.
    if(author_ID=="Hein"){
      
      mat <- t(mat)
      
    }
    
    intperind <- rowSums(mat)+colSums(mat)
    
    mat.red <- subset(mat,intperind!=0,intperind!=0)
    
    mat2 <- mat.red
    ids <- rownames(mat2)
    
    
    # blank dataset
    dom.data <- data.frame(interact.number=1:sum(mat2),
                           winner=NA,loser=NA)
    
    
    # fill it in
    count <- 1
    for (i in 1:nrow(mat2)) {
      for (j in 1:ncol(mat2)) {
        while (mat2[i,j] > 0) {
          dom.data$winner[count] <- ids[i]
          dom.data$loser[count] <- ids[j]
          mat2[i,j] <- mat2[i,j]-1
          count <- count + 1
        }
      }
    }
    
    
    # estimating randomized Elo-ratings
    scores <- elo_scores_2(winners=dom.data$winner,
                           losers=dom.data$loser,
                           identities = row.names(mat2),
                           randomise = TRUE,
                           n.rands = 1000,
                           return.as.ranks = TRUE)


    # ratio of interactions to individuals informs about
    # the uncertainty of the inferred hierarchy
    ratio <- length(dom.data$winner)/length(ids)


    # rank of the ratings
    rank <- rank(rowMeans(scores))


    # sparseness also informs about the uncertainty of
    # the inferred hierarchy
    unknowndyads<-rshps(mat.red)$unknowns/rshps(mat.red)$total


    # Triangle transitivity and its "significance" test
    # Are the dominance relationships transitive?
    Ttri <- ttri_test(mat.red)$ttri
    Ttri.p <- ttri_test(mat.red)$pval


    # Estimate of uncertainty based on "aniDom"
    rept <- estimate_uncertainty_by_repeatability(winners=dom.data$winner,
                                                  losers=dom.data$loser,
                                                  identities = row.names(mat2))

    # Plotting the shape of the hierarchy
    tiff(paste("Plots/hierarchies/",
               identifier,
               ".tiff",sep=""),
         height=10, width=10,
         units='cm', compression="lzw", res=600)

    plot_hierarchy_shape_2(fitted=TRUE,
                           ids,rank,
                           dom.data$winner,
                           dom.data$loser)

    lines(c(0,length(rank)),c(0.5,0.5),col="red",
          lty=3,lwd=1.5)

    text(1,0.02,paste0("ratio = ",round(ratio,1)),
         adj = 0,cex=0.8)

    text(1,0.12,paste0("sparseness = ",
                      round(unknowndyads,2)),adj = 0,cex=0.8)

    text(1,0.22,paste0("Ttri = ",
                      round(Ttri,2)),adj = 0,cex=0.8)

    dev.off()
    
    
  } else {
    
    seq.identifier <- paste0("sequence_",identifier)
    
    seq <- read.csv(paste0(seq.identifier,".csv"),header=TRUE, sep=",")
    seq$winner <- as.character(seq$winner)
    seq$loser <- as.character(seq$loser)
    
    
    # first, let's get rid off a special case of typo (in case it existed): 
    # winner = loser. Interactions involving so, will be deleted
    # from here on.
    seq <- seq[seq$winner!=seq$loser,]
    
    
    ids <- unique(c(as.character(seq$winner),as.character(seq$loser)))

    scores.seq <- elo_scores(winners=as.character(seq$winner),
                             losers=as.character(seq$loser),
                             identities = ids,
                             randomise = TRUE,
                             n.rands = 1000,
                             return.as.ranks = TRUE)

    ratio <- nrow(seq)/length(ids)


    rank <- rank(rowMeans(scores.seq))


    matrix.look <- get_wl_matrix(seq[,c("winner","loser")])


    unknowndyads<-rshps(matrix.look)$unknowns/rshps(matrix.look)$total


    Ttri <- ttri_test(matrix.look)$ttri
    Ttri.p <- ttri_test(matrix.look)$pval


    rept <- estimate_uncertainty_by_repeatability(winners=as.character(seq$winner),
                                                  losers=as.character(seq$loser),
                                                  identities = ids)

    # Plotting the shape of the hierarchy
    tiff(paste("Plots/hierarchies/",
               identifier,
               ".tiff",sep=""),
         height=10, width=10,
         units='cm', compression="lzw", res=600)

    ids <- unique(c(as.character(seq$winner),as.character(seq$loser)))

    plot_hierarchy_shape_2(fitted=TRUE,
                           ids,rank,
                           seq$winner,
                           seq$loser)

    lines(c(0,length(rank)),c(0.5,0.5),
          col="red",lty=3,lwd=1.5)

    text(1,0.02,paste0("ratio = ",
                       round(ratio,1)),adj = 0,cex=0.8)

    text(1,0.12,paste0("sparseness = ",
                       round(unknowndyads,2)),adj = 0,cex=0.8)

    text(1,0.22,paste0("Ttri = ",
                       round(Ttri,2)),adj = 0,cex=0.8)


    dev.off()
    
  }
  
  dom.db <- as.data.frame(cbind(names(rank),rank),
                          row.names=c(1,length(rank)))
  
  names(dom.db) <- c("ID","rank")

  bib.identifier <- paste0("bib_",identifier)

  bib.db <- read.csv(paste0(bib.identifier,".csv"),
                     header=TRUE, sep=",")

  bib.db.2 <- bib.db[bib.db$sex=="m",]
  bib.db.3 <- bib.db.2[!(is.na(bib.db.2$bib)),]

  bib.db.3$bib_type <- factor(bib.db.3$bib_type)

  if("width" %in% levels(bib.db.3$bib_type)){

    length.ID <- bib.db.3[bib.db.3$bib_type=="length",c("ID","bib")]
    length.ID <- length.ID[order(length.ID$ID),] #making sure everything is in the same order
    length <- length.ID[,c("bib")]
    
    width.ID <- bib.db.3[bib.db.3$bib_type=="width",c("ID","bib")]
    width.ID <- width.ID[order(width.ID$ID),]
    width <- width.ID[,c("bib")]

    #only if all IDs match, calculate area, otherwise,tell me!
    if(!("FALSE" %in% names(table(as.character(length.ID[,c("ID")])==
                                 as.character(width.ID[,c("ID")]))))){
      
      visible <- bib_Moller1987(length,width)
      
    } else{
      
      print("WARNING! Fucked up data")
      break
      
    }
    
    

    temp.db <- bib.db.3[bib.db.3$bib_type=="width",
                        c("ID","sex")]

    db.bibarea <- as.data.frame(cbind(as.character(temp.db$ID),
                                      as.character(temp.db$sex),
                                      as.numeric(visible),
                                      rep("mm2",length(visible)),
                                      rep("visible",length(visible))))

    names(db.bibarea) <- colnames(bib.db.3)

    bib.db.3 <- rbind(bib.db.3,db.bibarea)

    bib.db.3$bib <- as.numeric(bib.db.3$bib)

  }

  bibtype <- 0

  for (bibtype in levels(bib.db.3$bib_type)){

    if(bibtype!="width"){

      bib.db.4 <- bib.db.3[bib.db.3$bib_type==bibtype,]

      dom.bib.db <- merge(bib.db.4,dom.db,
                          by="ID",all.x=TRUE)

      dom.bib.db$rank <- as.integer(as.character(dom.bib.db$rank))

      dom.bib.db.2 <- dom.bib.db[!(is.na(dom.bib.db$bib)),]

      cor.S <- cor.test(dom.bib.db.2$bib,-dom.bib.db.2$rank,
                      method=c("spearman"),exact=F) #exact=F is to get rid off the warning:
                                                    #Cannot compute exact p-value with ties
                                                    #We know, but we are happy with approximate
                                                    #p-values in this case. No drama.

      cor.K <- cor.test(dom.bib.db.2$bib,-dom.bib.db.2$rank,
                        method=c("kendall"),exact=F) #same comment as above

      N <- nrow(dom.bib.db.2[!(is.na(dom.bib.db.2$rank)),])

      cor.S.list <- c(cor.S.list,cor.S$estimate[[1]])
      cor.S.p.list <- c(cor.S.p.list,round(cor.S$p.value,3))
      cor.K.list <- c(cor.K.list,cor.K$estimate[[1]])
      cor.K.p.list <- c(cor.K.p.list,round(cor.K$p.value,3))
      identifier.list <- c(identifier.list,identifier)
      Ttri.list <- c(Ttri.list,round(Ttri,2))
      Ttri.p.list <- c(Ttri.p.list,round(Ttri.p,3))
      ratio.list <- c(ratio.list,round(ratio,1))
      unknowndyads.list <- c(unknowndyads.list,round(unknowndyads,2))
      N.list <- c(N.list,N)
      bib.type.list <- c(bib.type.list,as.character(bibtype))
      rept.list <- c(rept.list,round(rept,3))
      year.list <- c(year.list,year)
      groupID.list <- c(groupID.list,group_ID)
      author.list <- c(author.list,author_ID)

      if(bibtype=="length"){

        x.mean.length <- c(x.mean.length,
                           mean(bib.db.3[bib.db.3$bib_type=="length",
                                         c("bib")]))
        x.sd.length <- c(x.sd.length,
                         sd(bib.db.3[bib.db.3$bib_type=="length",
                                     c("bib")]))

      }

    }

    else{

      x.mean.width <- c(x.mean.width,
                        mean(bib.db.3[bib.db.3$bib_type=="width",
                                      c("bib")]))
      x.sd.width <- c(x.sd.width,
                      sd(bib.db.3[bib.db.3$bib_type=="width",
                                  c("bib")]))

    }

  }
  
}

studies.dom <- as.data.frame(cbind(identifier.list,
                                   author.list,
                                   groupID.list,
                                   year.list,
                                   Ttri.list,
                                   Ttri.p.list,
                                   rept.list,
                                   ratio.list,
                                   unknowndyads.list,
                                   N.list,
                                   cor.S.list,
                                   cor.S.p.list,
                                   cor.K.list,
                                   cor.K.p.list,
                                   bib.type.list))

names(studies.dom) <- c("study","authorID","groupID","year",
                        "Ttri","Ttri_pvalue",
                        "repeatability.sim","ratio",
                        "prop.unknowndyads",
                        "N","r_spearman","r.pvalue_spearman",
                        "r_kendall","r.pvalue_kendall",
                        "bib_type")
                                

studies.dom$year <- as.integer(as.character(studies.dom$year))
studies.dom$Ttri <- as.numeric(as.character(studies.dom$Ttri))
studies.dom$Ttri_pvalue <- as.numeric(as.character(studies.dom$Ttri_pvalue))
studies.dom$ratio <- as.numeric(as.character(studies.dom$ratio))
studies.dom$prop.unknowndyads <- as.numeric(as.character(studies.dom$prop.unknowndyads))
studies.dom$N <- as.numeric(as.character(studies.dom$N))
studies.dom$repeatability.sim <- as.numeric(as.character(studies.dom$repeatability.sim))
studies.dom$r_spearman <- as.numeric(as.character(studies.dom$r_spearman))
studies.dom$r_kendall <- as.numeric(as.character(studies.dom$r_kendall))
studies.dom$r.pvalue_spearman <- as.numeric(as.character(studies.dom$r.pvalue_spearman))
studies.dom$r.pvalue_kendall <- as.numeric(as.character(studies.dom$r.pvalue_kendall))


# Creating a couple of variables to ease subsetting data later on
# Classifying all these estimates as our calculations from raw data
studies.dom$originalestimate <- "no"
studies.dom$originalestimate <- as.factor(studies.dom$originalestimate)

# original data published before?

unpub <- c("Lendvai","Toth","Sanchez-Tojar","Westneat")
studies.dom$published <- ifelse(as.character(studies.dom$authorID) %in% unpub,
                                "no",
                                ifelse(as.character(studies.dom$authorID)=="Bokony" &
                                         studies.dom$year==2010,
                                       "no",
                                       "yes"))
studies.dom$published <- as.factor(studies.dom$published)

# rawdata?
studies.dom$rawdata <- "yes"
studies.dom$rawdata <- as.factor(studies.dom$rawdata)


# reducing dataset to one type of bib measurement per study (see main text)
studies.dom.red <- studies.dom[studies.dom$bib_type!="visible2",]


new.data <- data.frame()
data.temp.2 <- data.frame()

for (study in levels(studies.dom.red$authorID)){
  
  data.temp <- studies.dom.red[studies.dom.red$authorID==study,]
  
  if(study=="Rojas"){
    
    data.temp.2 <- data.temp[data.temp$bib_type!="hidden",]
    
  } else if(study=="Hein"){
    
    data.temp.2 <- data.temp[data.temp$bib_type=="drawing" |
                               data.temp$bib_type=="hidden",]
    
  } else {
    
    if(study!="Sanchez-Tojar"){
    
    data.temp.2 <- data.temp[data.temp$bib_type!="length",]
    
    } else {
      
      data.temp.2 <- data.temp
      
    }
  }
  
  new.data <- rbind(new.data,data.temp.2)
  
}

# removing length measurements from Sanchez-Tojar 2014 (i.e. area was preferred)
new.data <- new.data[!(new.data$year==2014 &
                         new.data$bib_type=="length"),]

new.data[,c("study","N")]


# saving data
write.csv(new.data[order(new.data$year,new.data$groupID,new.data$authorID),],
          "processed_data/re-analysis_rElo_and_prioritybib.csv",row.names=FALSE)


new.data <- read.table("processed_data/re-analysis_rElo_and_prioritybib.csv",
                       header=TRUE,sep=",")



# Plotting histograms for each "dominance index"

new.data <- new.data[new.data$N>3,]

tiff("Plots/hist_dominance_indices.tiff",
     height=27, width=18,units='cm', compression="lzw", res=600) 

par(mfrow=c(3,2))

variables <- c("Ttri","repeatability.sim","ratio",
               "prop.unknowndyads","N","r_spearman")

for (varia in levels(as.factor(variables))){
  
  histogramming <- new.data[,varia]
  
  hist(histogramming,
       main="",col="grey85",breaks=15,xlab=varia,cex.lab=1.6)
  
  lines(c(mean(histogramming),mean(histogramming)),
        c(0,150),col="black",lty=3,lwd=2.5)
  
  lines(c(median(histogramming),median(histogramming)),
        c(0,150),col="red",lty=3,lwd=2.5)
  
}

dev.off()


par(mfrow=c(1,1))