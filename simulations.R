require(trendsegmentR)
require(IDetect)
require(not)
require(genlasso)
require(xtable)
require(dplyr)
require(cpop)

require(nsp)
require(dSTEM)


library(trendsegmentR)
library(IDetect)
library(genlasso)
library(not)



#######################################
#   Models
#######################################

model.wave1 <-  list(name = "wave1",
  sgnl.type = "PWLC",
  cpt = (1:9) * 150,
  chg.size = (-1)^{1:9} / 25,
  n = 150* 10,
  initial=c(-1,1/50))



model.wave2 <-  list(name = "wave2",
  sgnl.type = "PWL",
  cpt = (1:20) * 60,
  chg.size = matrix(c(rep(c(1, -1),10), (-1)^{1:20}/8), ncol=2),
  n = 20*60+60,
  initial=c(-1,1/16))



model.mix1 <-  list(name = "mix1",
  sgnl.type = "PWLC",
  cpt = (1:7) * 256,
  chg.size = c(1,-1,-1,1,1,-2,2)/64,
  n = 2048,
  initial=c(0,0))

model.mix2 <-  list(name = "mix2",
  sgnl.type = "PWL",
  cpt = (1:7) * 256,
  chg.size = matrix(c(c(0,-1,0,-1,1,1,0), c(1,-1,-1,1,1,-2,2)/64), ncol=2),
  n = 2048,
  initial=c(0,0))



model.mix3 <- list(name = "mix3",
  sgnl.type = "PWL",
  cpt = sort(c((1:7)*256, 512+12, 1280+9, 1792+6)),
  chg.size = matrix(c(c(-4, 6,-4,-1,1,-5,4,2,-7,7), c(0, 2,-3,2,-2,1,0,1,0,-2)/64), ncol=2),
  n = 2048,
  initial=c(2,0))



model.linsgmts <-  list(name = "linsgmts", # bump function length=10 / jump size=2
  sgnl.type = "PWL",
  cpt = sort(c((1:4)*2*256, (1:4)*2*256+5)),
  chg.size = matrix(c(rep(c(6,-6-4/64), 4), rep(c(1,-1)/64,4)), ncol=2),
  n = 256*9,
  initial=c(-1,0))

model.teeth <-  list(name = "teeth",
  sgnl.type = "PWC",
  cpt = (1:7) * 100,
  chg.size = 2*(-1)^{1:7},
  n = 100 * 8,
  initial = 1)

model.lin <-  list(name = "lin",
  sgnl.type = "PWLC",
  cpt = c(),
  chg.size = c(),
  n = 150* 10,
  initial=c(-1, 2/1500))

models <- list(model.wave1, model.wave2 , model.mix1,
               model.mix2, model.mix3, model.linsgmts,
               model.teeth, model.lin)
#par(mfrow=c(4,2),mar=rep(3,4))

#s=get_signal(model.wave7_random)
#plot(s)

#######################################
#   functions
#######################################

get_signal <- function(model){
  
  if(model$sgnl.type == "PWC"){
    
    sgnl <- rep(0, model$n)
    sgmts <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    sgnl[sgmts[1,1]:sgmts[1,2]] <- model$initial[1]
    
    for(j in 2:nrow(sgmts)){
      sgnl[sgmts[j,1]:sgmts[j,2]] <- sgnl[sgmts[j,1]-1] + model$chg.size[j-1]
    } 
    
  }else if(model$sgnl.type == "PWLC" & length(model$cpt)>0){
    
    sgnl <- rep(0, model$n)
    sgmts <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    slope <- model$initial[2]
    sgnl[sgmts[1,1]:sgmts[1,2]] <- model$initial[1] + sgmts[1,1]:sgmts[1,2] * model$initial[2]
    
    for(j in 2:nrow(sgmts)) {
      
      slope <- slope +  model$chg.size[j-1]
      for(k in sgmts[j,1]:sgmts[j,2]){
        sgnl[k] <- sgnl[k-1] + slope
      } 
    }
    
  }else if(model$sgnl.type == "PWLC" & length(model$cpt)==0){
    sgnl <- rep(0, model$n)
    sgmts <- matrix(c(1, model$n), ncol=2)
    slope <- model$initial[2]
    sgnl[sgmts[1,1]:sgmts[1,2]] <- model$initial[1] + sgmts[1,1]:sgmts[1,2] * model$initial[2]
    
  }else if(model$sgnl.type == "PWL"){
    
    sgnl <- rep(0, model$n)
    sgmts <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    slope <- model$initial[2]
    sgnl[sgmts[1,1]:sgmts[1,2]] <- model$initial[1] + sgmts[1,1]:sgmts[1,2] * slope
    
    for(j in 2:nrow(sgmts)) {
      slope <- slope +  model$chg.size[j-1,2]
      sgnl[sgmts[j,1]] <-  sgnl[sgmts[j-1,2]] + model$chg.size[j-1,1]
      
      if(sgmts[j,1]+1 < sgmts[j,2]){
        for(k in (sgmts[j,1]+1):sgmts[j,2]){
          sgnl[k] <- sgnl[k-1] + slope
        } 
      }
    }
  }
  
  return(sgnl)
  
}
finding_dH <- function(chp, modelnum, models){
  
  n <- models[[modelnum]]$n
  est.pnts <- sort(unique(c(0, chp, n)))
  true.pnts <- sort(unique(c(0, models[[modelnum]]$cpt, n)))
  
  d <- abs(matrix(rep(est.pnts, length(true.pnts)), nrow=length(est.pnts))
    -matrix(rep(true.pnts, length(est.pnts)), nrow=length(est.pnts), byrow=T))
  
  D1 <- max(apply(d, 2, min)) * 100 / n
  D2 <- max(apply(d, 1, min)) * 100 / n
  
  dH <- mean((abs(D1-D2) + D1 + D2)/2)
  
  return(dH)
  
}

assess <- function(object, modelnum, models){
  
  qdiff <- length(which(object$cpts>0)) - length(models[[modelnum]]$cpt)
  mse <- mean((get_signal(models[[modelnum]])-object$fit)^2)
  dh <- finding_dH(chp=object$cpts, modelnum=modelnum, models=models)
  
  return(list(qdiff=qdiff, mse=mse, dh=dh))
  
}
