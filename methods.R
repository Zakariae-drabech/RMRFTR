
#######################################
#   Methods
#######################################
sink("file")


require(trendsegmentR)
require(IDetect)
require(not)
require(genlasso)
require(xtable)
require(dplyr)
require(cpop)

require(nsp)
require(dSTEM)
require(ChangePointInference)


library(trendsegmentR)
library(IDetect)
library(genlasso)
library(not)
library(cpop) # apt-get install libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev
library(nsp)
library(dSTEM)
library(ChangePointInference) # remotes::install_github("gaviosha/ChangePointInference")

mSTEM <- function(x,gamma=40,typ="II-linear" ){
  
  tic <- proc.time()
  object <-  dSTEM::dstem(x,typ,gamma=gamma,alpha=0.05)
  #"I" if the change points are piecewise linear and continuous;
  #"II-step" if the change points are piecewise constant and noncontinuous; 
  #"II-linear" if the change points are piecewise linear and noncontinuous;
  #"mixture" if both type I and type II change points are include in data
  
  #cpts=c(object$vall,object$peak)
  cpts= dSTEM::est.pair(vall=object$vall,peak=object$peak,gamma=gamma)$cp
 
  
  est <- rep(0, length(x))
  if (length(cpts)==0) {
    lmfit=lm(x[1:length((x))]~c(1:length(x)))
    est[1:length((x))]=lmfit$fitted.values
    #cat("nbr cpts:", length(cpts), "\n")
  }
  else{  
  for(i in 1:length(cpts)){
    if (i == 1) {
      lmfit <- lm(x[1:cpts[1]]~c(1:cpts[1]))
      est[c(1:cpts[1])] <- lmfit$fitted.values
    
    } 
    else  {
      lmfit <- lm(x[(cpts[i-1]+1):cpts[i]]~c((cpts[i-1]+1):cpts[i]))
      est[c((cpts[i-1]+1):cpts[i])] <- lmfit$fitted.values
    }
  }
    lmfit <- lm(x[(cpts[length(cpts)] + 1):length(x)]~c((cpts[length(cpts)] + 1):length(x)))
    est[c((cpts[length(cpts)] + 1):length(x))] <- lmfit$fitted.values
    }
  
  toc <- proc.time()
  cpts=as.numeric(cpts)
  
  list(fit = est, cpts=cpts, elapsed=(toc-tic)[3])
  
}



diff_Inf <- function(x ,alpha=0.1, degree = 1){
  #x=as.numeric(x)
  tic <- proc.time()
  object <- ChangePointInference::diffInf(x, degree = 1,alpha=alpha)
  cpts=ChangePointInference::cpt(object)
  fit =predict(object)
  fit=as.numeric(unlist(fit))

  toc <- proc.time()
  cpts=as.numeric(cpts)
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
  
}
ts <- function(x, thr=1.3, p=0.04, bal=0,cont=F){
 
  tic <- proc.time()
  object <- trendsegment(x,th.const=thr, p = p, bal = bal,continuous = cont, indep = TRUE)
  toc <- proc.time()
  cpts=as.numeric(object$cpt)
  list(fit = object$est, cpts=cpts, elapsed=(toc-tic)[3])

}



not_sic <- function(x,contrast="pcwsLinMean"){
  # pcwsLinMean, pcwsLinMean

  tic <- proc.time()
  object <- not(x, method = "not", contrast = contrast, parallel = FALSE) # give an error in noiseless input
  cpts <- features(object, penalty="sic")$cpt
  fit <- predict(object, cpt=cpts)
  toc <- proc.time()
  cpts=as.numeric(cpts)
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}



id <- function(x,contrast ="slope",  Lambda=1,ht=F){ # adjusted cpts after fitting
  # e.g. x <- c(0:9, 5, 10:0)+rnorm(22)/3  gives three consecutive change-points
  tic <- proc.time()
  #object <- ID(x, contrast ="slope", ht=F) # gives a weird cpts in noiseless input
  object <- ID(x, contrast =contrast, ht=ht, lambda=Lambda) # this does not make any difference compared to the upper one in MODEL 5-7 which contain spikes
  if(length(object$cpt)==1 && object$cpt==c(0)){
    cpts <- c()
  } else{
    cpts <- object$cpt
  }

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()
  cpts=as.numeric(cpts)
  
  list(fit = object$fit, cpts=cpts, elapsed=(toc-tic)[3])

}



tf <- function(x,ord_TF=1){ # adjusted cpts after fitting
  # e.g. x <- c(rep(0,9), 2, 10:0) gives three consecutive cpts
  tic <- proc.time();
  object <-  trendfilter(y=x, ord=ord_TF);

  tf.cv <- cv.trendfilter(object);
  cpts <- which(abs(diff(object$fit[,tf.cv$i.min], differences=2)) > sqrt(.Machine$double.eps))+1;

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()

  list(fit = object$fit[,tf.cv$i.min], cpts=cpts, elapsed=(toc-tic)[3])

}
  # Code that generates warning messages 

c_pop <- function(y){ # adjusted cpts after fitting
  y= as.numeric(y);
    # e.g. x <- c(0:9, 5, 10:0)+rnorm(22)/3 returns three consecutive cpts
  tic <- proc.time();
  x=1:length(y) - 1
  object <- cpop::cpop(y,x = x,
                       grid = x,
                       beta = 2 * log(length(y)),
                       sd = sqrt(mean(diff(diff(y))^2)/6),
                       minseglen = 0,
                       prune.approx = FALSE) ;#CPOP.run(
  cpts <- cpop::changepoints(object);
  cpts=cpts$location;
  fit <-cpop::estimate(object);
  fit=fit$y_hat;

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()
  cpts=as.numeric(cpts)
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}


bup <- function(x, max.err=300){

  tic <- proc.time()

  left_x <- seq(1, length(x)-1, by=2)
  right_x <- left_x + 1
  right_x[length(right_x)] <- length(x)
  number_of_segments <- length(left_x)

  segment <- cbind(left_x, right_x, Inf)

  for(i in 1:(number_of_segments-1)){
    lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
    segment[i, 3] <- sum((lmfit$residuals)^2) # sse
  }

  while(min(segment[,3]) < max.err){ # max.err

    i <- which.min(segment[,3])

    if(i==1){
      lmfit <- lm(x[c(segment[i, 1]:segment[i+2, 2])]~ c(segment[i, 1]:segment[i+2, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse

      segment[i, 2] <- segment[i+1, 2]
      segment <- segment[-(i+1),,drop=F]
    } else if(i>1 && i < dim(segment)[1]-1){
      lmfit <- lm(x[c(segment[i, 1]:segment[i+2, 2])]~ c(segment[i, 1]:segment[i+2, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse

      segment[i, 2] <- segment[i+1, 2]
      segment <- segment[-(i+1),,drop=F]

      i <- i-1
      lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse
    } else{
      segment[i, 2] <- segment[i+1, 2]
      segment[i, 3] <- Inf
      segment <- segment[-(i+1),,drop=F]

      i <- i-1
      lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse
    }

  }

  ### change-points
  cpt <- c(segment[-dim(segment)[1], 2])
  ### estimated curve
  est <- rep(NA, length(x))
  for(i in 1:dim(segment)[1]){
    lmfit <- lm(x[segment[i,1]:segment[i,2]]~c(segment[i,1]:segment[i,2]))
    est[c(segment[i,1]:segment[i,2])] <- lmfit$fitted.values
  }


  toc <- proc.time()
  cpt=as.numeric(cpt)
  
  list(fit=est, cpts=cpt, elapsed=(toc-tic)[3])

}


