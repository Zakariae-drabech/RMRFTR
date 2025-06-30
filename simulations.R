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

#source("G:/My Drive/Doctorat/Offline_CPD_PWL/trendsegment-master/methods.R")


#######################################
#   Models
#######################################
#sample(1:1980, 20)
ll=rep(1,20)
ll[3]=0
ll[9]=0
ll[14]=0
model.wave7_random <- list(name = "wave_random4",
                           sgnl.type = "PWLC",
                           cpt = c(100,180,230,300,390,490,580,660,770,840,930,1100,1230,1320,1400,1550,1680,1770,1850,1940),
                           chg.size = ll* (-1)^{1:20} / 25,
                           n = 2000,
                           initial = c(-1, 0))
model.wave6_random <- list(name = "wave_random3",
                           sgnl.type = "PWLC",
                           cpt = c(100,180,230,300,390,490,580,660,770,840,930,1100,1230,1320,1400,1550,1680,1770,1850,1940),
                           chg.size =  (-1)^{1:20} / sample(40:50, 20,replace = TRUE),
                           n = 2000,
                           initial = c(-1, 1/50))
model.wave5_random <- list(name = "wave_random2",
                           sgnl.type = "PWLC",
                           cpt = c(100,180,230,300,390,490,580,660,770,840,930,1100,1230,1320,1400,1550,1680,1770,1850,1940),
                           chg.size =  (-1)^{1:20} / 25,
                           n = 2000,
                           initial = c(-1, 1/50))
model.wave_randomEachTime <- list(name = "wave_random1",
                           sgnl.type = "PWLC",
                           cpt = sort(sample(1:1980, 20)),
                           chg.size =  (-1)^{1:20} / 10,
                           n = 2000,
                           initial = c(-1, 1/50))
model.wave1_500 <-  list(name = "wave1",
                     sgnl.type = "PWLC",
                     cpt = (1:3) * 150,
                     chg.size = (-1)^{1:3} / 25,
                     n = 500,
                     initial=c(-1,1/50))
model.wave1_1000 <-  list(name = "wave1",
                     sgnl.type = "PWLC",
                     cpt = (1:6) * 150,
                     chg.size = (-1)^{1:6} / 25,
                     n = 1000,
                     initial=c(-1,1/50))
model.wave1 <-  list(name = "wave1",
  sgnl.type = "PWLC",
  cpt = (1:9) * 150,
  chg.size = (-1)^{1:9} / 25,
  n = 150* 10,
  initial=c(-1,1/50))
model.wave1_2000 <-  list(name = "wave1",
                          sgnl.type = "PWLC",
                          cpt = (1:13) * 150,
                          chg.size = (-1)^{1:13} / 25,
                          n = 2000,
                          initial=c(-1,1/50))
model.wave1_2500 <-  list(name = "wave1",
                          sgnl.type = "PWLC",
                          cpt = (1:16) * 150,
                          chg.size = (-1)^{1:16} / 25,
                          n = 2500,
                          initial=c(-1,1/50))
model.wave1_3000 <-  list(name = "wave1",
                          sgnl.type = "PWLC",
                          cpt = (1:19) * 150,
                          chg.size = (-1)^{1:19} / 25,
                          n = 3000,
                          initial=c(-1,1/50))
model.wave2 <-  list(name = "wave2",
  sgnl.type = "PWL",
  cpt = (1:20) * 60,
  chg.size = matrix(c(rep(c(1, -1),10), (-1)^{1:20}/8), ncol=2),
  n = 20*60+60,
  initial=c(-1,1/16))
model.wave3 <-  list(name = "wave3", 
                     sgnl.type = "PWLC",
                     cpt = c(256, 512, 768, 1024, 1152, 1280, 1344),
                     chg.size = (-1)^(1:7) *(1:7)/ 64, 
                     n = 1408,
                     initial=c(1,1/256))
model.wave4 <-  list(name = "wave4", 
                     sgnl.type = "PWLC",
                     cpt = (1:9) * 150,
                     chg.size = (-1)^{1:16} / 32, 
                     n = 150* 10,
                     initial=c(-1/2,1/64))

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

model.mix2_ <-  list(name = "mix2_",
  sgnl.type = "PWL",
  cpt = sort(c((1:7)*256, 257, 1793)),
  chg.size = matrix(c(c(-7,7,2,-2,1,-1,1,-7,7), c(0,-1,1,1,-2,2.5,-1.5,0,-1.5)/64), ncol=2),
  n = 2048,
  initial=c(2,0))

model.mix3 <- list(name = "mix3",
  sgnl.type = "PWL",
  cpt = sort(c((1:7)*256, 512+12, 1280+9, 1792+6)),
  chg.size = matrix(c(c(-4, 6,-4,-1,1,-5,4,2,-7,7), c(0, 2,-3,2,-2,1,0,1,0,-2)/64), ncol=2),
  n = 2048,
  initial=c(2,0))

model.mix3_old <- list(name = "mix3",
                   sgnl.type = "PWL",
                   cpt = sort(c((1:7)*256, 512+30, 1280+30, 1793)),
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
modelsCPU<-list(model.wave1_500, model.wave1_1000 ,model.wave1 ,model.wave1_2000,
                model.wave1_2500,model.wave1_3000)
models <- list(model.wave1, model.wave2 , model.mix1,
               model.mix2, model.mix3, model.linsgmts,
               model.teeth, model.lin,model.wave3, model.wave4,
               model.wave5_random)
models2 <- list(model.wave3, model.wave4 )
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



#######################################
#   Simulations
#######################################

#for(i in 1:length(models)){
#  truex <- get_signal(models[[i]])
#  x <- truex + rnorm(length(truex), sd=1)
#  plot(x, type="l", ylab="", xlab="t", ylim=range(c(x, truex)), col=8, main=paste("( M",i,")" ,models[[i]]$name))
#  lines(truex, col=1)
#}

sims<- function(N, modelnum, sgma, thr, p, bal, noise,maxerr, gamma,ht=FALSE){
  
  truex <- get_signal(models[[modelnum]]) # sgnl only
  n <- length(truex)
  
  par(mfrow=c(1,1),mar=c(2,2,1.5,1.5))
  
  ### Assesement
  result <- list(qdiff=matrix(NA, nrow=N, ncol=8), mse=matrix(NA, nrow=N, ncol=8), dh=matrix(NA, nrow=N, ncol=8), time=matrix(NA, nrow=N, ncol=8))
  result <- lapply(result, function(x) {colnames(x) <- c("TS", "NOT", "ID", "TF", "CPOP", "BUP","diffInf","mSTEM"); x})
  
  for(K in 1:N){
    print(K)
    set.seed(K)
    
    if(noise=="norm1"){
      x <- truex + rnorm(n, sd=sgma)
    } else if(noise=="t5"){
      x <- truex + rt(n, df=5) #* sqrt(3/5)
    } else if(noise=="laplace"){
      x <- truex + extraDistr::rlaplace(n, mu = 0, sigma = 1)
    } else if(noise=="ar103"){
      x <- truex + arima.sim(n=n, list(ar=c(0.3)), rand.gen = function(n) rnorm(n, sd=sgma))
    }
    
    plot(x, type="l", col=8, xlab="t", ylab="",
         ylim=range(x),
         main=paste(models[[modelnum]]$name, " ( seed = ", K, " )"))
    lines(truex, col=7, lwd=2)
    
    
    #####################################################################
    ################################## TS ###############################
    #####################################################################
    obj <- ts(x, thr=thr, p=p, bal=bal) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,1] <- a.obj$qdiff
    result[[2]][K,1] <- a.obj$mse
    result[[3]][K,1] <- a.obj$dh
    result[[4]][K,1] <- obj$elapsed
    
    # lines(obj$fit, col=1, lwd=1, lty=1)
    #####################################################################
    ################################# NOT ###############################
    #####################################################################
    obj <- not_sic(x) ## cpts is turned "NA" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,2] <- a.obj$qdiff
    result[[2]][K,2] <- a.obj$mse
    result[[3]][K,2] <- a.obj$dh
    result[[4]][K,2] <- obj$elapsed
    
    # lines(obj$fit, col=2, lwd=1, lty=2)
    #####################################################################
    ################################# ID ################################
    #####################################################################
    
    obj <- id(x, ht=ht) ## cpts is turned "0" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,3] <- a.obj$qdiff
    result[[2]][K,3] <- a.obj$mse
    result[[3]][K,3] <- a.obj$dh
    result[[4]][K,3] <- obj$elapsed
    
    # lines(obj$fit, col=3, lwd=1, lty=3)
    #####################################################################
    ########################## Trend Filtering ##########################
    #####################################################################
    obj <- tf(x) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,4] <- a.obj$qdiff
    result[[2]][K,4] <- a.obj$mse
    result[[3]][K,4] <- a.obj$dh
    result[[4]][K,4] <- obj$elapsed
    
    # lines(obj$fit, col=4, lwd=1, lty=4)
    #####################################################################
    ################################ CPOP ###############################
    #####################################################################
    obj <- c_pop(x) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,5] <- a.obj$qdiff
    result[[2]][K,5] <- a.obj$mse
    result[[3]][K,5] <- a.obj$dh
    result[[4]][K,5] <- obj$elapsed
    
    # lines(obj$fit, col=5, lwd=1, lty=5)
    
    #####################################################################
    ################################# BUP ###############################
    #####################################################################
    obj <- bup(x, max.err = maxerr) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,6] <- a.obj$qdiff
    result[[2]][K,6] <- a.obj$mse
    result[[3]][K,6] <- a.obj$dh
    result[[4]][K,6] <- obj$elapsed
    
    # lines(obj$fit, col=6, lwd=1, lty=6)
    
    
    #####################################################################
    ################################# diff_Inf ###############################
    #####################################################################
    obj <- diff_Inf(x, alpha=0.9) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,7] <- a.obj$qdiff
    result[[2]][K,7] <- a.obj$mse
    result[[3]][K,7] <- a.obj$dh
    result[[4]][K,7] <- obj$elapsed
    
    # lines(obj$fit, col=7, lwd=1, lty=7)
    #####################################################################
    ################################# diff_Inf ###############################
    #####################################################################
    obj <- mSTEM(x, gamma = gamma) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,8] <- a.obj$qdiff
    result[[2]][K,8] <- a.obj$mse
    result[[3]][K,8] <- a.obj$dh
    result[[4]][K,8] <- obj$elapsed
    
    # lines(obj$fit, col=8, lwd=1, lty=8)
    legend("topleft", c("obs","true","TS","NOT","ID","TF","CPOP", "BUP","diffInf","mSTEM"), col=c(8,8,1:8), lty=c(1,1,1:8), lwd=c(1,2,rep(1,8)), cex=0.7, bty="n", ncol=2)
  }
  
  return(result)
  
}

sims_one_method<- function(N, modelnum, sgma, num_CPs=2, thr=1.3, p=0.04, bal=0, noise="norm1", maxerr=300, gamma=40,ht=FALSE){
  truex <- get_signal(models[[modelnum]]) # sgnl only
  n <- length(truex)
  
  par(mfrow=c(1,1),mar=c(2,2,1.5,1.5))
  
  ### Assesement
  result <- list(qdiff=matrix(NA, nrow=N, ncol=8), mse=matrix(NA, nrow=N, ncol=8), dh=matrix(NA, nrow=N, ncol=8), time=matrix(NA, nrow=N, ncol=8))
  result <- lapply(result, function(x) {colnames(x) <- c("TS", "NOT", "ID", "TF", "CPOP", "BUP","diffInf","mSTEM"); x})
  
  for(K in 1:N){
    set.seed(K)
    #print(K)
    if(noise=="norm1"){
      x <- truex + rnorm(n, sd=sgma)
    } else if(noise=="t5"){
      x <- truex + rt(n, df=5) #* sqrt(3/5)
    } else if(noise=="laplace"){
      x <- truex + extraDistr::rlaplace(n, mu = 0, sigma = 1)
    } else if(noise=="ar103"){
      x <- truex + arima.sim(n=n, list(ar=c(0.3)), rand.gen = function(n) rnorm(n, sd=sgma))
    }
    obj <- ts(x, thr=thr, p=p, bal=bal,cont=F)
    #obj <- not_sic(x)
    #obj <- id(x, ht=FALSE)
    #obj <- tf(x) ## cpts is turned "integer(0)" if nothing is detected
    #obj <- c_pop(x)
    #obj <- bup(x, max.err = maxerr)
    #obj <- diff_Inf(x, alpha=0.9)
    #obj <-  Segmented(x,num_CPs=num_CPs)
    #obj <- mSTEM(x, gamma = gamma) 
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,1] <- a.obj$qdiff
    result[[2]][K,1] <- a.obj$mse
    result[[3]][K,1] <- a.obj$dh
    result[[4]][K,1] <- obj$elapsed
    
  }
  return(result)
}

Res_table<- function(m1n1,namefile ){
  qdiff=m1n1$qdiff
  MSE=m1n1$mse
  DH=m1n1$dh
  TIME=m1n1$time
  result <- matrix(0, nrow = 8, ncol = 10)
  
  for (K in 1:8) {
    for (n in 1:100) {
      result[K, 1] <- result[K, 1] + ifelse(qdiff[n,K] <= -3, 1, 0)
      result[K, 2] <- result[K, 2] + ifelse(qdiff[n,K] == -2, 1, 0)
      result[K, 3] <- result[K, 3] + ifelse(qdiff[n,K] == -1, 1, 0)
      result[K, 4] <- result[K, 4] + ifelse(qdiff[n,K] == 0, 1, 0)
      result[K, 5] <- result[K, 5] + ifelse(qdiff[n,K] == 1, 1, 0)
      result[K, 6] <- result[K, 6] + ifelse(qdiff[n,K] == 2, 1, 0)
      result[K, 7] <- result[K, 7] + ifelse(qdiff[n,K] >= 3, 1, 0)
      
      
    }
    result[K, 8] <- round(mean(MSE[,K]),3)
    result[K, 9] <-  round(mean(DH[,K]), 3)
    result[K, 10] <- round(mean(TIME[,K]), 2)
  }
  
  
  result=list(res=result)
  result <- lapply(result, function(x) {rownames(x) <- c("TS", "NOT", "ID", "TF", "CPOP", "BUP","diffInf","mSTEM"); x})
  result <- lapply(result, function(x) {colnames(x) <- c("<=-3","-2", "-1", "0","1","2",">=3","MSE" , "dh(x10^2)",  "time"); x})
  #write.table(result$res,paste("G:/My Drive/Doctorat/Offline_CPD_PWL/MiniorRevisionArticle2/MajorRevisionArticle2/num_res/",namefile,sep=""), 
  #                        append = FALSE, sep = " ", dec = ".",
  #                        row.names = TRUE, col.names = TRUE)
  return(result)
  
}

#m1n1 <- sims(N=1, modelnum=1, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1", maxerr=300, gamma=40)
#m2n1 <- sims(N=100, modelnum=2, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m3n1 <- sims(N=100, modelnum=3, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m4n1 <- sims(N=100, modelnum=4, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m5n1 <- sims(N=100, modelnum=5, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m6n1 <- sims(N=100, modelnum=6, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m7n1 <- sims(N=100, modelnum=7, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
#m8n1 <- sims(N=100, modelnum=8, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
 
 
#m2_TS <- sims_one_method(N=100, modelnum=4, sgma=1 )
#tablem2_TS=Res_table(m1n1 =m2_TS," ")
#print(tablem2_TS)
