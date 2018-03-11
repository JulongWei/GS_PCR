rm(list=ls())
library(pls)
source("cal.xx.R")
source("Mixed.efficient.R")
source("cv.BLUP.R")

##################
### Data input ###
##################
phe <- read.csv(file="./DataIn/sim.phe.csv")
gen <- read.csv(file="./DataIn/sim.gen.csv")
###pcr,blup
models <- "pcr"
###
y <- as.numeric(phe[,1])
y <- scale(y)
z <- as.matrix(t(gen))
z <- scale(z)
ny <- length(y)
sst <- sum(y^2)
##opfn
dir.create("Sim",showWarnings=FALSE)
opfn <- paste("./Sim/sim.",model,".hat.csv",sep="")
##


#######################
### PCR Hat methods ###
#######################
if(models=="pcr"){
   simfit <- pcr(y~z,validation="CV")
   ncomp <- simfit$ncomp
   scores <- simfit$scores
   DD <- diag(t(scores)%*%scores)
###
   predloo <- numeric(ncomp)
   H <- matrix(0,nindi,nindi)
   for(i in 1:ncomp){
      di <- as.numeric(1/DD[i])
      Ts <- scores[,i]  
      H <- H+Ts%*%t(Ts)*di
      h <- diag(H)
      yhat <- H%*%y
      res <- y-yhat    
      press <- sum((res/(1-h))^2)
      predloo[i] <- 1-press/sst 
   }
   r2.array <- as.data.frame(ncomps=0:ncomp,r2.pred=predloo)
   write.csv(r2.array,file=opfn,row.names=FALSE)
}


########################
### BLUP Hat methods ###
########################
if(models=="blup"){
   kk <- tcrossprod(z)
   kk.eigen <- eigen(kk)
   delta <- kk.eigen[[1]]
   uu <- kk.eigen[[2]]
###
   d <- data.frame(y=y,x=rep(1,ny))
   parms <- mixed.e(dataframe=d,kk.eigen=kk.eigen)
   yhat <- parms$random
   res <- y-as.numeric(parms$fixed)-yhat
   ress <- sum(res^2)
   sst1 <- sum((y-as.numeric(parms$fixed))^2)
   good <- 1-ress/sst1
###
   lambda <- as.numeric(parms$lambda)
   dd <- delta*lambda/(delta*lambda+1)
   H <- sweep(uu,2,dd,"*")%*%t(uu)
   h <- diag(H)
   press <- sum((res/(1-h))^2)
   pred <- 1-press/sst1
   predloo <- data.frame(goodness=good,preds=pred)
   write.csv(predloo,file=opfn,row.names=FALSE)
}




