#setwd("D:/Researching/Project12_PCA")
rm(list=ls())
library(pls)
library(parallel)
library(glmnet)
source("./R/cal.xx.R")
source("./R/cv.foldid.R")
source("./R/Mixed.efficient.R")
source("./R/cv.BLUP.R")
source("./R/cv.LASSO.R")


##################
### Data input ###
##################
phe<-read.csv(file="./DataIn/sim-phe-2000.csv")
gen<-read.csv(file="./DataIn/sim-gen-2000.csv")
###
y<-as.matrix(phe$y)
y<-scale(y)
z<-as.matrix(t(gen[,-c(1:4)]))
z<-scale(z)
ny<-nrow(y)
sst<-sum(y^2)
###arguments required from program
###pcr,blup
models<-"pcr"
mc<-32

if(models=="pcr"){
   s1<-Sys.time()
###################################################################
### main program 1, principal component regression, Hat methods ###
###################################################################
###   
###
   foldid<-as.vector(read.csv(file="./DataIn/sim.foldid.csv")[,1])
   foldsegmnt<-lapply(1:10,function(i){
      ii<-which(foldid==i)
      return(ii)
   })
###
   pls.options(parallel=mc)
   simfit<-pcr(y~z,validation="CV",segments=foldsegmnt)
   r2.seq<-as.numeric(R2(simfit)$val[1,,])
   ii<-which.max(r2.seq)-1
   ncomp<-simfit$ncomp
   scores<-simfit$scores
   DD<-diag(t(scores)%*%scores)
###
   goodloo<-numeric(ncomp)
   predloo<-numeric(ncomp)
   H<-matrix(0,nindi,nindi)
   for(i in 1:ncomp){
      di<-as.numeric(1/DD[i])
      Ts<-scores[,i]  
      H<-H+Ts%*%t(Ts)*di
      h<-diag(H)
      yhat<-H%*%y
      res<-y-yhat
      ress<-sum(res^2)      
      press<-sum((res/(1-h))^2)
      goodloo[i]<-1-ress/sst
      predloo[i]<-1-press/sst 
   }
   r2.array<-as.data.frame(ncomps=1:ncomp,r2.pred=predloo,r2=goodloo)
   write.csv(r2.array,file="./DataOut/Sim/sim.pcr.hat.csv",row.names=FALSE)
   s2<-Sys.time()
   cat(s1,s2,s2-s1,"\n")
}



if(models=="blup"){
#########################################
### main program 2, blup, Hat methods ###
#########################################
###
   s1<-Sys.time()
   load(file="./DataIn/sim.kk.eigen.RData")
   delta<-kk.eigen[[1]] 
   uu<-kk.eigen[[2]]
###
   d<-data.frame(y=y,x=rep(1,ny))
   parms<-mixed.e(dataframe=d,kk.eigen=kk.eigen)
   yhat<-parms$random
   res<-y-as.numeric(parms$fixed)-yhat
   ress<-sum(res^2)
   sst1<-sum((y-as.numeric(parms$fixed))^2)
   good<-1-ress/sst1
###predictability
   lambda<-as.numeric(parms$lambda)
   dd<-delta*lambda/(delta*lambda+1)
   H<-sweep(uu,2,dd,"*")%*%t(uu)
   h<-diag(H)
   press<-sum((res/(1-h))^2)
   pred<-1-press/sst1
   predloo<-data.frame(goodness=good,preds=pred)
   write.csv(predloo,file="./DataOut/Sim/sim.blup.hat.csv",row.names=FALSE)
s2<-Sys.time()
cat(s1,s2,s2-s1,"\n")
}




