#setwd("D:/Dropbox/Project5_PCA")
rm(list=ls())
library(pls)
library(parallel)
library(glmnet)
source("./R/cal.xx.R")
source("./R/cv.foldid.R")
source("./R/Mixed.efficient.R")
source("./R/cv.BLUP.R")
source("./R/cv.LASSO.R")

                               ##################################################################
                               #### Genomic selection (GS) program of cross validation (CV)  ####
                               ######         For pcr, plsr, blup, lasso methods          #######
                               #######                2017-12-01                      ###########
                               ##################################################################
                               ############   ###      #####      ###      ######################
                                               ###    ### ###    ###
                                                ###  ###   ###  ###
                                                 ######     ######
                                                  ####       ####
##################
### Data input ###
##################
##phe<-read.csv(file="./DataIn/RIL.phenew.scale.csv")
##phe<-t(phe[,-1])
##phe<-read.csv(file="./DataIn/RIL.meta.scale.csv")
phe<-read.csv(file="./DataIn/RIL.mrna.scale.csv")
gen<-read.csv(file="./DataIn/RIL.gen.csv")
z<-as.matrix(t(gen))
z<-scale(z)
nindi<-ncol(phe)
nphe<-nrow(phe)
##pcr,blup
#models<-"blup" 
##mrna,meta
datatype<-"mrna"
#phename<-c("yield","kgw","grain","tiller")
phename<-c(datatype,1:nphe,sep="")
mc<-32
models<-"pcr"


###################################################################
### main program 1, principal component regression, Hat methods ###
################################################################### 
if(models=="pcr"){
   s1<-Sys.time()   
   foldid<-as.vector(read.csv(file="./DataIn/RIL.foldid.csv")[,1])
   foldsegmnt<-lapply(1:10,function(i){
      ii<-which(foldid==i)
      return(ii)
   })
###
   y<-as.numeric(phe[1,])
   RILfit<-pcr(y~z,validation="CV",segments=foldsegmnt)
   r2.seq<-as.numeric(R2(RILfit)$val[1,,])
   ii<-which.max(r2.seq)-1
   ncomp<-RILfit$ncomp
   scores<-RILfit$scores
   DD<-diag(t(scores)%*%scores)

   r2.array<-mclapply(1:nphe,function(k){
      y<-as.numeric(phe[k,])
      sst<-sum(y^2)
###
      foldid<-foldid
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
      return(predloo)
   },mc.cores=mc)
   r2.array<-do.call(cbind,r2.array)
   r2.array<-as.data.frame(cbind(1:ncomp,r2.array))
   names(r2.array)<-c("ncomps",phename)
   opfn<-paste("./DataOut/Rice.",datatype,"/",datatype,".pcr.hat.csv",sep="")
   write.csv(r2.array,file=opfn,row.names=FALSE)
   s2<-Sys.time()
   cat(s1,s2,s2-s1,"\n")
   cat("done","\n")
}


models<-"blup"
if(models=="blup"){
#########################################
### main program 2, blup, Hat methods ###
#########################################
###
s1<-Sys.time()
load(file="./DataIn/RIL.kk.eigen.RData")
delta<-kk.eigen[[1]] 
uu<-kk.eigen[[2]]
##
r2.array<-mclapply(1:nphe,function(k){
   y<-as.numeric(phe[k,])
   sst<-sum(y^2)
   d<-data.frame(y=y,x=rep(1,nindi))
   parms<-mixed.e(dataframe=d,kk.eigen=kk.eigen)
   yhat<-parms$random
   res<-y-as.numeric(parms$fixed)-yhat
   ress<-sum(res^2)
   good<-1-ress/sst
##predictability
   lambda<-as.numeric(parms$lambda)
   sg2<-as.numeric(parms$sg2)
   se2<-as.numeric(parms$se2)
   dd<-delta*sg2/(delta*sg2+se2)
   H<-sweep(uu,2,dd,"*")%*%t(uu)
   h<-diag(H)
   press<-sum((res/(1-h))^2)
   predloo<-1-press/sst
   rr<-c(good,predloo)
  return(rr)
},mc.cores=mc)              
r2.array<-as.data.frame(do.call(rbind,r2.array))
names(r2.array)<-c("goodness","predloo")    
opfn<-paste("./DataOut/Rice.",datatype,"/",datatype,".blup.hat.csv",sep="")
write.csv(r2.array,file=opfn,row.names=FALSE)
s2<-Sys.time()
cat(s1,s2,s2-s1,"\n")
cat("done","\n")
}
