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
##pcr,plsr,blup,lasso
model<-"pcr"
mc<-64

#####################################################
### replicates cross validation,  pls package #######
#####################################################
if(model=="pcr"|model=="plsr"){
   s1<-Sys.time()
   ###                   
   foldids<-read.csv(file="./DataIn/sim.foldid.csv")
   nreps<-ncol(foldids)
   max.ncomp<-NULL
   r2.array<-NULL
   r2.array<-mclapply(1:nreps,function(i){
   ##foldid
      foldid<-foldids[,i] 
      pls.options(parallel=mc)
      foldsegmnt<-lapply(1:10,function(i){
         ii<-which(foldid==i)
         return(ii)
      })
   ##pcr
      if (model=="pcr"){
         simfit<-pcr(y~z,validation="CV",segments = foldsegmnt)
      }
   ##plsr
      if (model=="plsr"){
         simfit<-plsr(y~z,validation="CV",segments = foldsegmnt)
      }
   ##results
      r2.seq<-as.numeric(R2(simfit)$val[1,,])
      ncomp<-simfit$ncomp
      r2.seq<-c(ncomp,r2.seq)
      return(r2.seq)
   },mc.cores=mc)
   ###
   s2<-Sys.time()
   cat(s1,s2,s2-s1,"\n")
   ###
   r2.array<-do.call(cbind,r2.array)
   ncomp<-r2.array[1,1]
   r2.array<-r2.array[-1,]
   r2.array<-as.data.frame(cbind(1:ncomp,r2.array))
   names(r2.array)<-c("ncomp",paste("r2.pred",1:nreps,sep=""))
   opfn<-paste("./DataOut/Sim/sim.",model,".cv.csv",sep="")
   write.csv(r2.array,file=opfn,row.names=FALSE)
}


##########################################
### replicates cross validation, BLUP ####
##########################################
if(model=="blup"|model=="lasso"){
   s1<-Sys.time()
   d<-data.frame(y=y,x=rep(1,ny))
   load(file="./DataIn/sim.kk.eigen.RData")
   ###
   foldids<-read.csv(file="./DataIn/sim.foldid.csv")
   nfold<-ncol(foldids)
   ###
   r2.seq<-mclapply(1:nfold,function(i){
      foldid<-as.vector(foldids[,i])
   ##blup
      if (model=="blup"){
         parms<-cv.BLUP(dataframe=d,kk=kk,foldid=foldid)
      }
   ##lasso
   if(model=="lasso"){
      parms<-cv.LASSO(dataframe=d,gen=z,foldid=foldid)
   }
   r2<-as.numeric(parms$phe.r2)
   cat("folds:",i,"\n")
   return(r2)
   },mc.cores=mc)
   s2<-Sys.time()
   cat(s1,s2,s2-s1,"\n")
   r2.seq<-data.frame(x=as.numeric(r2.seq))
   opfn<-paste("./DataOut/Sim/sim.",model,".cv.r2.csv",sep="")
   write.csv(r2.seq,file=opfn,row.names=FALSE)
} 
