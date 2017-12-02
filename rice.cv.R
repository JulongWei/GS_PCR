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
phe<-read.csv(file="./DataIn/RIL.meta.scale.csv")
gen<-read.csv(file="./DataIn/RIL.gen.csv")
z<-as.matrix(t(gen))
z<-scale(z)
nindi<-ncol(phe)
nphe<-nrow(phe)

##meta,mrna
datatype<-"meta"
##pcr,plsr,blup,lasso
models<-"lasso"
mc<-64


#####################################################
### replicates cross validation,  pls package #######
#####################################################
if(models=="pcr"|models=="plsr"){
   s1<-Sys.time()
###                   
   folds<-read.csv(file="./DataIn/RIL.foldid.csv")
#pls.options(parallel=mc) 
   r2.array<-mclapply(1:nphe,function(k){
      y<-as.numeric(phe[k,])
      foldid<-folds[,1]
      foldsegmnt<-lapply(1:10,function(i){
         ii<-which(foldid==i)
         return(ii)
      }) 
  ##pcr
      if (models=="pcr"){
         RILfit<-pcr(y~z,validation="CV",segments = foldsegmnt)
      }
  ##plsr
      if (models=="plsr"){
         RILfit<-plsr(y~z,validation="CV",segments = foldsegmnt)
      }
  ##results
      ncomp<-as.numeric(RILfit$ncomp)
      r2.seq<-as.numeric(R2(RILfit)$val[1,,])
      r2.seq<-c(ncomp,r2.seq)
      #r2.max<-max(r2.seq)
      return(r2.seq)
   },mc.cores=mc)
###
   cat("\n","done","\n")
   r2.array<-do.call(cbind,r2.array)
   ncomp<-r2.array[1,1]
   r2.array<-r2.array[-1,]
   r2.array<-as.data.frame(cbind(0:ncomp,r2.array))
   names(r2.array)<-c("ncomp",paste(datatype,1:nphe,sep=""))
   opfn<-paste("./DataOut/Rice.",datatype,"/",datatype,".",models,".cv.csv",sep="")
   write.csv(r2.array,file=opfn,row.names=FALSE)
###
   s2<-Sys.time()
   cat(s1,s2,s2-s1,"\n")
}


####################################################
### replicates cross validation, BLUP and LASSO ####
####################################################
#models<-"blup"
y<-as.numeric(phe[1,])
RILfit<-pcr(y~z,validation="CV")
scores<-RILfit$scores
z<-as.matrix(scores)

if(models=="blup"|models=="lasso"){
  s1<-Sys.time()
  load(file="./DataIn/RIL.kk.eigen.RData")
  folds<-read.csv(file="./DataIn/RIL.foldid.csv")
  nfold<-ncol(folds)
  foldid<-folds[,1]
  
r2.array<-mclapply(1:nphe,function(k){
##
  y<-as.numeric(phe[k,])   
  d<-data.frame(y=y,x=rep(1,nindi)) 
  if (models=="blup"){
    parms<-cv.BLUP(dataframe=d,kk=kk,foldid=foldid)
  }
  ##lasso
  if(models=="lasso"){
    parms<-cv.LASSO(dataframe=d,gen=z,foldid=foldid)
  }
  r2<-as.numeric(parms$phe.r2)
  return(r2)
},mc.cores=mc)
###
  r2.array<-data.frame(r2=as.numeric(r2.array))
  opfn<-paste("./DataOut/Rice.",datatype,"/",datatype,".",models,".score.cv.csv",sep="")
  write.csv(r2.array,file=opfn,row.names=FALSE)
###     
s2<-Sys.time()
cat(s1,s2,s2-s1,"\n")
} 

