rm(list=ls())
library(pls)
source("cal.xx.R")
source("Mixed.efficient.R")

##################
### Data input ###
##################
phe <- read.csv(file="./DataIn/RIL.meta.scale.csv")
gen <- read.csv(file="./DataIn/RIL.gen.csv")
z <- as.matrix(t(gen))
z <- scale(z)
nindi <- ncol(phe)
nphe <- nrow(phe)
##pcr,blup
model <- "pcr"

##output filename
dir.create("Rice",showWarnings=FALSE)
opfn <- paste("./Rice/RIL.meta.", model, ".hat.csv", sep="")


#######################
### pcr Hat methods ###
#######################
if(model=="pcr"){
###
   y <- as.numeric(phe[1,])
   RILfit <- pcr(y~z,validation="CV")
   ncomp <- RILfit$ncomp
   scores <- RILfit$scores
   DD <- diag(t(scores)%*%scores)

   r2.array <- lapply(1:nphe,function(k){
      y <- as.numeric(phe[k,])
      sst <- sum(y^2)
###
      predloo <- numeric(ncomp)
      H <- matrix(0, nindi, nindi)
      for(i in 1:ncomp){
         di <- as.numeric(1/DD[i])
         Ts <- scores[,i]  
         H <- H + Ts%*%t(Ts)*di
         h <- diag(H)
         yhat <- H%*%y
         res <- y-yhat
         press <- sum((res/(1-h))^2)
         predloo[i] <- 1-press/sst
      }
      return(predloo)
   })
   r2.array <- do.call(cbind,r2.array)
   r2.array <- as.data.frame(cbind(1:ncomp, r2.array))
   names(r2.array) <- c("ncomps", paste("phe",1:nphe,sep=""))
   write.csv(r2.array, file=opfn, row.names=FALSE)
}

#######################
###blup Hat methods ###
#######################
if(model=="blup"){
   kk <- tcrossprod(z)
   kk.eigen <- eigen(kk)
   delta <- kk.eigen[[1]] 
   uu <- kk.eigen[[2]]
##
   r2 <- lapply(1:nphe,function(k){
      y <- as.numeric(phe[k,])
      sst <- sum(y^2)
      d <- data.frame(y=y,x=rep(1,nindi))
      parms <- mixed.e(dataframe=d,kk.eigen=kk.eigen)
      yhat <- parms$random
      res <- y-as.numeric(parms$fixed)-yhat
###
      lambda <- as.numeric(parms$lambda)
      sg2 <- as.numeric(parms$sg2)
      se2 <- as.numeric(parms$se2)
      dd <- delta*sg2/(delta*sg2+se2)
      H <- sweep(uu,2,dd,"*")%*%t(uu)
      h <- diag(H)
      press <- sum((res/(1-h))^2)
      predloo <- 1-press/sst
      return(predloo)
   })              
   r2 <- do.call(c, r2)
   r2 <- data.frame(r2=r2)
   write.csv(r2, file=opfn, row.names=FALSE)
}
