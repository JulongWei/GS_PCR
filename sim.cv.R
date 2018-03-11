rm(list=ls())
library(pls)
source("cal.xx.R")
source("Mixed.efficient.R")


##################
### Data input ###
##################
phe <- read.csv(file="./DataIn/sim.phe.csv")
gen <- read.csv(file="./DataIn/sim.gen.csv")
fold <- read.csv(file="./DataIn/sim.foldid.csv")  
##pcr,plsr
model <- "pcr"
###
y <- as.numeric(phe[,1])
y <- scale(y)
z <- as.matrix(t(gen))
z <- scale(z)
ny <- length(y)
sst <- sum(y^2)
nreps <- ncol(fold)
nreps<-1
###output file
dir.create("Sim",showWarnings=FALSE)
opfn <- paste("./Sim/sim.",model,".cv.csv",sep="")

#############################################################
### cross validation, pls package, including PCR and PLSR ###
#############################################################                   
r2.array <- lapply(1:nreps,function(i){
##foldid
   foldid <- fold[,i]
   foldsegmnt <- lapply(1:10,function(i){
      ii <- which(foldid==i)
      return(ii)
   })
##pcr
   if (model=="pcr"){
      simfit <- pcr(y~z,validation="CV",segments = foldsegmnt)
   }
##plsr
   if (model=="plsr"){
      simfit <- plsr(y~z,validation="CV",segments = foldsegmnt)
   }
##results
   r2.seq <- as.numeric(R2(simfit)$val[1,,])
   ncomp <- simfit$ncomp
   r2.seq <- c(ncomp,r2.seq)
   return(r2.seq)
})
###
r2.array <- do.call(cbind, r2.array)
ncomp <- r2.array[1,1]
r2.array <- r2.array[-1,]
r2.array <- as.data.frame(cbind2(0:ncomp,r2.array))
names(r2.array) <- c("ncomp",paste("r2.pred",1:nreps,sep=""))
write.csv(r2.array,file=opfn,row.names=FALSE)