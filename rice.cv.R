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
##pcr,plsr
model <- "pcr"

##output
dir.create("Rice",showWarnings=FALSE)
opfn <- paste("./Rice/RIL.meta.", model, ".cv.csv", sep="")

###fold id
fold <- as.numeric(read.csv(file="./DataIn/RIL.foldid.csv")[,1])
foldsegmnt <- lapply(1:10,function(i){
   ii <- which(fold==i)
   return(ii)
})
### 

######################################
### cross validation,  pls package ###
######################################
r2.array <- lapply(1:nphe, function(k){
   y <- as.numeric(phe[k,])
###pcr  
   if (model=="pcr"){
      RILfit <- pcr(y~z,validation="CV",segments = foldsegmnt)
   }
##plsr
   if (model=="plsr"){
      RILfit <- plsr(y~z,validation="CV",segments = foldsegmnt)
   }
##results
   ncomp <- as.numeric(RILfit$ncomp)
   r2.seq <- as.numeric(R2(RILfit)$val[1,,])
   r2.seq <- c(ncomp,r2.seq)
   return(r2.seq)
})
###
r2.array <- do.call(cbind,r2.array)
ncomp <- r2.array[1,1]
r2.array <- r2.array[-1,]
r2.array <- as.data.frame(cbind(0:ncomp,r2.array))
names(r2.array) <- c("ncomp",paste("phe", 1:nphe, sep=""))
write.csv(r2.array, file=opfn, row.names=FALSE)
###
