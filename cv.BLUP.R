#####prediction Exact kinship matrix
#########
cv.BLUP<-function(dataframe,gen,kk=NULL,foldid=NULL,seed=NULL){
#
   d<-dataframe
#   gen<-as.matrix(gen)
   indi<-nrow(d)
   cv<-max(foldid)
#
###
   p.total<-1:indi
   cv.var<-NULL
   tmp<-NULL
   rr<-NULL
   for ( i in 1:cv){

#data divided into two parts, one for reference and the other for test.
###
#   
       p1<-which(foldid==i)
       p0<-p.total[-p1]
       n0<-length(p0)
       n1<-length(p1)
#       
       d0<-d[p0,]
       x0<-as.matrix(d0[,-1])
       d1<-d[p1,]
       x1<-as.matrix(d1[,-1])
#       
#       gen0<-gen[p0,]
#       gen1<-gen[p1,]
#
       kk0<-kk[p0,p0]
       kk1<-kk[p1,p1]
       cov.TR<-kk[p1,p0]               #T short for test and R short for reference 

#parameter estimated in the reference
###
       kk0.eigen<-eigen(kk0,symmetric=TRUE)
       parms<-mixed.e(dataframe=d0,kk.eigen=kk0.eigen)     
   
       beta0<-as.matrix(parms$fixed)
       sg2<-parms$sg2
       se2<-parms$se2
       lambda<-parms$lambda 
       var0<-c(sg2,se2,lambda)
       cv.var<-rbind(cv.var,var0)                   
#           
       v<-kk0*sg2+diag(n0)*se2

######################################################################################################
### There are two type of R-squared, #################################################################
### (1) the correlation between phenotype and predicted prenotype ####################################
### (2) the correlation between phenotype-fixed effects and breeding value, just polygenic effects ###
######################################################################################################      
       predict1<-sg2*cov.TR%*%solve(v,tol=1e-50)%*%(d0[,1]-x0%*%beta0)
#       predic1<-x1[,-1,drop=FALSE]%*%beta0[-1]       
       mu<-as.matrix(x1)%*%as.matrix(beta0)
       y1<-d1[,1]
       yhat<-y1
       predict1<-predict1+mu
###       
###
       tmp0<-cbind(p1,mu,yhat,predict1)
       tmp<-rbind(tmp,tmp0)
       r0<-cor(y1,predict1)
       rr<-c(rr,r0)    
#                       
   }
   
#data output
###
   phe.r2<-cor(tmp[,3],tmp[,4])^2
#
   rownames(tmp)<-NULL
   tmp<-as.data.frame(tmp)
   names(tmp)<-c("id","fixed","yhat","predict")
#
   rownames(cv.var)<-NULL
   cv.var<-as.data.frame(cv.var)
   names(cv.var)<-c("sg2","se2","lambda")
   RR<-list(phe.r2=phe.r2,cv.var=cv.var,predic=tmp,foldid=foldid)
   return(RR)
}
