cv.LASSO<-function(dataframe,gen,foldid=NULL,cv,lambda,seed=NULL){
#
   d<-dataframe
   gen<-as.matrix(gen)
   indi<-nrow(d)
   if ( is.null(foldid)) foldid<-cv.foldid(indi=indi,cv=cv,seed=seed)
   cv<-max(foldid)
#
   y<-as.matrix(d[,1])
   x<-as.matrix(d[,-1])
   s<-ncol(x)
   w<-cbind(x,gen)[,-1]
   nw<-ncol(w)
#
#   fit<-cv.glmnet(x=w,y=y,foldid=foldid)
#   lambda<-fit$lambda.min
#
   p.total<-1:indi
   tmp<-NULL
   rr<-NULL
   for (i in 1:cv){
#       
       p1<-which(foldid==cv)
       p0<-p.total[-p1]
       np0<-length(p0)
       np1<-length(p1)
#
       y0<-d[p0,1]
       x0<-x[p0,]
       w0<-w[p0,]
#
       y1<-d[p1,1]
       x1<-x[p1,]
       w1<-w[p1,]
#
       ffit<-cv.glmnet(x=w0,y=y0)
       coeff<-coef(ffit,s="lambda.min")
       beta0<-as.matrix(coeff[1:s])
       gamma<-as.matrix(coeff[-c(1:s)])

######################################################################################################
### There are two type of R-squared, #################################################################
### (1) the correlation between phenotype and predicted prenotype ####################################
### (2) the correlation between phenotype-fixed effects and breeding value, just polygenic effects ###
######################################################################################################         
       mu<-as.vector(as.matrix(x1)%*%beta0)
       yhat<-y1
#
       z<-w1[,s:nw]
       predict1<-z%*%gamma+mu
       
###       
###
       tmp0<-cbind(p1,mu,yhat,predict1)
       tmp<-rbind(tmp,tmp0) 
       r0<-cor(yhat,predict1)
       rr<-c(rr,r0)             
   }

# output
#
   phe.r2<-cor(tmp[,3],tmp[,4])^2
#
   tmp<-as.data.frame(tmp)
   names(tmp)<-c("id","mu","y","predict")
#
   RR<-list(rr=rr,phe.r2=phe.r2,predic=tmp,foldid=foldid)
   return(RR)
}