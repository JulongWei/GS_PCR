###function-mixed.e
mixed.e<-function(dataframe,gen=NULL,kk.eigen,stability=TRUE){
#
   d<-dataframe
   y<-as.matrix(d[,1])
   x<-as.matrix(d[,-1])
   n<-nrow(d)
   s<-ncol(x)
#  
   qq<-kk.eigen
   delta<-qq[[1]]
   uu<-qq[[2]]
## 
   yu<-t(uu)%*%y
   xu<-t(uu)%*%x 
##loglike() function   
   loglike<-function(theta){
      lambda<-exp(theta)
      d0<-sum(log(abs(lambda*delta+1)))      
      h<-1/(lambda*delta+1)
      yy<-cal.xx(x=yu,y=yu,H=diag(h))
      xy<-cal.xx(x=xu,y=yu,H=diag(h))
      xx<-cal.xx(x=xu,y=xu,H=diag(h))
      
      V<-abs(yy-t(xy)%*%solve(xx,tol=1e-100)%*%xy)
      d1<-unlist(determinant(xx))[1]
      loglike<--0.5*(d0+d1+(n-s)*log(V))
      return(-loglike)
   }


##solve parms 
if(stability){

   intervals<-c(-100,-15,-10,-5,0,20,40,80)
   intials<-c(-16,-12,-7,-2,12,30,60)
   ni<-length(intervals)-1
   thetas<-NULL
   for ( i in 1:ni){
       theta<-intials[i]
       parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=intervals[i],upper=intervals[i+1])
       theta<-parms$par
       parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=-100,upper=100)
       fi<-parms$par
       thetas<-c(thetas,fi)
   }
   fns<-sapply(1:ni,function(i){
       theta<-thetas[i]
       fn<-loglike(theta)
       return(fn)
       } )
   ii<-which.min(fns)
   lambda<-exp(thetas[ii])
   fn1<-fns[ii]
}else{
   theta<-0
   parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-100,upper=100)
   lambda<-exp(parm$par)
   conv<-parm$convergence
   fn1<-parm$value
}

##LRT 
   fn0<-loglike(-Inf)
   lrt<-2*(fn0-fn1)
   if (lrt<=1e-03){
      lrt.p<-0.5+pchisq(lrt,df=1,lower.tail=FALSE)*0.5
   }else{
      lrt.p<-1-(0.5+pchisq(lrt,df=1)*0.5)   
   }
##
   h<-1/(lambda*delta+1)
   xx<-cal.xx(x=xu,y=xu,H=diag(h))
   xy<-cal.xx(x=xu,y=yu,H=diag(h))   
   fixed<-solve(xx,xy,tol=1e-50)       
   yt<-yu-xu%*%fixed
   yPy<-sum(yt*h*yt)
   se2<-yPy/(n-s)
   sg2<-lambda*se2      
   vfixed<-solve(xx,tol=1e-50)*sg2
   #random<-uu%*%((delta*lambda)*h*yt)
   dd<-(delta*sg2)/(delta*sg2+se2)
   random<-uu%*%(dd*yt)
   #
   rr<-list(fixed=fixed,vfixed=vfixed,random=random,sg2=sg2,se2=se2,lambda=lambda,lrt=lrt,lrt.p=lrt.p) 
   class(rr)<-"mixed" 
   if(!is.null(gen)){ 
     cc<-qq[[3]]
     markers<-t(z)%*%uu%*%((1/cc)*lambda*h*yt)
     rr$markers<-markers
   }
   return(rr)
}