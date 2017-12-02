##function-cv.foldid
cv.foldid<-function(indi,seed=NULL,cv){
   set.seed(seed=seed)
   
   tmp<-1:indi
   size<-indi%/%cv
   foldid<-rep(0,indi)

   for ( i in 1:cv){
       if (i == cv) size=length(tmp)
       p.test<-sample(x=tmp,size=size) 
       foldid[p.test]<-i
       tmp<-setdiff(tmp,p.test)
   }
   foldid<-as.vector(foldid)
   return(foldid)
}
