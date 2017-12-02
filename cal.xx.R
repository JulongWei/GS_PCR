##function-cal.xx
cal.xx<-function(x,y,H=NA){
##
  x<-as.matrix(x)
  y<-as.matrix(y)
  nr<-ncol(x)
  nc<-ncol(y)
  n0<-nrow(x)
##  
  if (is.na(H)) H<-diag(n0)
  xx<-t(x)%*%H%*%y
  xx<-as.matrix(xx)
  return(xx)
}