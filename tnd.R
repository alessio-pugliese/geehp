### Caricare il pacchetto 'pracma'
### IMPORTANTE!! la FDC deve essere già standardizzata
fdc.tnd<-function(fdc,norm,maxd,log=TRUE){
  l<-length(fdc)
  d<-1:l/(l+1) ##### weibull
  #d<-(1-0.375):(l-0.375)/(l+0.25) ### Blom plotting position
  
  d1<-d[which(fdc<=1)[1]]  
  
  if (norm==TRUE){
  f<-approxfun(qnorm(d),fdc,rule=2)  
  tnd<-(qnorm(maxd)-qnorm(d1))-quad(f,xa=qnorm(d1),xb=qnorm(maxd))
  }
  
  else { 
  f<-approxfun(d,fdc,rule=2)
  tnd<-(maxd-d1)-quad(f,xa=d1,xb=maxd)
  }

return(tnd)
}