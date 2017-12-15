resample.FDC <- function(fdc,sam,norm){
  l<-length(fdc)
  d<-1:l/(l+1)                    ### Weibull plotting position
  #d<-(1-0.375):(l-0.375)/(l+0.25) ### Blom plotting position 
  
  FDCr<-c()
  if (norm==TRUE){
    FDCr=approx(qnorm(d),fdc,rule=2,method="linear",xout=qnorm(sam))
  }  else {
    FDCr=approx(d,fdc,rule=2,method="linear",xout=sam)
  }
return(FDCr)
}
