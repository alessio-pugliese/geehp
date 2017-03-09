seasonal.FDC<-function(S,Q,cod){
##### S is a matrix or data.frame with the first column as a voctor expressing  
#     the months in a year (1:12) and the second column expresses the season to assign at    
##    the generic i-th month.  

#  S<-data.frame(month=1:12,season=c(1,1,1,2,2,3,3,3,3,3,4,1)) ## ex. several seasons in a year 
#  S<-data.frame(month=1:12,season=c(1,1,1,1,1,1,1,1,1,1,1,1)) ## ex. unique season in year
#  Q<-read.table("QMG_ABR_MAR_COMP.txt")
#  cod<-scan("bacini_preliminare.txt")  ## is a vector containing a subset of catchment codes
 
  num_seas<-length(unique(S[,2])) ### computes the no. of seasons 
  tot_sites_code<-unique(Q[,1]) ### vector with the overall no. of stations 
  N<-length(tot_sites_code) 
  n<-length(cod)
##### variable initialization
  Qs<-list()  
  years<-list()
  l<-c()
##### start computation
  for (i in 1:n){ ### the cycle accounts the subset of catchments
    Qs[[i]]<-as.list(1:num_seas)
    for (k in 1:num_seas){
      Qs[[i]][[k]]<-vector()
    }
    dummy<-Q[which(Q[,1]==cod[i]),]
    r<-dim(dummy)[1]  ### no of rows = no of years * 12 (months)
    years[[i]]<-unique(dummy[,2])
    l[i]<-length(years[[i]])
    for (j in 1:r){
        dis<-dummy[j,4:34]
        dis<-as.vector(t(dis))
        dis<-dis[-which(dis==-9999)]
        if (length(dis)==29)  ### discar record from leap years
        {dis<-dis[-29]} 
        read.month<-dummy[j,3]
        season<-S[read.month,2]
        if (length(Qs[[i]][[season]])==0)
          Qs[[i]][[season]]<-dis
        else
          Qs[[i]][[season]]<-c(Qs[[i]][[season]],dis)
    }    
    for (k in 1:num_seas){
      Qs[[i]][[k]]<-sort(Qs[[i]][[k]],decreasing=TRUE)
    }
  }
return(list(Qs=Qs,Y=years))  ### output disharges (Qs) and a vector with the years of record (Y)
}
  
