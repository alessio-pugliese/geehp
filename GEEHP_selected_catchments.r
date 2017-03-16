#source('GEEHP_FDCs_LOOCV.R')
source("extract_series.r")
library(zoo)
library(parallel)

##### Load streamflow data and catchment boundaries
set.seed(1)
predictionLocations <- readOGR("shp","tyrol_ehype")
load("Q.tyrol.RData") 
selected.catchments<-read.csv("Code_D_Area_selected_wo_removed.csv")


#######  Varaibles initialization 
#######
ID.ehype<-selected.catchments$EHYPE
l<-length(ID.ehype)
pos.ID.ehype<-c()
pos.ID.obs<-c()
est.FDC<-matrix(NA,nrow=np,ncol=l)
err.FDC<-matrix(NA,nrow=np,ncol=l)
ehype.FDC<-matrix(NA,nrow=np,ncol=l)
ehype.MAF<-c()
A.ehype<-c()

######   Top-kriging procedure for FDCs predictions at selected catchements
######
for (i in 1:l){
  print(paste("loop cycle no=",i))
  pos.ID.ehype[i]<-which(predictionLocations$macroid1==ID.ehype[i])
  pos.ID.obs[i]<-which(selected.catchments$OBS[i]==cod)
  
  ehype.catchments<-predictionLocations
  obs.to.ehype<-observations[-pos.ID.obs[i],]
  obs.to.ehype$obs<-TND[-pos.ID.obs[i]]
  tyrolRtopObj.FDC<-createRtopObject(obs.to.ehype,ehype.catchments,
                                   formulaString = obs~1,
                                   params = list(gDist = TRUE,debug.level=0,
                                                 rresol = 500,
                                                 nmax=vic,
                                                 #nclus=nr_cores,
                                                 wlim=1,
                                                 partialOverlap=TRUE
                                                 ))
  
                                              

  tyrolRtopObj.FDC<- rtopVariogram(tyrolRtopObj.FDC)
  tyrolRtopObj.FDC<- rtopFitVariogram(tyrolRtopObj.FDC)
  #tyrolRtopObj.FDC<- checkVario(tyrolRtopObj.FDC, cloud = TRUE, identify = TRUE,
  #                              acor = 0.01,log="")
  tyrolRtopObj.FDC<- rtopKrige(tyrolRtopObj.FDC,wret=TRUE)

  weights.obs.to.ehype<-tyrolRtopObj.FDC$weight[pos.ID.ehype[i],]
  
  #est.dimless.FDC<-t(weights.obs.to.ehype)%*%t(y[,-pos.ID.obs[i]])
  est.dimless.FDC<-as.vector(y[,-pos.ID.obs[i]]%*%(weights.obs.to.ehype))

  ##### MAF prediction at ehype locations
  A.ehype[i]<-as.numeric(gArea(ehype.catchments[pos.ID.ehype[i],],byid=TRUE)/(1000*1000))
  # A.obs<-as.numeric(gArea(observations,byid=TRUE)/(1000*1000))
  # A.obs<-A.obs[pos.ID.obs]
  obs.to.ehype$obs.MAF<-MAF[-pos.ID.obs[i]]/(A[-pos.ID.obs[i]]^c2)
  tyrolRtopObj.MAF<-createRtopObject(obs.to.ehype,ehype.catchments,
                                   formulaString = obs.MAF~1,
                                   params = list(gDist = TRUE,debug.level=0,
                                                 rresol = 500,
                                                 #nclus=nr_cores,
                                                 partialOverlap=TRUE,
                                                 nmax=vic)
                                   )

  tyrolRtopObj.MAF<- rtopVariogram(tyrolRtopObj.MAF)
  tyrolRtopObj.MAF<- rtopFitVariogram(tyrolRtopObj.MAF)
  #tyrolRtopObj.MAF<- checkVario(tyrolRtopObj.MAF, cloud = TRUE, identify = TRUE,
  #                              acor = 0.01,log="")
  tyrolRtopObj.MAF<- rtopKrige(tyrolRtopObj.MAF)
  pred.MAF<-tyrolRtopObj.MAF$predictions$var1.pred[pos.ID.ehype[i]]*(A.ehype[i]^c2)

  #### storing the result of FDC prediction in a matrix
est.FDC[,i]<-est.dimless.FDC*pred.MAF
}

#######################################################################################
#######################################################################################
##### Time series manipulation and preparation
A.obs<-as.numeric(gArea(observations,byid=TRUE)/(1000*1000))
A.obs<-A.obs[pos.ID.obs]
streamflow.series<-list()
for (i in 1:l){
  if (cod[pos.ID.obs[i]]%in%Qn$V1){
    dummy.site<-Qn[which(Qn$V1==cod[pos.ID.obs[i]]),2:3]
    dummy.site[,1]<-as.Date(dummy.site[,1])
    rownames(dummy.site)<-1:length(dummy.site[,1])
    colnames(dummy.site)<-c("Date","Q")
    extract.series<-dummy.site
    
  }
  else {
    dummy.site<-Qs[which(Qs$ID==cod[pos.ID.obs[i]]),]
    extract.series<-extract_series(Q=dummy.site)
    negative.values<-which(extract.series[,2]<0)
    extract.series[negative.values,2]<-NA
    }
    
  streamflow.series[[i]]<-extract.series  
}

EHYPE.streamflow.series<-list()
dummy.ehype<-Q.ehype.tyrol[,pos.ID.ehype]
#colnames(EHYPE.streamflow.series)<-ID.ehype
t.ehype<-seq.Date(from=as.Date("1979-01-01","%Y-%m-%d"),to=as.Date("2010-12-31","%Y-%m-%d"),by="days")
for (j in 1:l){
  dummy.ehype[,j]<-dummy.ehype[,j]*A.obs[j]/A.ehype[j]
  EHYPE.streamflow.series[[j]]<-data.frame(Date=t.ehype,Q=dummy.ehype[,j])
  
}

###### assimilation algorithm - povides corrected EHYPE series
######
EHYPE.streamflow.series.corrected<-list()
EHYPE.streamflow.series.cut<-list()
for (j in 1:l){
  date.init<-streamflow.series[[j]]$Date[1]
  date.end<-tail(streamflow.series[[j]]$Date,1)
  pos.dates.init.ehype<-which(EHYPE.streamflow.series[[j]]$Date<date.init)
  pos.dates.end.ehype<-which(EHYPE.streamflow.series[[j]]$Date>date.end)
  dummy.ehype.streamflow<-EHYPE.streamflow.series[[j]]$Q[-c(pos.dates.init.ehype,pos.dates.end.ehype)]
  EHYPE.streamflow.series.cut[[j]]<-data.frame(Date=streamflow.series[[j]]$Date,
                                              Q=dummy.ehype.streamflow)
  complete.ehype.FDC<-sort(dummy.ehype.streamflow,decreasing = TRUE)
  ehype.MAF[j]<-mean(dummy.ehype.streamflow)
  d.ehype<-1:(length(complete.ehype.FDC))/(length(complete.ehype.FDC)+1) 
  ##### resampling FDCs to 20 points  
  resamp<-resample.FDC(complete.ehype.FDC,sam=sam,norm=FALSE)
  ehype.FDC[,j]<-resamp$y
  #### computing error duration curve
  err.FDC[,j]<-est.FDC[,j]-ehype.FDC[,j]
  f_error<-approxfun(sam,err.FDC[,j],rule=2)
  #### computing corrected series
  corrected.series<-c()
  for (i in 1:length(streamflow.series[[j]]$Q)){
    if (is.na(streamflow.series[[j]]$Q[i])){
      corrected.series[i]<-NA
      }
    else {
    #Q_star<-EHYPE.streamflow.series[[j]]$Q[i]
    Q_star<-dummy.ehype.streamflow[i]
    k<-which(complete.ehype.FDC<=Q_star)[1]
    d_star<-d.ehype[k]
    error_star<-f_error(d_star)
    corrected.series[i]<-Q_star+error_star
    }
  }
  
  EHYPE.streamflow.series.corrected[[j]]<-data.frame(Date=streamflow.series[[j]]$Date,
                                                     Q=corrected.series)
  print(j)
}
#############################################################################
#############################################################################
##### Performances analsysis and diagnostics

NSE.global.EHYPE<-c()
NSE.global.EHYPE.corrected<-c()

for (j in 1:l){
  remove.na<-which(is.na(streamflow.series[[j]]$Q))
  if (length(remove.na)>0){
  NSE.global.EHYPE[j]<-1-sum((streamflow.series[[j]]$Q[-remove.na]-EHYPE.streamflow.series.cut[[j]]$Q[-remove.na])^2)/
            sum((streamflow.series[[j]]$Q[-remove.na]-mean(streamflow.series[[j]]$Q[-remove.na]))^2)
  NSE.global.EHYPE.corrected[j]<-1-sum((streamflow.series[[j]]$Q[-remove.na]-EHYPE.streamflow.series.corrected[[j]]$Q[-remove.na])^2)/
    sum((streamflow.series[[j]]$Q[-remove.na]-mean(streamflow.series[[j]]$Q[-remove.na]))^2)
  }
  else{
    NSE.global.EHYPE[j]<-1-sum((streamflow.series[[j]]$Q-EHYPE.streamflow.series.cut[[j]]$Q)^2)/
      sum((streamflow.series[[j]]$Q-mean(streamflow.series[[j]]$Q))^2)
    NSE.global.EHYPE.corrected[j]<-1-sum((streamflow.series[[j]]$Q-EHYPE.streamflow.series.corrected[[j]]$Q)^2)/
      sum((streamflow.series[[j]]$Q-mean(streamflow.series[[j]]$Q))^2)
  }
}


plot(1:11,NSE.global.EHYPE,col="red",
     ylim=c(-3,1),xlab="Catchment ID",ylab="NSE [-]",
     axes=F)
points(1:11,NSE.global.EHYPE.corrected,col="blue")
axis(side = 1,at = 1:11,labels =paste(cod[pos.ID.obs]),cex.axis=0.75,las=2)
axis(side = 2)
abline(h=mean(NSE.global.EHYPE),col="red",lty="dashed")
abline(h=mean(NSE.global.EHYPE.corrected),col="blue",lty="dashed")
title(paste("Avg. NSE ehype =", round(mean(NSE.global.EHYPE),3),"\n",
            "Avg. NSE ehype corrected =", round(mean(NSE.global.EHYPE.corrected),3)))
legend(x="bottomleft",legend = c("EHYPE","cEHYPE","Avg.NSE-EHYPE","Avg.NSE-cEHYPE"),pch=c(1,1,NA,NA),
       lty = c(NA,NA,"dashed","dashed"),col=c("red","blue","red","blue"))
box()
########################################################################
# LNSE
LNSE.global.EHYPE<-c()
LNSE.global.EHYPE.corrected<-c()

# EHYPE.streamflow.series.obs.LNSE<-list()
# EHYPE.streamflow.series.cut.LNSE<-list()
# EHYPE.streamflow.series.obs.LNSE
for (j in 1:l){
  remove.na.zero.ehype<-unique(c(which(streamflow.series[[j]]$Q==0),
                            which(EHYPE.streamflow.series.cut[[j]]$Q==0),
                            which(is.na(streamflow.series[[j]]$Q))))
  remove.na.zero.ehype.corrected<-unique(c(which(streamflow.series[[j]]$Q==0),
                            which(EHYPE.streamflow.series.corrected[[j]]$Q==0),
                            which(is.na(streamflow.series[[j]]$Q))))
  if (length(remove.na.zero.ehype)!=0) {
    LNSE.global.EHYPE[j]<-1-sum((log(streamflow.series[[j]]$Q[-remove.na.zero.ehype])-
                    log(EHYPE.streamflow.series.cut[[j]]$Q[-remove.na.zero.ehype]))^2)/
                   sum((log(streamflow.series[[j]]$Q[-remove.na.zero.ehype])-
                  mean(log(streamflow.series[[j]]$Q[-remove.na.zero.ehype])))^2)
   }
  else{
    LNSE.global.EHYPE[j]<-1-sum((log(streamflow.series[[j]]$Q)-log(EHYPE.streamflow.series.cut[[j]]$Q))^2)/
     sum((log(streamflow.series[[j]]$Q)-mean(log(streamflow.series[[j]]$Q)))^2)
   }
  if (length(remove.na.zero.ehype.corrected)!=0) {
     LNSE.global.EHYPE.corrected[j]<-1-sum((log(streamflow.series[[j]]$Q[-remove.na.zero.ehype.corrected])-
                        log(EHYPE.streamflow.series.corrected[[j]]$Q[-remove.na.zero.ehype.corrected]))^2)/
                        sum((log(streamflow.series[[j]]$Q[-remove.na.zero.ehype.corrected])-
                               mean(log(streamflow.series[[j]]$Q[-remove.na.zero.ehype.corrected])))^2)  
  }
  else {
    LNSE.global.EHYPE.corrected[j]<-1-sum((log(streamflow.series[[j]]$Q)-log(EHYPE.streamflow.series.corrected[[j]]$Q))^2)/
    sum((log(streamflow.series[[j]]$Q)-mean(log(streamflow.series[[j]]$Q)))^2)
  }

}

plot(1:11,LNSE.global.EHYPE,col="red",ylim=c(-3,1),xlab="Catchment ID",
     ylab="LNSE [-]",axes=F)
axis(2)
axis(1,at=1:11,labels = paste(selected.catchments$OBS),
     las=2,cex.axis=0.8)
points(1:11,LNSE.global.EHYPE.corrected,col="blue")
abline(h=mean(LNSE.global.EHYPE),col="red",lty="dashed")
abline(h=mean(LNSE.global.EHYPE.corrected),col="blue",lty="dashed")
title(paste("Avg. LNSE ehype =", round(mean(LNSE.global.EHYPE),3),"\n",
            "Avg. LNSE ehype corrected =", round(mean(LNSE.global.EHYPE.corrected),3)))
legend(x="bottomleft",legend = c("EHYPE","cEHYPE","Avg.LNSE-EHYPE","Avg.LNSE-cEHYPE"),pch=c(1,1,NA,NA),
       lty = c(NA,NA,"dashed","dashed"),col=c("red","blue","red","blue"))

box()

#save.image(file="GEEHP_final_results.RData")
##### plot simulated vs. corrected vs. observed streamflow series 
for (j in 1:l){
  plot(streamflow.series[[j]],type="l",xlab="Days",ylab=expression("Streamflow [m"^3*"/s]"))
  lines(EHYPE.streamflow.series.cut[[j]],col="red")
  lines(EHYPE.streamflow.series.corrected[[j]],col="blue")
  title(paste("ID.obs =",cod[pos.ID.obs[j]],"-","ID.ehype =",ID.ehype[j]))
  png(filename = paste("./pics/comparison_pair_",j,".png",sep=""),
      width = 550, height = 550)
  plot(streamflow.series[[j]]$Q,EHYPE.streamflow.series.cut[[j]]$Q,col="red",
       xlab=expression("Observed streamflow [m"^3*"/s]"),
       ylab=expression("Predicted streamflow [m"^3*"/s]"),
       pch=19,cex=0.5,log="xy",asp=1)
  points(streamflow.series[[j]]$Q,EHYPE.streamflow.series.corrected[[j]]$Q,
         col="blue",pch=19,cex=0.5,log="xy")
  title(paste("ID.obs =",cod[pos.ID.obs[j]],"-","ID.ehype =",ID.ehype[j],"\n",
              "LNSE.EHYPE = ", round(LNSE.global.EHYPE[j],3),"\n",
              "LNSE.EHYPE.corrected = ", round(LNSE.global.EHYPE.corrected[j],3)))
  abline(0,1)
  dev.off()
}

######### year-to-year comparison with focus onto the site with best improvement 
NSE.year.EHYPE<-list()
NSE.year.EHYPE.corrected<-list()
LNSE.year.EHYPE<-list()
LNSE.year.EHYPE.corrected<-list()

for (j in 1:l){
  anni<-unique(format(streamflow.series[[j]]$Date,"%Y"))
  num.anni<-length(unique(format(streamflow.series[[j]]$Date,"%Y")))
  dummy.1<-c()
  dummy.2<-c()
  for (i in 1:num.anni){
    pos.series.per.year<-which(format(streamflow.series[[j]]$Date,"%Y")==anni[i] & !is.na(streamflow.series[[j]]$Q))
    obs.streamflow.per.year<-streamflow.series[[j]]$Q[pos.series.per.year]
    EHYPE.streamflow.per.year<-EHYPE.streamflow.series.cut[[j]]$Q[pos.series.per.year]
    EHYPE.corrected.streamflow.per.year<-EHYPE.streamflow.series.corrected[[j]]$Q[pos.series.per.year]
    dummy.1[i]<-1-sum((obs.streamflow.per.year-EHYPE.streamflow.per.year)^2)/sum((obs.streamflow.per.year-mean(obs.streamflow.per.year))^2)
    dummy.2[i]<-1-sum((obs.streamflow.per.year-EHYPE.corrected.streamflow.per.year)^2)/sum((obs.streamflow.per.year-mean(obs.streamflow.per.year))^2)
  if (j==8){ ##### site with best improvement 
    plot(streamflow.series[[j]][pos.series.per.year,],type="l",ylim=c(0,max(EHYPE.streamflow.series.cut[[j]]$Q)))
    lines(EHYPE.streamflow.series.cut[[j]][pos.series.per.year,],col="red")
    lines(EHYPE.streamflow.series.corrected[[j]][pos.series.per.year,],col="blue")
    legend(x="topleft",legend = paste("yr. =",anni[i],"NSE EHYPE =",round(dummy.1[i],3),"NSE cEHYPE =",round(dummy.2[i],3)),cex=0.75)
    title(paste("ID.obs =",cod[pos.ID.obs[j]],"-","ID.ehype =",ID.ehype[j]))
    
    }
  }
  NSE.year.EHYPE[[j]]<-dummy.1
  NSE.year.EHYPE.corrected[[j]]<-dummy.2
}





























