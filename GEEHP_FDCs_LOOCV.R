#### leave-one-out cross-validation of flow duration curves @ Tyrol dataset
####
#loocv_TNDTK_tyrol<-function(vic){
source('tnd.R')
source('my_SEAS_FDC.r')
source('resample_FDC.r')
source('boxplot_w_negatives_logy.R')
source('lseq.R')
##### packages 
require(rtop)
require(maptools)
require(rgdal)
require(pracma)
require(nortest)
require(rgeos)
#### Load shp files layer for Top-kriging 
observations <- readOGR("shp","tyrol")
#predictionLocations <- readOGR("Dataset_Tyrol", "eHYPE_subbasin_TYROL")

#### FDC construction from daily streamflow data
Qn<-read.table("data_at.asc")
Qs<-read.table("Qgiorn_AltoAdige.txt",header = TRUE,sep=",")

#### gauges selected after GB and JP evaluations 

gauges.removed<-c(201889,201780,201772,201681,201525,201194,201178) ### flagged with red
gauges.removed<-c(gauges.removed,201699,201715) #### flagged as influenced by small headwater dams 
#### gauges removal  
for (i in 1:length(gauges.removed)){
  Qn<-Qn[-which(Qn[,1]==gauges.removed[i]),] ### from dataset
  observations<-observations[-which(observations$ID==gauges.removed[i]),] ### baoundaries removal from shapefile
}

nt.length<-length(unique(Qn[,1]))
st.length<-length(observations$ID)-nt.length

# # Tolgo inoltre 4 (201699) e 10 (201715)
# Qn=Qn[-which(Qn[,1]==201699),]
# Qn=Qn[-which(Qn[,1]==201715),]
# #
##################################
### POR-FDCs builder
S<-data.frame(month=1:12,season=c(1,1,1,1,1,1,1,1,1,1,1,1))
cod<-observations$ID



cod.nt<-observations$ID[1:nt.length]
cod.st<-observations$ID[(nt.length+1):end(observations$ID)[1]]
FDC.st<-seasonal.FDC(S,Qs,cod.st)
FDC.nt<-list(Qs=list(),Y=list())
for (i in 1:length(cod.nt)){
  dummy<-which(Qn[,1]==cod.nt[i])
  FDC.nt$Qs[[i]]<-sort(Qn[dummy,3],decreasing=TRUE)
  #FDC.nt$Y[[i]]<-
}  
FDC<-list(Q=list(),d=list())
years<-c()
for (i in 1:length(cod)){
  #if (i<=32){
  if (i<=nt.length){
    FDC$Q[[i]]<-FDC.nt$Qs[[i]]
    FDC$d[[i]]<-1:length(FDC$Q[[i]])/length(FDC$Q[[i]]+1) #### weibull plotting position
    years[i]<-30
  } else {
    #FDC[[i]]<-FDC.st$Qs[[i-32]][[1]]
    FDC$Q[[i]]<-FDC.st$Qs[[i-nt.length]][[1]]
    FDC$d[[i]]<-1:length(FDC$Q[[i]])/length(FDC$Q[[i]]+1)
    years[i]<-length(FDC.st$Y[[i-nt.length]])
  }
}
#############################################
#############################################
###### FDCs plot
plot(FDC$d[[1]],FDC$Q[[1]]/mean(FDC$Q[[1]]), type="l",log="y",ylim=c(0.001,100),xlim=c(0,1),
     xlab="Duration [-]", ylab="Dimensionless streamflow [m3/s]",col="gray60",axes=FALSE)
axis(1,at=seq(0,1,0.1),labels = paste(seq(0,1,0.1)))
axis(2,at=lseq(0.01,100,5),labels = paste(lseq(0.01,100,5)))
for (i in 2:length(cod)){
  lines(FDC$d[[i]],FDC$Q[[i]]/mean(FDC$Q[[i]]),col="gray60")
}
box()

##############################################
###### MAF 
MAF<-c()
##### Mean annual flow 
for (i in 1:length(cod)){
  MAF[i]<-mean(FDC$Q[[i]])
}
#### analysis of the scaling relationship between MAF and drainage area 
A<-as.numeric(gArea(observations,byid=TRUE)/(1000*1000))
plot(MAF~A,log="xy",xlab="Drainage Area [km^2]",ylab="MAF [m^3/s]",
     xlim=c(1,100000),ylim=c(0.1,1000))
loglin.mod<-lm(log(MAF)~log(A))
c1<-exp(loglin.mod$coefficients[1])
c2<-loglin.mod$coefficients[2]
curve(c1*x^c2,add=TRUE)
text(50,100,paste("R2=",round(summary(loglin.mod)$r.squared,digits=3)))
text(10000,1,paste("MAF=",round(c1,digit=2),"*A^",round(c2,digit=2)))
#################################################
### Minimum record lenght for maximum duration to be used in TND computatins

min.year<-min(years)
maxd<-min.year*365/(min.year*365+1) #limite inferiore delle durate massime

#### TND computed using a gaussian transformation of the duration axis
TND<-c()
for (i in 1:length(cod)){
  TND[i]<-fdc.tnd(FDC$Q[[i]]/MAF[i],norm=TRUE,maxd=maxd)
}
###### Hypotheses test for normality of TNDs
ad.test(TND)
hist(TND,freq=FALSE,ylim=c(0,1),xlim=c(1,4.5))
lines(density(TND))
curve(dnorm(x,mean=mean(TND),sd=sd(TND)),add=TRUE,col="red")
qqnorm(TND)
qqline(TND)

hist((TND-mean(TND))/sd(TND),freq=FALSE)
lines(density((TND-mean(TND))/sd(TND)))
curve(dnorm,from=-4,to=4,add=TRUE,col="red")
qqnorm((TND-mean(TND))/sd(TND))
qqline((TND-mean(TND))/sd(TND))
###########################################

#### Top-kriging TND through LOOCV
set.seed(1)
observations$obs<-TND
#vic<-vic
vic<-6
rtopObj<-createRtopObject(observations, observations,
                          formulaString = obs~1,
                          params = list(gDist = TRUE, rresol = 500,wlim=1,
                                        nmax=vic))

rtopObj <- rtopVariogram(rtopObj) #crea l'oggetto variogramma
rtopObj <- rtopFitVariogram(rtopObj) ## fitta il variogramma teorico 
#rtopObj <- checkVario(rtopObj) 
rtopObj <- rtopKrige(rtopObj,cv=TRUE,wret=TRUE)
pred.TND<-rtopObj$predictions$var1.pred ### TND stimati con il TK
NSE.tnd<-1-sum((observations$obs-pred.TND)^2)/sum((observations$obs-mean(observations$obs))^2)
weights<-rtopObj$weight  #### geostatistical TK weights 

########################################################
###### plotting empirical vs predicted TND values 
plot(observations$obs,rtopObj$predictions$var1.pred,xlim=c(1,4),ylim=c(1,4),
     xlab="Empirical TND [-]",ylab="Predicted TND [-]",pch=19, cex=1.2)
abline(0,1)
text(1.5,3.5,paste("NSE=",round(NSE.tnd,digit=2)),cex=1.2)
text(2,3.7, "Reference standardization value Q*=MAF")
title(paste("LOOCV TNDTK vic=",vic,sep=""))

#### Top-kriging MAF through LOOCV 
observations$obs<-MAF/(A^c2)
rtopObj.MAF<-createRtopObject(observations, observations,
                              formulaString = obs~1,
                              params = list(gDist = TRUE, rresol = 500,
                                            nmax=vic))
rtopObj.MAF <- rtopVariogram(rtopObj.MAF)
rtopObj.MAF <- rtopFitVariogram(rtopObj.MAF)
#rtopObj.MAF <- checkVario(rtopObj.MAF, cloud = TRUE, identify = TRUE,acor = 0.01,log="")
rtopObj.MAF <- rtopKrige(rtopObj.MAF,cv=TRUE)
pred.MAFTK<-rtopObj.MAF$predictions$var1.pred*A^c2
NSE.MAF<-1-sum((MAF-pred.MAFTK)^2)/sum((MAF-mean(MAF))^2)
LNSE.MAF<-1-sum((log(MAF)-log(pred.MAFTK))^2)/sum((log(MAF)-mean(log(MAF)))^2)

### plotting empirical vs predicted MAF
plot(MAF,pred.MAFTK,xlim=c(0.5,500),ylim=c(0.5,500),log="xy",
     xlab="Empirical MAF [m3/s]",ylab="Predicted MAF [m3/s]",pch=19, cex=1.2)
abline(0,1)
text(2,250,paste("NSE=",round(NSE.MAF,digit=2)),cex=1.2)
text(2,150,paste("LNSE=",round(LNSE.MAF,digit=2)),cex=1.2)
title(paste("LOOCV MAFTK vic=",vic,sep=""))

### plotting empirical vs predicted TND
plot(TND,pred.TND,xlim=c(1,4),ylim=c(1,4),
     xlab="Empirical TND [-]",ylab="Predicted TND [-]",pch=4)
text(3.5,1.5,paste("NSE=",round(NSE.tnd,digit=2)),cex=1.2)
abline(0,1)
title("LOOCV TND")



##### Cross-validation of dimensionless FDC 
np<-20 
sam<-pnorm(seq(qnorm(1-maxd),qnorm(maxd),(qnorm(maxd)-qnorm(1-maxd))/(np-1)))
y<-matrix(NA,nrow=np,ncol=length(cod)) ## portate empiriche adimensionali   
Y<-matrix(NA,nrow=np,ncol=length(cod)) ## portate empiriche dimensionali
est<-matrix(NA,nrow=np,ncol=length(cod)) ## portate stimate adimensionali
EST<-matrix(NA,nrow=np,ncol=length(cod)) ## portate stimate adimensionali
for(j in 1:length(cod)){
  resamp<-resample.FDC(FDC$Q[[j]]/MAF[j],sam=sam,norm=FALSE)  
  y[,j]<-resamp$y
}
est<-y%*%t(weights)  ### estimates
NSE.FDC.site<-c()
LNSE.FDC.site<-c()
NSE.FDC.duration<-c()
LNSE.FDC.duration<-c()
for (i in 1:length(cod)){
  EST[,i]<-est[,i]*pred.MAFTK[i]
  Y[,i]<-y[,i]*MAF[i]
  NSE.FDC.site[i]<-1-sum((Y[,i]-EST[,i])^2)/sum((Y[,i]-mean(Y[,i]))^2)
  zeros<-which(as.vector(Y[,i])<=0 | as.vector(EST[,i])<=0)
  if (length(zeros)!=0){
    logy<-log(as.vector(Y[-zeros,i]))
    logest<-log(as.vector(EST[-zeros,i]))
  } else {
    logy<-log(Y[,i])
    logest<-log(EST[,i])
  }
  LNSE.FDC.site[i]<- 1-sum((logy-logest)^2)/sum((logy-mean(logy))^2)
}

for (i in 1:np){
  NSE.FDC.duration[i]<- 1-sum((Y[i,]-EST[i,])^2)/sum((Y[i,]-mean(Y[i,]))^2)
  zeros<-which(as.vector(Y[i,])<=0 | as.vector(EST[i,])<=0)
  if (length(zeros)!=0){
    logy<-log(as.vector(Y[i,-zeros]))
    logest<-log(as.vector(EST[i,-zeros]))
  } else {
    logy<-log(Y[i,])
    logest<-log(EST[i,])
  }
  LNSE.FDC.duration[i]<- 1-sum((logy-logest)^2)/sum((logy-mean(logy))^2)
}

def.par <- par(no.readonly = TRUE) # save default, for resetting...
layout(mat=matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(8,4),heights = c(5,5))
par(mar=c(5,4,4,2),oma=c(3,3,3,3),mgp=c(2.5,0.75,0))
plot(Y,EST,log="xy",xlim=c(0.01,1500),ylim=c(0.01,1500),cex=0.7,pch=19,
     xlab=expression("Empirical streamflow [m"^3*"/s]"),ylab=expression("Predicted streamflow [m"^3*"/s]"),
     xaxt="n",yaxt="n")
axis(side=1,at=c(0.01,1,100),c("0.01","1","100"))
axis(side=2,at=c(0.01,1,100),c("0.01","1","100"))
abline(0,1)
mtext(paste("LOOCV TNDTK vic=",vic,sep=""),side=3,line = 1)
text(paste("Avg. NSE = ",round(mean(NSE.FDC.site),digits=3)),x=0.1,y=100)
text(paste("Avg. LNSE = ",round(mean(LNSE.FDC.site),digits=3)),x=0.1,y=50)

NSE.matrix<-data.frame(NSE=NSE.FDC.site,LNSE=LNSE.FDC.site)
boxplot.w.neg.logy(NSE.matrix,ymin=-5,boxwex=c(0.5,0.5),ylab="Nash-Sutcliffe Efficiencies [-]")

plot(sam,NSE.FDC.duration,type="l",lwd=3,ylim=c(0.5,1),xlab="Duration [-]",ylab = "Nash-Sutcliffe Efficiencies [-]")
lines(sam,LNSE.FDC.duration,lty="dashed",lwd=3)
legend(x="bottomleft",legend=c("NSE","LNSE"),
       lty=c(1,2),lwd=c(3,3),cex=0.8,seg.len = 6,
       bty = "n",y.intersp = 0.05)

mtext(paste("Number of gauged catchments = ",nt.length+st.length,sep=""),side=3,line = 1,outer=TRUE,cex=1.5)

par(def.par)










