require(rtop)
require(maptools)
require(rgdal)
require(pracma)
require(nortest)
require(rgeos)
library(ggplot2)
library(ggmap)

OBS.boundaries <- readOGR("shp","tyrol")
EHYPE.boundaries <- readOGR("shp","tyrol_ehype")
distance.centroids<-read.csv("shp/distance_matrix.csv",header = TRUE,row.names = 1)
distance.centroids<-distance.centroids/1000


Aobs<-gArea(OBS.boundaries,byid = TRUE)/(1000*1000) ### in km2
Aehype<-gArea(EHYPE.boundaries,byid = TRUE)/(1000*1000) ### in km2

Amatrix<-matrix(NA,nrow=52,ncol=55)
Dmatrix<-matrix(NA,nrow=52,ncol=55)
colnames(Amatrix)<-paste(OBS.boundaries$ID)
rownames(Amatrix)<-paste(EHYPE.boundaries$macroid1)
colnames(Dmatrix)<-paste(OBS.boundaries$ID)
rownames(Dmatrix)<-paste(EHYPE.boundaries$macroid1)
comb.matrix<-matrix(NA,nrow = 52*55,ncol=6)
colnames(comb.matrix)<-c("EHYPE","OBS","A_OBS","A_EHYPE","Diff_A_(%)","Dist_(km)")
k<-1
for (i in 1:52){
  for (j in 1:55){
    if (gOverlaps(OBS.boundaries[j,],EHYPE.boundaries[i,])==TRUE){
    diff_area<-abs((Aobs[j]-Aehype[i]))/Aobs[j]
    #diff_area<-abs(gArea(OBS.boundaries[j,])-gArea(EHYPE.boundaries[i,]))/gArea(OBS.boundaries[j,])
    diff_dist<-gDistance(OBS.boundaries[j,],EHYPE.boundaries[i,],hausdorff = TRUE)/1000
    overlap<-gOverlaps(OBS.boundaries[j,],EHYPE.boundaries[i,])
    Amatrix[i,j]<-diff_area
    Dmatrix[i,j]<-diff_dist
    comb.matrix[k,]<-c(EHYPE.boundaries$macroid1[i],
                          OBS.boundaries$ID[j],round(Aobs[j],1),round(Aehype[i],1),
                          round(diff_area,3),round(diff_dist,1))
        }
    k<-k+1
  }
}

### write the result on a csv file
### in a separte file, we filter the catchements on the basis of distance and areas  

write.csv(comb.matrix,"Code_D_Area_list.csv",row.names = FALSE)

### print on a background googls map EHYPE and obs catchements 

pairs<-read.csv("Code_D_Area_selected.csv")
the_map<-list()
for (i in 1:dim(pairs)[1]){
  OBS.boundaries.latlong<-spTransform(OBS.boundaries[which(OBS.boundaries$ID==pairs$OBS[i]),],CRS("+proj=longlat +datum=WGS84"))
  OBS.boundaries.latlong<-fortify(OBS.boundaries.latlong)
  EHYPE.boundaries.latlong<-spTransform(EHYPE.boundaries[which(EHYPE.boundaries$macroid1==pairs$EHYPE[i]),],CRS("+proj=longlat +datum=WGS84"))
  EHYPE.boundaries.latlong<-fortify(EHYPE.boundaries.latlong)
  box.obs<-make_bbox(lon = OBS.boundaries.latlong$long,lat = OBS.boundaries.latlong$lat,f = 0.2)
  box.ehype<-make_bbox(lon = EHYPE.boundaries.latlong$long,lat = EHYPE.boundaries.latlong$lat,f = 0.2)
  shift<-0.25
  maxNE<-c(max(box.ehype[3],box.obs[3])+shift,max(box.ehype[4],box.obs[4])+shift)
  minSO<-c(min(box.ehype[1],box.obs[1])-shift,min(box.ehype[2],box.obs[2])-shift)
  mybox<-c(minSO,maxNE)
  map<-get_map(location = mybox, maptype = "satellite",source = "google")
  the_map[[i]]<-ggmap(map,darken = c(0.2,"white"),legend="topright")+
  geom_polygon(aes(x=long,y=lat,group=group),fill=NA,size=1.2,color="black",data=OBS.boundaries.latlong)+
  geom_polygon(aes(x=long,y=lat,group=group),fill=NA,size=1.2,color="red",data=EHYPE.boundaries.latlong)+
  ggtitle(paste(i," - ","EHYPE:",pairs$EHYPE[i]," - ","Meas. Catch.:",pairs$OBS[i]," - ","Diff. Area:",pairs$Diff_A_...[i]*100,"%"))
}
pdf("catch_pairs.pdf",paper = "a4")
the_map
dev.off()

 




  
