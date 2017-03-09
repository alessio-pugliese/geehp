boxplot.w.neg.logy<-function(x,ymin=FALSE,...){
  if (is.null(dim(x))) {
    dimens<-c(0,length(x))
    labels.x <-NA  
  } else {
    dimens<-dim(x)
    labels.x<-colnames(x)
  }
  bp<-boxplot(x,at=seq(1,dimens[2]),plot=FALSE)
  if (ymin==FALSE){
    low.bound<-ceiling(log10(abs(min(bp$out))+1))
    f<-round(min(bp$out),digits=-1)-10    
  } else {
    low.bound<-ceiling(log10(abs(ymin)+1))
    f<-round(ymin,digits=-1)-10
  }

  bxp(z=bp,outline=FALSE,ylim=c(-low.bound,1),axes="FALSE",...)
  box()
  abline(0,0,lty="dashed")
  for (i in 1:dimens[2]){
  outliers<-bp$out[which(bp$group==i)]
  outliers01<-which(outliers>=0)
  neg.outliers<-which(outliers<0)
  points(rep(i,length(outliers01)),outliers[outliers01])
  points(rep(i,length(neg.outliers)),-log10(abs(outliers[neg.outliers])+1))
  }
  axis.ticks.neg<- c(seq(-f,5,-5),1,0.5)
  axis.ticks.pos<-seq(0,1,0.1)
  axis(2,at=c(-log10(axis.ticks.neg+1),axis.ticks.pos),
       labels=paste(c(-axis.ticks.neg,axis.ticks.pos)),cex.axis=0.8,tcl=-0.3)
  axis(1,at=seq(1,dimens[2]),labels=labels.x,cex.axis=0.8,tcl=-0.3)
}