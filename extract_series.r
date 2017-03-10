extract_series<-function(Q){
# Installazione del pacchetto "date"
library(date)
#
#####################################################################################
#######     Importazione dei dati e costruzione del data frame contenente     #######
#######                     date e valori di precipitazione                   #######
#
# Inserire il nome della stazione
#nome.stazione="Anversa"
#nome.stazione="Capistrello"
#
### Importazione dati da file .csv
#dati=read.csv(file=paste(nome.stazione,".csv",sep=""),header=TRUE,sep=",")
dati<-Q
# Numero di giorni del mese
giorni.mese<-c(31,28,31,30,31,30,31,31,30,31,30,31)

# Data frame contenente date e valori di precipitazione
serie=data.frame(Date=as.Date(character()),Q=numeric(),stringsAsFactors=FALSE,rownames=NULL)

# Contatore
conta=1
# Produce i due vettori
for (i in 1:length(dati$ID))
{
  # Riconoscimento anni bisestili
  if (dati$M[i]==2 & ((dati$Y[i]-1904)/4-round((dati$Y[i]-1904)/4,0)==0))
  {
    n.giorni=29  
  }
  else
  {
    n.giorni=giorni.mese[dati$M[i]]}
  for (j in 1:n.giorni)
  {
    # prima colonna: date di accadimento
    serie[conta,1]=as.Date(paste(dati$Y[i],"/",dati$M[i],"/",j,sep = ""),"%Y/%m/%d")
    # seconda colonna: dati di pioggia
    serie[conta,2]=as.numeric(as.vector(dati[i,3+j])) # as.vector forza il dataframe a vettore, poi
    # as.numeric forza il vettore ad uno scalare
    conta=conta+1 # aggiorna il contatore
  }
}
return(serie)
}


