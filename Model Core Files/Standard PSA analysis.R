screen_strategies<-c(0,1,2,3,4,9)
screen_strategy<-0
wtp<-seq(from=0,to=100000,by=5000)
alt_names<-c("No Screening","Risk-1","Risk-2","3 yearly","2 yearly","Risk-3")
fullqalyresults<-data.frame(matrix(ncol=6,nrow=2966851))
colnames(fullqalyresults)<-alt_names
fullcostresults<-data.frame(matrix(ncol=6,nrow=2966851))
colnames(fullcostresults)<-alt_names

for (j in 1:6){
  screen_strategy<-screen_strategies[1]
  load(paste("PSA results/PSA_",screen_strategy,"_",1,".Rdata",sep = ""))
  psaqalyresults<-data.frame(results$QALY)
  psacostresults<-data.frame(results$Cost) 
  for (i in 2:10){
    
    #name of saved files needed
    load(paste("PSA results/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
    qalyresults<-data.frame(results$QALY)
    costresults<-data.frame(results$Cost)  
    psaqalyresults<-rbind(psaqalyresults,qalyresults)
    psacostresults<-rbind(psacostresults,costresults)
  }
  fullqalyresults[,j]<-psaqalyresults
  fullcostresults[,j]<-psacostresults
  rm(psaqalyresults)
  rm(psacostresults)
}

#name of saved files needed
load(paste("Risksample/risksample_",1,".Rdata",sep = ""))
names(splitsample) <- sub(paste('^X',"1",".",sep=""), '', names(splitsample))
risksample<-splitsample

for (i in 2:10){
  load(paste("Risksample/risksample_",i,".Rdata",sep = ""))
  names(splitsample) <- sub(paste('^X',i,".",sep=""), '', names(splitsample))
  risksample<-rbind(risksample,splitsample)
}

psa_obj <- make_psa_obj(cost = fullcostresults,
                        effectiveness = fullqalyresults,
                        parameters = risksample,
                        strategies = alt_names,
                        currency = "£")

ceac_obj<- ceac(wtp,psa_obj)
plot(ceac_obj)


