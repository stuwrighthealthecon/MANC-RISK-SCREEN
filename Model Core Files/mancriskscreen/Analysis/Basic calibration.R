calibruns<-50

calibresults<-as.data.frame(matrix(nrow=calibruns,ncol=6))
colnames(calibresults)<-c("c","Mean","SD","P","LCbump","RMSE")

plotlist<-list()

for (c in 1:calibruns){

  source("MANC-RISK-SCREEN fast parallelised.R")
  source("Analysis/Validation.R")
  calibresults[c,]<-c(c,gamma_mean,gamma_sd,gamma_p,lc_bump,RMSE)
  
plotlist[[c]]<-incidenceplot
print(paste(c," of ",calibruns," complete"))
}

calibresults %>% arrange(RMSE)

plotlist[[7]]

plot(density(risksample$liferisk),main="Distribution of Lifetime Breast Cancer Risk at Birth")
count(((risksample$liferisk-2.2-lc_bump)>40))/nrow(risksample)
