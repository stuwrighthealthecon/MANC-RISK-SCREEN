screen_strategy<-3
load(paste("Deterministic results/misclassification_and_preventative_drug/Determ_",screen_strategy,"_",1,".Rdata",sep = ""))

#Remove women with cancer diagnosed before age 50
results<-results %>% filter(results[,4]>50 | results[,4]==0)
detresults<-results

for (i in 2:10){
  #name of saved files needed
  load(paste("Deterministic results/misclassification_and_preventative_drug/Determ_",screen_strategy,"_",i,".Rdata",sep = ""))
  
  #Remove women with cancer diagnosed before age 50
  results<-results %>% filter(results[,4]>50 | results[,4]==0)
  detresults<-rbind(detresults,results)
}

screendet<-detresults %>% filter(detresults$`screen detected`==1)
clindet<-detresults %>% filter(detresults$`screen detected`==0)

tabulate(screendet$Stage)
tabulate(detresults$Stage)

mean(clindet$`Cancer Size`)
quantile(clindet$`Cancer Size`,probs=c(0.25,0.75))
plot(density(clindet$`Cancer Size`))
mean(screendet$`Cancer Size`)
quantile(screendet$`Cancer Size`,probs=c(0.25,0.75))
plot(density(screendet$`Cancer Size`))


cancersbyscreen<-tabulate(screendet$`Cancer Screen Number`)
median(cancersbyscreen)



