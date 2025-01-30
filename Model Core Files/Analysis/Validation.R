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