library("parallel")
library("mgcv")
library("tidyverse")

screen_strategies<-c(0,1,2,3,4,9)
screen_strategy<-0
load(paste("PSA/PSA_",screen_strategy,"_",1,".Rdata",sep = ""))
results<-results %>% filter(results[,4]>50 | results[,4]==0)
results<-results[-c(3:4)]
psaresults<-results

for (i in 2:10){
  #name of saved files needed
  load(paste("PSA/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
  results<-results %>% filter(results[,4]>50 | results[,4]==0)
  results<-results[-c(3:4)]
  psaresults<-rbind(psaresults,results)
}

for (j in 2:6){
screen_strategy<-screen_strategies[j]

for (i in 1:10){
  #name of saved files needed
  load(paste("PSA/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
  results<-results %>% filter(results[,4]>50 | results[,4]==0)
  results<-results[-c(3:4)]
 psaresults<-rbind(psaresults,results)
}
}


psaresults[,29][psaresults[,29]==0]<-"noscreening"
psaresults[,29][psaresults[,29]==1]<-"procas"
psaresults[,29][psaresults[,29]==2]<-"tertiles"
psaresults[,29][psaresults[,29]==3]<-"3yr"
psaresults[,29][psaresults[,29]==4]<-"2yr"
psaresults[,29][psaresults[,29]==9]<-"fullstrat"

psaresults[,29]<-as.factor(psaresults[,29])

save(psaresults,file = paste("PSA/PSA_","psaresults",".Rdata",sep = "")) 
psaresults<-psaresults[-c(15,16)]

modQ <- bam(data = psaresults,
            formula = QALY ~
              s(PSA_util_1to3, by = alternative, bs = "cr") +
              s(PSA_util_4, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_1, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_2, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_3, by = alternative, bs = "cr") +
              s(PSA_meta_survival_54, by = alternative, bs = "cr") +
              s(PSA_meta_survival_74, by = alternative, bs = "cr") +
              s(PSA_meta_survival_99, by = alternative, bs = "cr") +
              s(PSA_beta_1, by = alternative, bs = "cr") +
              s(PSA_beta_2, by = alternative, bs = "cr") +
              s(PSA_VDG1_sen, by = alternative, bs = "cr") +
              s(PSA_VDG2_sen, by = alternative, bs = "cr") +
              s(PSA_VDG3_sen, by = alternative, bs = "cr") +
              s(PSA_VDG4_sen, by = alternative, bs = "cr") +
              s(PSA_log_norm_mean, by = alternative, bs = "cr") +
              s(PSA_log_norm_sd, by = alternative, bs = "cr") +
              alternative)
summary(modQ)

modQ[2:43]<-NULL
saveRDS(modQ,file="QALYmodelslim.RDS")

modC <- bam(data = psaresults,
            formula = Cost ~ 
              s(PSA_cost_strat, by = alternative, bs = "cr") +
              s(PSA_costvar, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_1, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_2, by = alternative, bs = "cr") +
              s(PSA_gamma_survival_3, by = alternative, bs = "cr") +
              s(PSA_meta_survival_54, by = alternative, bs = "cr") +
              s(PSA_meta_survival_74, by = alternative, bs = "cr") +
              s(PSA_meta_survival_99, by = alternative, bs = "cr") +
              s(PSA_beta_1, by = alternative, bs = "cr") +
              s(PSA_beta_2, by = alternative, bs = "cr") +
              s(PSA_VDG1_sen, by = alternative, bs = "cr") +
              s(PSA_VDG2_sen, by = alternative, bs = "cr") +
              s(PSA_VDG3_sen, by = alternative, bs = "cr") +
              s(PSA_VDG4_sen, by = alternative, bs = "cr") +
              s(PSA_log_norm_mean, by = alternative, bs = "cr") +
              s(PSA_log_norm_sd, by = alternative, bs = "cr") +
              s(PSA_cost_follow_up, by = alternative, bs = "cr") +
              s(PSA_cost_biop, by = alternative, bs = "cr") +
              s(PSA_costscreen, by = alternative, bs = "cr") +
              alternative)
summary(modC)

modC[2:43]<-NULL
saveRDS(modC,file="costmodelslim.RDS")
