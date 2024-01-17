library("parallel")
library("mgcv")
library("tidyverse")

#Load PSA results for each strategy
screen_strategies<-c(0,1,2,3,4,9)
screen_strategy<-0
load(paste("PSA results/PSA_",screen_strategy,"_",1,".Rdata",sep = ""))
results<-results %>% filter(results[,4]>50 | results[,4]==0)
results<-results[-c(3:4)]
psaresults<-results

for (i in 2:10){
  #name of saved files needed
  load(paste("PSA results/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
  results<-results %>% filter(results[,4]>50 | results[,4]==0)
  results<-results[-c(3:4)]
  psaresults<-rbind(psaresults,results)
}

for (j in 2:6){
screen_strategy<-screen_strategies[j]

for (i in 1:10){
  #name of saved files needed
  load(paste("PSA results/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
  results<-results %>% filter(results[,4]>50 | results[,4]==0)
  results<-results[-c(3:4)]
 psaresults<-rbind(psaresults,results)
}
}

#Replace alternative name with string
psaresults$alternative[psaresults$alternative==0]<-"No Screening"
psaresults$alternative[psaresults$alternative==1]<-"Risk-1"
psaresults$alternative[psaresults$alternative==2]<-"Risk-2"
psaresults$alternative[psaresults$alternative==3]<-"3 yearly"
psaresults$alternative[psaresults$alternative==4]<-"2 yearly"
psaresults$alternative[psaresults$alternative==9]<-"Risk-3"

psaresults$alternative<-as.factor(psaresults$alternative)

#Slim down psaresults and garbage clean to save space
psaresults<-psaresults[-c(3:4)]
psaresults<-psaresults[-c(4:5)]
psaresults<-psaresults[-c(16:17)]
rm(results)
gc()

#Save combined psaresults as backup
save(psaresults,file = paste("PSA results/PSA_","psaresults",".Rdata",sep = "")) 

#Remove cost variables for QALY GAM
psaresults<-psaresults[-c(22:26)]
psaresults<-psaresults[-c(18:19)]
gc()

#Estimate QALY GAM
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

#Slim down QALY GAM and save
modQ[2:43]<-NULL
saveRDS(modQ,file="QALYmodelslim.RDS")

#Re-load full PSA results
load(paste("PSA results/PSA_psaresults",".Rdata",sep = ""))

#Remove QoL data from psaresults
psaresults<-psaresults[-c(20:21)]
gc()

#Estimate cost GAM
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

#Slim down cost GAM and save
modC[2:43]<-NULL
saveRDS(modC,file="costmodelslim.RDS")
