library("parallel")

screen_strategies<-c(0,1,2,3,4,9)
screen_strategy<-0
load(paste("PSA/PSA_",screen_strategy,"_",1,".Rdata",sep = ""))
for (i in 2:10){
  #name of saved files needed
  load(paste("PSA/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
  psaresults<-rbind(psaresults,results)
}
psaresults<-results

for (j in 2:6){
screen_strategy<-screen_strategies[j]

for (i in 1:10){
  #name of saved files needed
  load(paste("PSA/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
 psaresults<-rbind(psaresults,results)
}
}

names(psaresults)[1] <- 'QALY'
names(psaresults)[2] <- 'Cost'
names(psaresults)[3] <- 'Screens'
names(psaresults)[4:23]<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                        "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                        "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                        'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                        "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                        "PSA_util_1to3","PSA_util_4")
names(psaresults)[24]<-"alternative"

alternatives<-c("noscreening","procas","tertiles","3yr","2yr","fullstrat")
psaresults[,25]<-c(1:length(psaresults[,1]))
names(psaresults)[25]<-"id"


psaresults[,24][psaresults[,24]==9]<-"fullstrat"
psaresults[,24]<-as.factor(psaresults[,24])

save(psaresults,file = paste("PSA/PSA_","psaresults",".Rdata",sep = "")) 

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
              s(PSA_MRI_cdr, by = alternative, bs = "cr") +
              s(PSA_US_cdr, by = alternative, bs = "cr") +
              s(PSA_log_norm_mean, by = alternative, bs = "cr") +
              s(PSA_log_norm_sd, by = alternative, bs = "cr") +
              alternative)
summary(modQ)

save(modQ,file="QALYmodel.Rdata")

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
              s(PSA_MRI_cdr, by = alternative, bs = "cr") +
              s(PSA_US_cdr, by = alternative, bs = "cr") +
              s(PSA_log_norm_mean, by = alternative, bs = "cr") +
              s(PSA_log_norm_sd, by = alternative, bs = "cr") +
              alternative)
summary(modC)