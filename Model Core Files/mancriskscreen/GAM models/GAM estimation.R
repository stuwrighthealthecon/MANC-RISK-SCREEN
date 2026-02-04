library("parallel")
library("mgcv")
library("tidyverse")

alternatives <- c(0, 1, 2, 3, 4, 9)

filenames <- list.files("PSA results/misclassification", full.names = TRUE)
alldata <- lapply(filenames, function(x) {
  get(load(x, .GlobalEnv))
})
alldata <- do.call("rbind", alldata)

psaresults <- alldata

#Replace alternative name with string
psaresults$alternative[psaresults$alternative == 0] <- "No Screening"
psaresults$alternative[psaresults$alternative == 1] <- "Risk-1"
psaresults$alternative[psaresults$alternative == 2] <- "Risk-2"
psaresults$alternative[psaresults$alternative == 3] <- "3 yearly"
psaresults$alternative[psaresults$alternative == 4] <- "2 yearly"
psaresults$alternative[psaresults$alternative == 9] <- "Risk-3"

psaresults$alternative <- as.factor(psaresults$alternative)

#Slim down psaresults and garbage clean to save space
psaresults <- psaresults[-c(3, 4, 6, 8, 9, 10, 11, 12, 13)]
rm(results, alldata)
gc()

#Save combined psaresults as backup
save(
  psaresults,
  file = paste("PSA results/PSA_", "psaresults", ".Rdata", sep = "")
)

#Remove cost variables for QALY GAM
load(paste("PSA results/PSA_psaresults", ".Rdata", sep = ""))

psaresultscan <- psaresults %>% filter(psaresults$Cancer == 1)

psaresultscan <- psaresultscan[
  -c(21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36)
]
gc()

rm(psaresults)

#Estimate QALY GAM
modQ <- bam(
  data = psaresultscan,
  discrete = TRUE,
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
    alternative
)
summary(modQ)

#Slim down QALY GAM and save
#modQ[2:43]<-NULL
saveRDS(modQ, file = "GAM models/QALYmodelcan.RDS")

#Re-load full PSA results
load(paste("PSA results/PSA_psaresults", ".Rdata", sep = ""))

psaresultscan <- psaresults %>% filter(psaresults$Cancer == 1)

#Remove QoL data from psaresults
psaresultscan <- psaresultscan[-c(23:24)]
rm(modQ)
gc()

#Estimate cost GAM
modC <- bam(
  data = psaresultscan,
  discrete = TRUE,
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
    alternative
)
summary(modC)

#Slim down cost GAM and save
#modC[2:43]<-NULL
saveRDS(modC, file = "GAM models/costmodelcan.RDS")

load(paste("PSA results/PSA_psaresults", ".Rdata", sep = ""))

psaresultscan <- psaresults %>% filter(psaresults$Cancer == 0)

psaresultscan <- psaresultscan[
  -c(21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36)
]
gc()

rm(psaresults)

#Estimate QALY GAM
modQ <- bam(
  data = psaresultscan,
  discrete = TRUE,
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
    alternative
)
summary(modQ)

#Slim down QALY GAM and save
#modQ[2:43]<-NULL
saveRDS(modQ, file = "GAM models/QALYmodelnocan.RDS")

#Re-load full PSA results
load(paste("PSA results/PSA_psaresults", ".Rdata", sep = ""))

psaresultscan <- psaresults %>% filter(psaresults$Cancer == 0)

#Remove QoL data from psaresults
psaresultscan <- psaresultscan[-c(23:24)]
rm(modQ)
gc()

#Estimate cost GAM
modC <- bam(
  data = psaresultscan,
  discrete = TRUE,
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
    alternative
)
summary(modC)

#Slim down cost GAM and save
#modC[2:43]<-NULL
saveRDS(modC, file = "GAM models/costmodelnocan.RDS")
