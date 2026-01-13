#Stage I to III survival parameters
PSA_gamma_survival1 <- c(-6, -3)
PSA_gamma_survival2 <- c(-5, -1)
PSA_gamma_survival3 <- c(-3, -0.01)

#Metastatic survival parameters
PSA_meta_survival1 <- c(-2.5, -0.5)
PSA_meta_survival2 <- c(-2, -0.3)
PSA_meta_survival3 <- c(-1.5, -0.1)

#Draw Mammography with sensitivity conditional on tumour diameter parameters W-F
PSA_beta1 <- c(1.47, 0.1)
PSA_beta2 <- c(6.51, 0.5)

#Sensitivity by VDG group
PSA_Sen_VDG1 <- c(48, 16)
PSA_Sen_VDG2 <- c(208, 75)
PSA_Sen_VDG3 <- c(113, 76)
PSA_Sen_VDG4 <- c(40, 38)
Sen_VDG_av <- 0.757

#Draw supplemental Screening CDRs
PSA_MRI_cdr <- c(99.495, 19799.5)
PSA_US_cdr <- c(35.89, 11927)

#Draw tumour growth rate parameters
PSA_log_norm_mean <- c(0.8, 1.2)
PSA_log_norm_sd <- c(1.31, 0.11)

# Chemoprevention Drug parameters

# First bring in log hazard ratios from networked analysis
loghaz_ests <- readRDS("Data/PreventionOutputs.RDS")
efficacy_ests <- loghaz_ests[1]
dropout_ests <- loghaz_ests[4]

# Extract parameters for multivariate normal draws
efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

#Drug uptake parameters
PSA_uptake_1 <- c(.71, .1)
PSA_uptake_2 <- c(.71, .1)

#Draw costs
PSA_cost_strat <- c(1.8826894, 0.1015175)
cost_inflator <- 0.2
cost_dist_sd <- 0.1020408
PSA_costvar <- c(0, cost_dist_sd)
PSA_costscreen <- c(0, cost_dist_sd)
PSA_cost_follow_up <- c(0, cost_dist_sd)
PSA_cost_biop <- c(0, cost_dist_sd)
PSA_cost_US <- c(0, cost_dist_sd)
PSA_cost_MRI <- c(0, cost_dist_sd)
PSA_cost_drug <- c(0, cost_dist_sd)

#Generate utility draws
PSA_util_1 <- c(0.6, 0.9)
PSA_util_4 <- c(0.5, 0.8)
