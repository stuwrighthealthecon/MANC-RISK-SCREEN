#Stage I to III survival parameters
survmvn<-data.frame(c(-5.46208,-5.2077,-5.8016),
                    c(-3.8163,-3.75901,-3.8811),
                    c(-2.72264,-2.66053,-2.78617))
survcovmat<-cov(survmvn)
survmeans<-c(survmvn[1,1],survmvn[1,2],survmvn[1,3])

#Metastatic survival parameters
metmvn<-data.frame(c(-1.78723,-1.67922,-1.89434),
                   c(-1.38762,-1.33512,-1.49956),
                   c(-1.01051,-0.93338,-1.08304))
metmat<-cov(metmvn)
metmeans<-c(metmvn[1,1],metmvn[1,2],metmvn[1,3])

#Draw Mammography with sensitivity conditional on tumour diameter parameters W-F
PSA_beta1<-c(1.47,0.1)
PSA_beta2<-c(6.51,0.5)

#Sensitivity by VDG group
PSA_Sen_VDG1<-c(96,16)
PSA_Sen_VDG2<-c(298,86)
PSA_Sen_VDG3<-c(212,93)
PSA_Sen_VDG4<-c(61,39)
Sen_VDG_av<- 0.757

#Draw supplemental Screening CDRs
PSA_MRI_cdr<-c(99.495,19799.5)
PSA_US_cdr<-c(35.89,11927)

#Draw tumour growth rate parameters
PSA_log_norm_mean <- c(1.07,0.09)
PSA_log_norm_sd <- c(1.31,0.11)

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
PSA_uptake_1<-c(.71, .1)
PSA_uptake_2<-c(.71, .1)

#Draw costs
PSA_cost_strat<-c(2.13387381,0.06349671)
cost_inflator<-0.1
cost_dist_sd<-(cost_inflator/1.96)^2
PSA_costvar<-c(0,cost_dist_sd)
PSA_costscreen<-c(0,cost_dist_sd)
PSA_cost_follow_up<-c(0,cost_dist_sd)
PSA_cost_biop<-c(0,cost_dist_sd)
PSA_cost_US<-c(0,cost_dist_sd)
PSA_cost_MRI<-c(0,cost_dist_sd)
PSA_cost_drug<-c(0,cost_dist_sd)

#Generate utility draws
utilmat<-data.frame(c(1-0.82,1-0.81,1-0.83),c(1-0.75,1-0.73,1-0.77))
lnutilmat<-log(utilmat)
covutil<-cov(lnutilmat)
utilmeans<-c(log(1-0.82),log(1-0.75))