install.packages("MASS")

library("MASS")

#Define number of draws
mcruns<-1000

#Set distributions for breast cancer survival by stage
survmvn<-data.frame(c(-5.46208,-5.2077,-5.8016),c(-3.8163,-3.75901,-3.8811),c(-2.72264,-2.66053,-2.78617))
survcovmat<-cov(survmvn)
survmeans<-c(survmvn[1,1],survmvn[1,2],survmvn[1,3])
PSA_gamma_survival<-mvrnorm(mcruns,survmeans,survcovmat)

#Metatstatic survival parameters
metmvn<-data.frame(c(-1.78723,-1.67922,-1.89434),c(-1.38762,-1.33512,-1.49956),c(-1.01051,-0.93338,-1.08304))
metmat<-cov(metmvn)
metmeans<-c(metmvn[1,1],metmvn[1,2],metmvn[1,3])
PSA_meta_survival<-mvrnorm(mcruns,metmeans,metmat)

#Mammography with sensitivity conditional on tumour diameter parameters W-F
PSA_beta1 <- rnorm(mcruns,1.47,0.1)
PSA_beta2 <- rnorm(mcruns,6.51,0.5)

#Mammography sensitivity by volpara density grade from PREVENTICON
PSA_Sen_VDG <- data.frame(rbeta(mcruns,96,16),rbeta(mcruns,298,86),rbeta(mcruns,212,93),rbeta(mcruns,61,39))
Sen_VDG_av <- 0.757

#Supplemental Screening
PSA_MRI_cdr <- rbeta(mcruns,99.495,19799.5) #CDR for MRI in Mammo negative women (incremental)
PSA_US_cdr <- rbeta(mcruns,35.89,11927) #CDR for US in Mammo negative women (incremental)

#Set tumour growth rate parameters
PSA_log_norm_mean <- rnorm(mcruns,1.07,0.09)
PSA_log_norm_sd <- rnorm(mcruns,1.31,0.11)

#Costs
PSA_cost_strat<-rlnorm(mcruns,2.36932871,0.05601143)
PSA_costvar<-rnorm(mcruns,0,0.1020408)

#Generate utility draws
utilmat<-data.frame(c(1-0.82,1-0.81,1-0.83),c(1-0.75,1-0.73,1-0.77))
lnutilmat<-log(utilmat)
covutil<-cov(lnutilmat)
utilmeans<-c(log(1-0.82),log(1-0.75))
PSA_util<-1-exp(mvrnorm(mcruns,utilmeans,covutil))

PSA_all_p<-cbind(PSA_gamma_survival,PSA_meta_survival,PSA_beta1,PSA_beta2,
                 PSA_Sen_VDG,PSA_MRI_cdr,PSA_US_cdr,PSA_log_norm_mean,
                 PSA_log_norm_sd,PSA_cost_strat,PSA_costvar,PSA_util)
PSA_all_p<-as.data.frame(PSA_all_p)
colnames(PSA_all_p)<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                       "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                       "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                       'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                       "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                       "PSA_util_1to3","PSA_util_4")

save(PSA_all_p,file="PSA_values.RData")
