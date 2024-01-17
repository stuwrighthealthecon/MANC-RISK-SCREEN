library("mgcv")
library("dampack")
library(gt)
library(tidyverse)
library(magrittr)
library("MASS")

#Set number of PSA runs to estimate over
mcruns<-1000000

#Set alternatives
alternative<-c(0,1,2,3,4,9)

#Set wtp thresholds to estimate over
wtp<-seq(from=0,to=100000,by=1000)

#Load GAMs
modQ<-readRDS("QALYmodelslim.RDS")
modC<-readRDS("costmodelslim.RDS")

#Draw stage I to III survival parameters
survmvn<-data.frame(c(-5.46208,-5.2077,-5.8016),
                    c(-3.8163,-3.75901,-3.8811),
                    c(-2.72264,-2.66053,-2.78617))
survcovmat<-cov(survmvn)
survmeans<-c(survmvn[1,1],survmvn[1,2],survmvn[1,3])
PSA_gamma_survival<-mvrnorm(mcruns,survmeans,survcovmat)

# Draw Metatstatic survival parameters
metmvn<-data.frame(c(-1.78723,-1.67922,-1.89434),c(-1.38762,-1.33512,-1.49956),c(-1.01051,-0.93338,-1.08304))
metmat<-cov(metmvn)
metmeans<-c(metmvn[1,1],metmvn[1,2],metmvn[1,3])
PSA_meta_survival<-mvrnorm(mcruns,metmeans,metmat)

#Draw Mammography with sensitivity conditional on tumour diameter parameters W-F
PSA_beta1 <- rnorm(mcruns,1.47,0.1)
PSA_beta2 <- rnorm(mcruns,6.51,0.5)

#Draw Mammography sensitivity by volpara density grade from PREVENTICON
PSA_Sen_VDG <- data.frame(rbeta(mcruns,96,16),rbeta(mcruns,298,86),rbeta(mcruns,212,93),rbeta(mcruns,61,39))
Sen_VDG_av <- 0.757

#Draw supplemental Screening CDRs
PSA_MRI_cdr <- rbeta(mcruns,99.495,19799.5) #CDR for MRI in Mammo negative women (incremental)
PSA_US_cdr <- rbeta(mcruns,35.89,11927) #CDR for US in Mammo negative women (incremental)

#Draw tumour growth rate parameters
PSA_log_norm_mean <- rnorm(mcruns,1.07,0.09)
PSA_log_norm_sd <- rnorm(mcruns,1.31,0.11)

#Draw costs
PSA_cost_strat<-(rlnorm(mcruns,2.13387381,0.06349671)*1.0272)
PSA_costvar<-rnorm(mcruns,0,0.1020408)
PSA_costscreen<-rnorm(mcruns,0,0.1020408)
PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
PSA_cost_US<-rnorm(mcruns,0,0.1020408)
PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)

#Generate utility draws
utilmat<-data.frame(c(1-0.82,1-0.81,1-0.83),c(1-0.75,1-0.73,1-0.77))
lnutilmat<-log(utilmat)
covutil<-cov(lnutilmat)
utilmeans<-c(log(1-0.82),log(1-0.75))
PSA_util<-1-exp(mvrnorm(mcruns,utilmeans,covutil))

#Bind monte carlo draws
PSA_all_p<-cbind(PSA_gamma_survival,PSA_meta_survival,PSA_beta1,PSA_beta2,
                 PSA_Sen_VDG,PSA_MRI_cdr,PSA_US_cdr,PSA_log_norm_mean,
                 PSA_log_norm_sd,PSA_cost_strat,PSA_costvar,PSA_util,PSA_costscreen,
                 PSA_cost_follow_up,PSA_cost_biop,PSA_cost_US,PSA_cost_MRI)
PSA_all_p<-as.data.frame(PSA_all_p)
colnames(PSA_all_p)<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                       "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                       "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                       'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                       "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                       "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                       "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI")

alt_names<-c("No Screening","Risk-1","Risk-2","3 yearly","2 yearly","Risk-3")

#Create data.frames to estimate costs and QALYs
output_costs<-data.frame(matrix(nrow=mcruns,ncol=length(alternative)))
output_qalys<-data.frame(matrix(nrow=mcruns,ncol=length(alternative)))
colnames(output_costs)<-alt_names
colnames(output_qalys)<-alt_names

#For each alternative
for (i in 1:length(alternative)){
  PSA_all_p$alternative<-as.factor(alt_names[i])

#Predict QALYs and Costs using GAMS
output_qalys[,i]<-predict.bam(modQ,PSA_all_p)
output_costs[,i]<-predict.bam(modC,PSA_all_p)
}

#Rename columns
alt_names<-c("No Screening","Risk-1","Risk-2","Three yearly","Two yearly","Risk-3")
colnames(output_costs)<-alt_names
colnames(output_qalys)<-alt_names

#Produce PSA object
psa_obj <- make_psa_obj(cost = output_costs,
                        effectiveness = output_qalys,
                        parameters = PSA_all_p,
                        strategies = alt_names,
                        currency = "£")

#Create CEAC object
ceac_obj<- ceac(wtp,psa_obj)

#Plot CEAC
plot(ceac_obj,frontier="FALSE",points="FALSE",xlab="Willingness to Pay (Thousand £ / QALY")

#Creat PSA results table and save
psa_sum <- summary(psa_obj, 
                   calc_sds = TRUE)
psa_sum
write.csv(psa_sum,file="psa results summary.csv")

#Plot cost-effectiveness plane for PSA results
icers <- calculate_icers(cost = psa_sum$meanCost, 
                         effect = psa_sum$meanEffect, 
                         strategies = psa_sum$Strategy)
plot(icers,labels="all")

#One-way sensitivity analysis for costs, effectiveness and net benefit
sacost<-owsa(psa_obj,outcome="cost")
owsa_tornado(sacost,return="plot")

saeff<-owsa(psa_obj,outcome="eff")
owsa_tornado(saeff,return="plot")

sanb<-owsa(psa_obj,outcome="nmb",wtp=20000)
owsa_tornado(sanb,return="plot")

#Plot chart of impact of changes in parameters on optimal strategy
owsa_opt_strat(sanb,return="plot",col="full",plot_const=FALSE)

#Plot ceac omitting infeasible strategies (too high number of scans)
feasqalys<-output_qalys[-c(3,5)]
feascosts<-output_costs[-c(3,5)]
feas_psa_obj <- make_psa_obj(cost = feascosts,
                        effectiveness = feasqalys,
                        parameters = PSA_all_p,
                        strategies = alt_names[-c(3,5)],
                        currency = "£")
feas_ceac_obj<- ceac(wtp,feas_psa_obj)
plot(feas_ceac_obj,frontier="FALSE",points="FALSE",xlab="Willingness to Pay (Thousand £ / QALY")
