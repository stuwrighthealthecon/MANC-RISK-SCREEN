create_sample<-function(PSA=0,intervals=0,seed=1,screen_strategy){

  #Import synthetic dataset derived from PROCAS2 study
  risk_mat<-read.csv("Data/synthetic_risk_data.csv")[,2:4]
  colnames(risk_mat)<-c("VBD","tenyrrisk","liferisk")
  
  #Creat column for risk group to be entered
  risk_mat$risk_group<-numeric(length(risk_mat$VBD))
  
  #Set VDG based on breast density
  risk_mat$VDG<-1+findInterval(risk_mat[,1],VDG_interval)
  
  #Breast density cut-offs for supplemental screening
  density_cutoff <-3
  
  #Set up data frame of women's lifetimes to simulate
  risksample<-risk_mat[sample(nrow(risk_mat),inum*ifelse(PSA==1,mcruns,1),replace=TRUE),]
  risksample[,c("MRI_screen","US_screen","risk_predicted","feedback",
                "interval_change","life_expectancy","cancer",
                "clinical_detect_size",
                "growth_rate")]<-numeric(length=length(risksample$VBD))
  
  #If risk-stratified screening used then determine if each woman chooses to have
  #risk predicted, attends risk consultation, and changes interval
  
    risksample$risk_predicted<-rbinom(length(nrow(risksample)),1,risk_uptake)
    risksample$feedback<-ifelse(risksample$risk_predicted==1 & 
                                   rbinom(length(nrow(risksample)),1,(risk_feedback))==1,1,0)
    risksample$interval_change<-ifelse(risksample$feedback==1 & 
                                          rbinom(length(nrow(risksample)),1,risk_feedback)==1,1,0)
  
  
  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy<- rweibull(n = length(risksample$life_expectancy),
                                        shape = acmmortality_wb_a, 
                                        scale = acmmortality_wb_b)
  risksample$life_expectancy<-ifelse(risksample$life_expectancy<=start_age,
                                     qweibull(p = dqrunif(n = 1,
                min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),
                max = 1),
    shape = acmmortality_wb_a, scale = acmmortality_wb_b),risksample$life_expectancy)
  risksample$life_expectancy<-ifelse(risksample$life_expectancy>=rep(100,length=length(risksample$life_expectancy)),
                                     100,
                                     risksample$life_expectancy)
  
  #Determine if a cancer will develop
  risksample$cancer<-ifelse(dqrunif(length(risksample$cancer),0,1)
                                <(risksample$liferisk/100),1,0)
  
  #Set clinical detection size for cancer
  risksample$clinical_detect_size<-risksample$cancer*(dqrnorm(
    n = length(risksample$cancer),
    mean = clin_detection_m,sd = clin_detection_sd))
  
  #Set growth rate
  risksample$growth_rate<-risksample$cancer*qlnorm(
    dqrunif(length(risksample$cancer),0,1),
    meanlog=log_norm_mean,sdlog=sqrt(log_norm_sd))
  
  if(PSA==1){
  if(intervals==0){
  #########################Add Monte Carlo Draws into Sample##############  
  
  #Draw stage I to III survival parameters
  survmvn<-data.frame(c(-5.46208,-5.2077,-5.8016),
                      c(-3.8163,-3.75901,-3.8811),
                      c(-2.72264,-2.66053,-2.78617))
  survcovmat<-cov(survmvn)
  survmeans<-c(survmvn[1,1],survmvn[1,2],survmvn[1,3])
  PSA_gamma_survival<-mvrnorm(mcruns,survmeans,survcovmat) %>%
    data.frame() %>%
    transpose()
  
  # Draw Metatstatic survival parameters
  metmvn<-data.frame(c(-1.78723,-1.67922,-1.89434),c(-1.38762,-1.33512,-1.49956),c(-1.01051,-0.93338,-1.08304))
  metmat<-cov(metmvn)
  metmeans<-c(metmvn[1,1],metmvn[1,2],metmvn[1,3])
  PSA_meta_survival<-mvrnorm(mcruns,metmeans,metmat) %>%
                     data.frame() %>%
                     transpose()
  
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
  
  #Drug parameters
  
  # First bring in log hazard ratios from networked analysis
  loghaz_ests <- readRDS("PreventionOutputs.RDS")
  efficacy_ests <- loghaz_ests[1]
  dropout_ests <- loghaz_ests[4]
  
  # Extract parameters for multivariate normal draws
  efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
  efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
  dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
  dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()
  
  PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
            data.frame() %>%
            transpose()
  PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
            data.frame() %>%
            transpose()
  
  PSA_uptake_1 <- rnorm(mcruns, .71, .1)
  PSA_uptake_2 <- rnorm(mcruns, .71, .1)
  
  #Draw costs
  PSA_cost_strat<-(rlnorm(mcruns,2.13387381,0.06349671)*1.0272)
  PSA_costvar<-rnorm(mcruns,0,0.1020408)
  PSA_costscreen<-rnorm(mcruns,0,0.1020408)
  PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
  PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
  PSA_cost_US<-rnorm(mcruns,0,0.1020408)
  PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
  PSA_cost_drug <- rnorm(mcruns,0,0.1020408)
  
  #Generate utility draws
  utilmat<-data.frame(c(1-0.82,1-0.81,1-0.83),c(1-0.75,1-0.73,1-0.77))
  lnutilmat<-log(utilmat)
  covutil<-cov(lnutilmat)
  utilmeans<-c(log(1-0.82),log(1-0.75))
  PSA_util<-(1-exp(mvrnorm(mcruns,utilmeans,covutil))) %>%
            data.frame() %>%
            transpose()
  
  }else{
    
  PSA_gamma_survival<-data.frame(runif(mcruns,-6,-3),runif(mcruns,-5,-1),runif(mcruns,-3,-0.01))
  PSA_meta_survival<-data.frame(runif(mcruns,-2.5,-0.5),runif(mcruns,-2,-0.3),runif(mcruns,-1.5,-0.1))
    
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
  PSA_log_norm_mean <- runif(mcruns,0.8,1.2)
  PSA_log_norm_sd <- rnorm(mcruns,1.31,0.11)
  
  #Drug parameters
  
  # First bring in log hazard ratios from networked analysis
  loghaz_ests <- readRDS("PreventionOutputs.RDS")
  efficacy_ests <- loghaz_ests[1]
  dropout_ests <- loghaz_ests[4]
  
  # Extract parameters for multivariate normal draws
  efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
  efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
  dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
  dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()
  
  PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
            data.frame() %>%
            transpose()
  PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
            data.frame() %>%
            transpose()
  
  PSA_uptake_1 <- rnorm(mcruns, .71, .1)
  PSA_uptake_2 <- rnorm(mcruns, .71, .1)
    
  #Costs
  PSA_cost_strat<-rlnorm(mcruns,2.13387381,0.06349671)
  PSA_costvar<-rnorm(mcruns,0,0.1020408)
  PSA_costscreen<-rnorm(mcruns,0,0.1020408)
  PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
  PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
  PSA_cost_US<-rnorm(mcruns,0,0.1020408)
  PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
  PSA_cost_drug <- rnorm(mcruns,0,0.1020408)
    
  PSA_util<-data.frame(runif(mcruns,0.6,0.9),runif(mcruns,0.5,0.8))
  }
  
  #Generate id for monte carlo set
  mcid<-c(1:mcruns)
  
  #Bind monte carlo draws
  PSA_all_p<-cbind(PSA_gamma_survival,
                   PSA_meta_survival,
                   PSA_beta1,
                   PSA_beta2,
                   PSA_Sen_VDG,
                   PSA_MRI_cdr,
                   PSA_US_cdr,
                   PSA_log_norm_mean,
                   PSA_log_norm_sd,
                   PSA_eff,
                   PSA_dropout,
                   PSA_uptake_1,
                   PSA_uptake_2,
                   PSA_cost_strat,
                   PSA_costvar,
                   PSA_util,
                   PSA_costscreen,
                   PSA_cost_follow_up,
                   PSA_cost_biop,
                   PSA_cost_US,
                   PSA_cost_MRI,
                   PSA_cost_drug,
                   mcid)
  PSA_all_p<-as.data.frame(PSA_all_p)
  colnames(PSA_all_p)<-c("PSA_gamma_survival_1",
                          "PSA_gamma_survival_2",
                          "PSA_gamma_survival_3",
                          "PSA_meta_survival_54",
                          "PSA_meta_survival_74",
                          "PSA_meta_survival_99",
                          "PSA_beta_1",
                          "PSA_beta_2",
                          'PSA_VDG1_sen',
                          'PSA_VDG2_sen',
                          'PSA_VDG3_sen',
                          'PSA_VDG4_sen',
                          "PSA_MRI_cdr",
                          "PSA_US_cdr",
                          "PSA_log_norm_mean",
                          "PSA_log_norm_sd",
                          "PSA_eff_ana",
                          "PSA_eff_tam",
                          "PSA_dropout_ana",
                          "PSA_dropout_tam",
                          "PSA_uptake_1",
                          "PSA_uptake_2",
                          "PSA_cost_strat",
                          "PSA_costvar",
                          "PSA_util_1to3",
                          "PSA_util_4",
                          "PSA_costscreen",
                          "PSA_cost_follow_up",
                          "PSA_cost_biop",
                          "PSA_cost_US",
                          "PSA_cost_MRI",
                          "PSA_cost_drug",
                          "mcid")
  
  #Bind individual level parameters and monte carlo draws
  masterframe<-data.frame(matrix(nrow=inum*mcruns,ncol=length(risksample[1,])+length(PSA_all_p[1,])))
  masterframe[,1:14]<-risksample
  masterframe[,15:47]<-PSA_all_p
  colnames(masterframe)[1:14]<-colnames(risksample)
  colnames(masterframe)[15:47]<-colnames(PSA_all_p)
  
  #Split the dataframe into chunks for easier computation
  masterframe$split<-(rep(1:chunks,times=round(length(masterframe$VBD)/chunks)))
  masterframe<-masterframe %>% filter(masterframe$life_expectancy>=50)
  
  negsample<-masterframe %>% filter(masterframe$cancer==0)
  save(negsample,file = paste("Risksample/negsample.Rdata",sep=""))
  masterframe<-masterframe %>% filter(masterframe$cancer==1)
  
  risksplit<-split(masterframe,masterframe$split)
  
  #Clean up redundant inputs
  rm(masterframe,risksample,PSA_all_p,risk_mat)
  
  #Save risk sample in chunks
    for(i in 1:chunks){
      cancer_col <- paste("X",i,".cancer",sep="") %>% as.name()
      splitsample <- risksplit[i] %>% as.data.frame() %>% filter(!!cancer_col==1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
      save(splitsample,file = paste("Risksample/possample_",i,".Rdata",sep=""))
    }
  } else {
    risksample$split<-(rep(1:chunks,times=round(length(risksample$VBD)/chunks)))
    risksample<-risksample %>% filter(risksample$life_expectancy>=50)
    
      negsample<-risksample %>% filter(risksample$cancer==0)
      save(negsample,file = paste("Risksample/negsample.Rdata",sep=""))
      risksample<-risksample %>% filter(risksample$cancer==1)}
    risksplit<-split(risksample,risksample$split)
    
    rm(risksample,risk_mat)
    
    #Save risk sample in chunks
      for(i in 1:chunks){
        cancer_col <- paste("X",i,".cancer",sep="") %>% as.name()
        splitsample <- risksplit[i] %>% as.data.frame() %>% filter(!!cancer_col==1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
        save(splitsample,file = paste("Risksample/possample_",i,".Rdata",sep=""))
    }
  
}
cmp_create_sample<-cmpfun(create_sample)

# Version of create_sample that stores separate "true" and predicted ten year
# risk:
create_sample_with_misclass<-function(PSA=0,intervals=0,seed=1,screen_strategy){
  
  #Import synthetic dataset derived from PROCAS2 study
  risk_mat<-read.csv("Data/synthetic_risk_data_with_misclassification.csv")[,2:7] %>%
    dplyr::select(!syn.Age)
  colnames(risk_mat)<-c("VBD",
                        "tenyrrisk_est",
                        "liferisk_est",
                        "tenyrrisk_true",
                        "liferisk_true")
  
  #Creat column for risk group to be entered
  risk_mat$risk_group<-numeric(length(risk_mat$VBD))
  
  #Set VDG based on breast density
  risk_mat$VDG<-1+findInterval(risk_mat[,1],VDG_interval)
  
  #Breast density cut-offs for supplemental sreening
  density_cutoff <-3
  
  #Set up data frame of women's lifetimes to simulate
  risksample<-risk_mat[sample(nrow(risk_mat),inum*ifelse(PSA==1,mcruns,1),replace=TRUE),]
  risksample[,c("MRI_screen","US_screen","risk_predicted","feedback",
                "interval_change","life_expectancy","cancer",
                "clinical_detect_size",
                "growth_rate")]<-numeric(length=length(risksample$VBD))
  
  #If risk-stratified screening used then determine if each woman chooses to have
  #risk predicted, attends risk consultation, and changes interval
  
    risksample$risk_predicted<-rbinom(length(risksample$VBD),1,risk_uptake)
    risksample$feedback<-ifelse(risksample$risk_predicted==1 & 
                                  rbinom(length(risksample$VBD),1,(risk_feedback))==1,1,0)
    risksample$interval_change<-ifelse(risksample$feedback==1 & 
                                         rbinom(length(risksample$VBD),1,risk_feedback)==1,1,0)
  
  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy<- rweibull(n = length(risksample$life_expectancy),
                                        shape = acmmortality_wb_a, 
                                        scale = acmmortality_wb_b)
  risksample$life_expectancy<-ifelse(risksample$life_expectancy<=start_age,
                                     qweibull(p = dqrunif(n = 1,
                                                          min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),
                                                          max = 1),
                                              shape = acmmortality_wb_a, scale = acmmortality_wb_b),risksample$life_expectancy)
  risksample$life_expectancy<-ifelse(risksample$life_expectancy>=rep(100,length=length(risksample$life_expectancy)),
                                     100,
                                     risksample$life_expectancy)
  #Determine if a cancer will develop
  risksample$cancer<-ifelse(dqrunif(length(risksample$cancer),0,1)
                            <(risksample$liferisk_true/100),1,0)
  
  #Set clinical detection size for cancer
  risksample$clinical_detect_size<-risksample$cancer*(dqrnorm(
    n = length(risksample$cancer),
    mean = clin_detection_m,sd = clin_detection_sd))
  
  #Set growth rate
  risksample$growth_rate<-risksample$cancer*qlnorm(
    dqrunif(length(risksample$cancer),0,1),
    meanlog=log_norm_mean,sdlog=sqrt(log_norm_sd))
  
  if(PSA==1){
    if(intervals==0){
      #########################Add Monte Carlo Draws into Sample##############  
      
      #Draw stage I to III survival parameters
      survmvn<-data.frame(c(-5.46208,-5.2077,-5.8016),
                          c(-3.8163,-3.75901,-3.8811),
                          c(-2.72264,-2.66053,-2.78617))
      survcovmat<-cov(survmvn)
      survmeans<-c(survmvn[1,1],survmvn[1,2],survmvn[1,3])
      PSA_gamma_survival<-mvrnorm(mcruns,survmeans,survcovmat) %>%
        data.frame() %>%
        transpose()
      
      # Draw Metatstatic survival parameters
      metmvn<-data.frame(c(-1.78723,-1.67922,-1.89434),c(-1.38762,-1.33512,-1.49956),c(-1.01051,-0.93338,-1.08304))
      metmat<-cov(metmvn)
      metmeans<-c(metmvn[1,1],metmvn[1,2],metmvn[1,3])
      PSA_meta_survival<-mvrnorm(mcruns,metmeans,metmat) %>%
        data.frame() %>%
        transpose()
      
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
      
      #Drug parameters
      
      # First bring in log hazard ratios from networked analysis
      loghaz_ests <- readRDS("PreventionOutputs.RDS")
      efficacy_ests <- loghaz_ests[1]
      dropout_ests <- loghaz_ests[4]
      
      # Extract parameters for multivariate normal draws
      efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
      efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
      dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
      dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()
      
      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame() %>%
        transpose()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame() %>%
        transpose()
      
      PSA_uptake_1 <- rnorm(mcruns, .71, .1)
      PSA_uptake_2 <- rnorm(mcruns, .71, .1)
      
      #Draw costs
      PSA_cost_strat<-(rlnorm(mcruns,2.13387381,0.06349671)*1.0272)
      PSA_costvar<-rnorm(mcruns,0,0.1020408)
      PSA_costscreen<-rnorm(mcruns,0,0.1020408)
      PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
      PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
      PSA_cost_US<-rnorm(mcruns,0,0.1020408)
      PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
      PSA_cost_drug <- rnorm(mcruns,0,0.1020408)
      
      #Generate utility draws
      utilmat<-data.frame(c(1-0.82,1-0.81,1-0.83),c(1-0.75,1-0.73,1-0.77))
      lnutilmat<-log(utilmat)
      covutil<-cov(lnutilmat)
      utilmeans<-c(log(1-0.82),log(1-0.75))
      PSA_util<-(1-exp(mvrnorm(mcruns,utilmeans,covutil))) %>%
        data.frame() %>%
        transpose()
      
    }else{
      
      PSA_gamma_survival<-data.frame(runif(mcruns,-6,-3),runif(mcruns,-5,-1),runif(mcruns,-3,-0.01))
      PSA_meta_survival<-data.frame(runif(mcruns,-2.5,-0.5),runif(mcruns,-2,-0.3),runif(mcruns,-1.5,-0.1))
      
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
      PSA_log_norm_mean <- runif(mcruns,0.8,1.2)
      PSA_log_norm_sd <- rnorm(mcruns,1.31,0.11)
      
      #Drug parameters
      
      # First bring in log hazard ratios from networked analysis
      loghaz_ests <- readRDS("PreventionOutputs.RDS")
      efficacy_ests <- loghaz_ests[1]
      dropout_ests <- loghaz_ests[4]
      
      # Extract parameters for multivariate normal draws
      efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
      efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
      dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
      dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()
      
      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame() %>%
        transpose()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame() %>%
        transpose()
      
      PSA_uptake_1 <- rnorm(mcruns, .71, .1)
      PSA_uptake_2 <- rnorm(mcruns, .71, .1)
      
      #Costs
      PSA_cost_strat<-rlnorm(mcruns,2.13387381,0.06349671)
      PSA_costvar<-rnorm(mcruns,0,0.1020408)
      PSA_costscreen<-rnorm(mcruns,0,0.1020408)
      PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
      PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
      PSA_cost_US<-rnorm(mcruns,0,0.1020408)
      PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
      PSA_cost_drug <- rnorm(mcruns,0,0.1020408)
      
      PSA_util<-data.frame(runif(mcruns,0.6,0.9),runif(mcruns,0.5,0.8))
    }
    
    #Generate id for monte carlo set
    mcid<-c(1:mcruns)
    
    #Bind monte carlo draws
    PSA_all_p<-cbind(PSA_gamma_survival,
                     PSA_meta_survival,
                     PSA_beta1,
                     PSA_beta2,
                   PSA_Sen_VDG,
                   PSA_MRI_cdr,
                   PSA_US_cdr,
                   PSA_log_norm_mean,
                   PSA_log_norm_sd,
                   PSA_eff,
                   PSA_dropout,
                   PSA_uptake_1,
                   PSA_uptake_2,
                   PSA_cost_strat,
                   PSA_costvar,
                   PSA_util,
                   PSA_costscreen,
                   PSA_cost_follow_up,
                   PSA_cost_biop,
                   PSA_cost_US,
                   PSA_cost_MRI,
                   PSA_cost_drug,
                   mcid)
  PSA_all_p<-as.data.frame(PSA_all_p)
  colnames(PSA_all_p)<-c("PSA_gamma_survival_1",
                         "PSA_gamma_survival_2",
                         "PSA_gamma_survival_3",
                         "PSA_meta_survival_54",
                         "PSA_meta_survival_74",
                         "PSA_meta_survival_99",
                         "PSA_beta_1",
                         "PSA_beta_2",
                         'PSA_VDG1_sen',
                         'PSA_VDG2_sen',
                         'PSA_VDG3_sen',
                         'PSA_VDG4_sen',
                         "PSA_MRI_cdr",
                         "PSA_US_cdr",
                         "PSA_log_norm_mean",
                         "PSA_log_norm_sd",
                         "PSA_eff_ana",
                         "PSA_eff_tam",
                         "PSA_dropout_ana",
                         "PSA_dropout_tam",
                         "PSA_uptake_1",
                         "PSA_uptake_2",
                         "PSA_cost_strat",
                         "PSA_costvar",
                         "PSA_util_1to3",
                         "PSA_util_4",
                         "PSA_costscreen",
                         "PSA_cost_follow_up",
                         "PSA_cost_biop",
                         "PSA_cost_US",
                         "PSA_cost_MRI",
                         "PSA_cost_drug",
                         "mcid")
    
    #Bind individual level parameters and monte carlo draws
    masterframe<-data.frame(matrix(nrow=inum*mcruns,ncol=length(risksample[1,])+length(PSA_all_p[1,])))
    masterframe[,1:16]<-risksample
    masterframe[,17:49]<-PSA_all_p
    colnames(masterframe)[1:16]<-colnames(risksample)
    colnames(masterframe)[17:49]<-colnames(PSA_all_p)
    
    #Split the dataframe into chunks for easier computation
    masterframe$split<-(rep(1:chunks,times=round(length(masterframe$VBD)/chunks)))
    masterframe<-masterframe %>% filter(masterframe$life_expectancy>=50)
    
      negsample<-masterframe %>% filter(masterframe$cancer==0)
      save(negsample,file = paste("Risksamplewithmisclass/negsample.Rdata",sep=""))
    risksplit<-split(masterframe,masterframe$split)
    
    #Clean up redundant inputs
    rm(masterframe,risksample,PSA_all_p,risk_mat)
    
    #Save risk sample in chunks
      for(i in 1:chunks){
        cancer_col <- paste("X",i,".cancer",sep="") %>% as.name()
        splitsample <- risksplit[i] %>% as.data.frame() %>% filter(!!cancer_col==1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
        save(splitsample,file = paste("Risksamplewithmisclass/possample_",i,".Rdata",sep=""))
      }
  } else {
    risksample$split<-(rep(1:chunks,times=round(length(risksample$VBD)/chunks)))
    risksample<-risksample %>% filter(risksample$life_expectancy>=50)
    
      negsample<-risksample %>% filter(risksample$cancer==0)
      save(negsample,file = paste("Risksamplewithmisclass/negsample.Rdata",sep=""))
    risksplit<-split(risksample,risksample$split)
    
    rm(risksample,risk_mat)
    
    #Save risk sample in chunks
      for(i in 1:chunks){
        cancer_col <- paste("X",i,".cancer",sep="") %>% as.name()
        splitsample <- risksplit[i] %>% as.data.frame() %>% filter(!!cancer_col==1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
        save(splitsample,file = paste("Risksamplewithmisclass/possample_",i,".Rdata",sep=""))
      }
  }
}

# Work out new version of IncidenceMortality adjusted for effect of drug on
# hazard ratios for a single individual.
get_drug_adj_IM <- function(ind_from_risksample,
                          uptake,
                          persistence,
                          risk_red){
  
  drug_IM <- Incidence_Mortality
  # If takes drug assign time taking, otherwise set to zero
  if (dqrunif(1,0,1) < uptake[ind_from_risksample$risk_group, ind_from_risksample$starting_menses_status]){
    time_taking_drug <- min(rexp(1, rate = persistence[ind_from_risksample$risk_group, ind_from_risksample$starting_menses_status]),
                                                                   course_length)
    hazard_ratio <- (1 - risk_red[ind_from_risksample$risk_group, ind_from_risksample$starting_menses_status]) *
      time_taking_drug * 
      log(1/completion_prob[ind_from_risksample$starting_menses_status])
    new_weibull_scale <- inc_scale / hazard_ratio
    
    
    drug_IM$Cond.on.getting.BC..prob.of.getting.cancer.at.age.t <- dweibull(Incidence_Mortality$age,
                                                                            shape=inc_shape,
                                                                            scale=new_weibull_scale)
  }
  else{
    time_taking_drug <- 0
  }
  
  return(list(time_taking_drug, drug_IM))
}

# The following function is used during PSA to get new risk reduction and
# dropout rate matrices for the preventative drug courses based on log hazard
# ratios drawn during the sample generation.
redraw_drug_pars <- function(risksample){
  
  # New efficacies
  ana_eff <- exp(risksample$PSA_eff_ana)
  tam_eff <- exp(risksample$PSA_eff_tam)
  
  full_course_len <- 5
  
  # Assume constant drop out rate with 77% of individuals reaching 5yr mark
  tam_completion_prob <- .77
  
  # New dropout estimates
  ana_dropout_rate <- risksample$PSA_dropout_ana
  tam_dropout_rate <- risksample$PSA_dropout_tam
  
  # Estimate Anastrozole completion probability from Tamoxifen estimate and hazard ratios:
  ana_completion_prob <- exp(ana_dropout_rate - tam_dropout_rate) * tam_completion_prob
  
  completion_prob <- c(ana_completion_prob, tam_completion_prob)
  
  # Now estimate per-unit-time dropout rates based on exponential time to dropout:
  tam_dropout_rate <- (1. / full_course_len) * log(1 / (tam_completion_prob))
  ana_dropout_rate <- (1. / full_course_len) * log(1 / (ana_completion_prob))
  
  # Work out mean time taking each drug
  mean_tam_length <- full_course_len * (1 - tam_completion_prob) / log(1 / tam_completion_prob)
  mean_ana_length <- full_course_len * (1 - ana_completion_prob) / log(1 / ana_completion_prob)
  
  # Quick check: quantities below give efficacy of taking full five-year course
  # assuming linear relationship between time taking and hazard ratio. If both of
  # these are less than one then the linear model is safe to use because no one
  # goes past this point.
  tam_full_course_eff <- 1 - (1 - tam_eff) * log(1 / tam_completion_prob) / (1 - tam_completion_prob)
  ana_full_course_eff <- 1 - (1 - ana_eff) * log(1 / ana_completion_prob) / (1 - ana_completion_prob)
  if ((tam_full_course_eff>1)|(ana_full_course_eff>1)){
    print("Assumption of linear change to hazard ratio with time taking drug will
        not work for one or both drugs being simulated.")
  }
  
  course_length <- c(5., 5.)
  
  #Assign women to risk groups based on 10yr risk if using risk-stratified approach  
  if(screen_strategy==1 | screen_strategy==9) {
    risk_red <- matrix(c(ana_eff, tam_eff,
                         ana_eff, tam_eff,
                         ana_eff, tam_eff,
                         ana_eff, tam_eff,
                         ana_eff, tam_eff),
                       nrow = 5,
                       ncol = 2)
    
    course_length <- c(5., 5.)
    
    uptake <-rbind(c(0., 0.),
                   c(0., 0.),
                   c(0., 0.),
                   c(.71, .71),
                   c(.71, .71))
    
    persistence <- matrix(c(ana_dropout_rate, tam_dropout_rate,
                            ana_dropout_rate, tam_dropout_rate,
                            ana_dropout_rate, tam_dropout_rate,
                            ana_dropout_rate, tam_dropout_rate,
                            ana_dropout_rate, tam_dropout_rate),
                          nrow = 5,
                          ncol = 2)
  } else
    if(screen_strategy==2) {
      risk_red <- matrix(c(ana_eff, tam_eff,
                           ana_eff, tam_eff,
                           ana_eff, tam_eff),
                         nrow = 3,
                         ncol = 2)
      
      uptake <-rbind(c(0., 0.),
                     c(0., 0.),
                     c(.71, .71))
      
      persistence <- matrix(c(ana_dropout_rate, tam_dropout_rate,
                              ana_dropout_rate, tam_dropout_rate,
                              ana_dropout_rate, tam_dropout_rate),
                            nrow = 3,
                            ncol = 2)
    } else
      if(screen_strategy==7 | screen_strategy==8) {
        risk_red <- matrix(c(ana_eff, tam_eff,
                             ana_eff, tam_eff),
                           nrow = 2,
                           ncol = 2)
        
        uptake <-rbind(c(0., 0.),
                       c(.71, .71))
        
        persistence <- matrix(c(ana_dropout_rate, tam_dropout_rate,
                                ana_dropout_rate, tam_dropout_rate),
                              nrow = 2,
                              ncol = 2)
      }  else{
        risk_red <- matrix(c(ana_eff, tam_eff),
                           nrow = 1,
                           ncol = 2)
        
        uptake <-matrix(c(0., 0.),
                        nrow = 1,
                        ncol = 2)
        
        persistence <- matrix(c(ana_dropout_rate, tam_dropout_rate),
                              nrow = 1,
                              ncol = 2)
      }
  return(list(risk_red, uptake, persistence))
}