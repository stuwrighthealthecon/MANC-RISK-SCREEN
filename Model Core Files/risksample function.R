create_sample<-function(PSA=0,intervals=0,seed=1){

  #Import synthetic dataset derived from PROCAS2 study
  risk_mat<-read.csv("synthetic_risk_data.csv")[,2:4]
  colnames(risk_mat)<-c("VBD","tenyrrisk","liferisk")
  
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
  
  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy<- rweibull(n = length(risksample$life_expectancy),
                                        shape = acmmortality_wb_a, 
                                        scale = acmmortality_wb_b)
  for (i in 1:length(risksample$life_expectancy)) {
    if(risksample[i,"life_expectancy"] <= start_age){risksample[i,"life_expectancy"]<-qweibull(
      p = dqrunif(n = 1,
                  min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),
                  max = 1),
      shape = acmmortality_wb_a, scale = acmmortality_wb_b)}}
  
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
    
  #Costs
  PSA_cost_strat<-rlnorm(mcruns,2.13387381,0.06349671)
  PSA_costvar<-rnorm(mcruns,0,0.1020408)
  PSA_costscreen<-rnorm(mcruns,0,0.1020408)
  PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
  PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
  PSA_cost_US<-rnorm(mcruns,0,0.1020408)
  PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
    
  PSA_util<-data.frame(runif(mcruns,0.6,0.9),runif(mcruns,0.5,0.8))
  }
  
  #Generate id for monte carlo set
  mcid<-c(1:mcruns)
  
  #Bind monte carlo draws
  PSA_all_p<-cbind(PSA_gamma_survival,PSA_meta_survival,PSA_beta1,PSA_beta2,
                   PSA_Sen_VDG,PSA_MRI_cdr,PSA_US_cdr,PSA_log_norm_mean,
                   PSA_log_norm_sd,PSA_cost_strat,PSA_costvar,PSA_util,PSA_costscreen,
                   PSA_cost_follow_up,PSA_cost_biop,PSA_cost_US,PSA_cost_MRI,mcid)
  PSA_all_p<-as.data.frame(PSA_all_p)
  colnames(PSA_all_p)<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                         "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                         "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                         'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                         "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                         "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                         "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","mcid")
  
  #Bind individual level parameters and monte carlo draws
  masterframe<-data.frame(matrix(nrow=inum*mcruns,ncol=length(risksample[1,])+length(PSA_all_p[1,])))
  masterframe[,1:14]<-risksample
  masterframe[,15:40]<-PSA_all_p
  colnames(masterframe)[1:14]<-colnames(risksample)
  colnames(masterframe)[15:40]<-colnames(PSA_all_p)
  
  #Split the dataframe into chunks for easier computation
  masterframe$split<-(rep(1:chunks,times=round(length(masterframe$VBD)/chunks)))
  masterframe<-masterframe %>% filter(masterframe$life_expectancy>=50)
  risksplit<-split(masterframe,masterframe$split)
  
  #Clean up redundant inputs
  rm(masterframe,risksample,PSA_all_p,risk_mat)
  
  #Save risk sample in chunks
  for(i in 1:chunks){
    splitsample<-as.data.frame(risksplit[i])
    save(splitsample,file = paste("Risksample/risksample_",i,".Rdata",sep=""))
  }
  } else {
    risksample$split<-(rep(1:chunks,times=round(length(risksample$VBD)/chunks)))
    risksample<-risksample %>% filter(risksample$life_expectancy>=50)
    risksplit<-split(risksample,risksample$split)
    
    rm(risksample,risk_mat)
    
    #Save risk sample in chunks
    for(i in 1:chunks){
      splitsample<-as.data.frame(risksplit[i])
      save(splitsample,file = paste("Risksample/risksample_",i,".Rdata",sep=""))
    }
    
  }
}

# Version of create_sample that stores separate "true" and predicted ten year
# risk:
create_sample_with_misclass<-function(PSA=0,intervals=0,seed=1){
  
  #Import synthetic dataset derived from PROCAS2 study
  risk_mat<-read.csv("synthetic_risk_data_with_misclassification.csv")[,2:7] %>%
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
  
  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy<- rweibull(n = length(risksample$life_expectancy),
                                        shape = acmmortality_wb_a, 
                                        scale = acmmortality_wb_b)
  for (i in 1:length(risksample$life_expectancy)) {
    if(risksample[i,"life_expectancy"] <= start_age){risksample[i,"life_expectancy"]<-qweibull(
      p = dqrunif(n = 1,
                  min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),
                  max = 1),
      shape = acmmortality_wb_a, scale = acmmortality_wb_b)}}
  
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
      
      #Costs
      PSA_cost_strat<-rlnorm(mcruns,2.13387381,0.06349671)
      PSA_costvar<-rnorm(mcruns,0,0.1020408)
      PSA_costscreen<-rnorm(mcruns,0,0.1020408)
      PSA_cost_follow_up<-rnorm(mcruns,0,0.1020408)
      PSA_cost_biop<-rnorm(mcruns,0,0.1020408)
      PSA_cost_US<-rnorm(mcruns,0,0.1020408)
      PSA_cost_MRI<-rnorm(mcruns,0,0.1020408)
      
      PSA_util<-data.frame(runif(mcruns,0.6,0.9),runif(mcruns,0.5,0.8))
    }
    
    #Generate id for monte carlo set
    mcid<-c(1:mcruns)
    
    #Bind monte carlo draws
    PSA_all_p<-cbind(PSA_gamma_survival,PSA_meta_survival,PSA_beta1,PSA_beta2,
                     PSA_Sen_VDG,PSA_MRI_cdr,PSA_US_cdr,PSA_log_norm_mean,
                     PSA_log_norm_sd,PSA_cost_strat,PSA_costvar,PSA_util,PSA_costscreen,
                     PSA_cost_follow_up,PSA_cost_biop,PSA_cost_US,PSA_cost_MRI,mcid)
    PSA_all_p<-as.data.frame(PSA_all_p)
    colnames(PSA_all_p)<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                           "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                           "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                           'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                           "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                           "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                           "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","mcid")
    
    #Bind individual level parameters and monte carlo draws
    masterframe<-data.frame(matrix(nrow=inum*mcruns,ncol=length(risksample[1,])+length(PSA_all_p[1,])))
    masterframe[,1:14]<-risksample
    masterframe[,15:40]<-PSA_all_p
    colnames(masterframe)[1:14]<-colnames(risksample)
    colnames(masterframe)[15:40]<-colnames(PSA_all_p)
    
    #Split the dataframe into chunks for easier computation
    masterframe$split<-(rep(1:chunks,times=round(length(masterframe$VBD)/chunks)))
    masterframe<-masterframe %>% filter(masterframe$life_expectancy>=50)
    risksplit<-split(masterframe,masterframe$split)
    
    #Clean up redundant inputs
    rm(masterframe,risksample,PSA_all_p,risk_mat)
    
    #Save risk sample in chunks
    for(i in 1:chunks){
      splitsample<-as.data.frame(risksplit[i])
      save(splitsample,file = paste("Risksample/risksample_",i,".Rdata",sep=""))
    }
  } else {
    risksample$split<-(rep(1:chunks,times=round(length(risksample$VBD)/chunks)))
    risksample<-risksample %>% filter(risksample$life_expectancy>=50)
    risksplit<-split(risksample,risksample$split)
    
    rm(risksample,risk_mat)
    
    #Save risk sample in chunks
    for(i in 1:chunks){
      splitsample<-as.data.frame(risksplit[i])
      save(splitsample,file = paste("Risksamplewithmisclass/risksample_",i,".Rdata",sep=""))
    }
    
  }
}

# Work out new version of IncidenceMortality adjusted for effect of drug on hazard ratios.
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