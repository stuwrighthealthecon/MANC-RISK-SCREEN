#Lookup function for costs
fnLookupBase <- function(iStage, iAge, iLE) {
  as.numeric(tblLookup$CDCost.p.i.d[tblLookup$Stage==iStage & tblLookup$Age==iAge & tblLookup$Yr==iLE])
}

#Load functions required for model
Incidence_function <- function(){
  
  #x<- cumsum(Incidence_Mortality[,2][21:101])
  #Sample an incidence time (based on vector of probabilities of getting cancer at age t conditional on getting cancer and surviving to age t)
  incidence_time_1 <- sample(x = Incidence_Mortality[,1][start_age:101],size = 1,prob = Incidence_Mortality[,2][start_age:101])
  #add within year time (as uniform)
  incidence_time <- incidence_time_1+ dqrunif(1,0,1)
  #back calculate the genesis time of the tumour
  
  #First determine if screen detected or clinical detected in current data 
  
  detect_mode <- 1 #clinical detected
  if (incidence_time <= screen_endage & incidence_time >= screen_startage & dqrunif(1,0,1)<prop_screen_detected){detect_mode <- 0} #screen detected
  
  #size at detection - as number of tumour doublings in diameter from a 0.25mm diameter
  clin_detect_size_g <- risk_data[13]
  clin_detect_size_g <- start_size*2^clin_detect_size_g
  ca_size_incidence <- clin_detect_size_g
  if(detect_mode == 0){
    screen_detect_size_g <- dqrnorm(n = 1,mean = screen_detection_m,sd = screen_detection_sd)
    screen_detect_size_g[screen_detect_size_g < 3.5] <- 3.5 # to prevent unrealistic left tail
    screen_detect_size_g[screen_detect_size_g >= 9] <- 8.99 # to prevent unrealistic right tail
    screen_doubles <- screen_detect_size_g #size as number of doublings
    screen_detect_size_g <- start_size*2^screen_detect_size_g 
    ca_size_incidence <- screen_detect_size_g
    #generate a valid potential clinical detection time - ie after screen detection
    if(clin_detect_size_g< screen_detect_size_g){
      
      clin_detect_size_g <- qnorm(dqrunif(n = 1,min = pnorm(screen_doubles,clin_detection_m,clin_detection_sd),max = 1),mean = clin_detection_m,sd = clin_detection_sd)
      clin_detect_size_g[clin_detect_size_g >= 9] <- 8.999 # to prevent unrealistic right tail
      clin_detect_size_g <- start_size*2^clin_detect_size_g
    }
  }
  
#Incidence Simulator
#Load necessary ONS data

  #generate (non-cancer) mortality time conditional on survivng to time t(peroid covered by BC survival function)
  if (incidence_time<90){
    mort_time <- sample(x = Incidence_Mortality[,1][incidence_time_1:101],size = 1,prob = Incidence_Mortality[,3][incidence_time_1:101])
  }else{mort_time <- 101} # if get cancer at >90 then all time in simulation is covered 
  # to replace mort_age
  #How does this interact with the calculation of overall survival given breast cancer diagnosis?
  #Breast Cancer Survial reverts to mort_age after 10 years
  
  result<-c(incidence_time,detect_mode,ca_size_incidence, clin_detect_size_g,mort_time)
  return(result)
}
cmp_incidence_function<-cmpfun(Incidence_function)

#stage calculator 

#matrix with proporiton in each ?stage group for each size category
#Updated 1010 include DCIS
stage_by_size <- function(Ca_size){
  stage_cat <- 0
  
  #first determine if advanced cancer or not based on metastatic prob by size (categorical)
  if(Ca_size<= 25){m_size <- 25}else{m_size <- ceiling((Ca_size-25)/10)*10+25}
  if (m_size > 85){m_size <- 85}
  if(dqrunif(1,0,1) < metastatic_prob[match(m_size, metastatic_prob[,1]),2] && stage_cat == 0){stage_cat <- 4} 
  
  #sample from categories 1,2 & 3 with probability of each based on the correct row of stagebysize matrix 
  # Ca_size is continuous, need to match to closest larger value in stage_by_size column 1
  if(stage_cat == 0){
    size_cat <- findInterval(Ca_size,ca_size_cut)
    stage_cat <- sample(x=c(1,2,3,5),size = 1,prob = c(stage_by_size_mat[size_cat,])) #1 best 3 worst prognosis
  }
  #return the stage category
  result <- stage_cat
  
  return(result)
}
cmp_stage_by_size<-cmpfun(stage_by_size)

#Screen test results simulation
#Updated 0915 with weedon-fekjaer 2008 estimates
#args (inputs) are tumour diameter, VDG, MRI_screening(0/1), US_screening(0/1) 
screening_result <- function(Ca_size,VDG,MRI_screening,US_screening){
  
  #Caculate size/density specific sensitivity  
  Sensitivity <- if(
    exp((Ca_size - beta2)/beta1)/(1+exp((Ca_size-beta2)/beta1))>sensitivity_max){sensitivity_max}
  else{exp((Ca_size - beta2)/beta1)/(1+exp((Ca_size-beta2)/beta1))} #use to set max sensitivity 0.95
  
  dense_OR <- (Sen_VDG[VDG]/(1-Sen_VDG[VDG]))/(Sen_VDG_av/(1-Sen_VDG_av))
  Sensitivity <- ((Sensitivity/(1-Sensitivity))*dense_OR)/(1+((Sensitivity/(1-Sensitivity))*dense_OR))
  
  rnd_1 <- dqrunif(1,0,1) # random number used to compare to Sensitivity with and without supplemental screening
  if(rnd_1<Sensitivity){
    Screen_detected_ca <- 1
    Mammo_detected_ca <- 1 #keep track of which stage it is detected
  }else{Screen_detected_ca <-0
  Mammo_detected_ca <- 0} #test if cancer is detected by screening mammo
  #Additional screening depends on density/risk
  if(Screen_detected_ca == 0){
    if(MRI_screening == 1){
      MRI_supp_odds <- (Sensitivity/(1-Sensitivity))*((MRI_cdr+Mammo_cdr)/Mammo_cdr)
      MRI_supp_sens <- MRI_supp_odds/(MRI_supp_odds+1)
      if(rnd_1 < MRI_supp_sens){
        Screen_detected_ca <- 1
        MRI_detected_ca <- 1
      }else{MRI_detected_ca <- 0}}else{MRI_detected_ca <- 0}
    
    if(US_screening==1){
      US_supp_odds <- (Sensitivity/(1-Sensitivity))*((US_cdr+Mammo_cdr)/Mammo_cdr)
      US_supp_sens <- US_supp_odds/(US_supp_odds+1)
      if(rnd_1 < US_supp_sens){
        Screen_detected_ca <-1
        US_detected_ca <- 1
      }else{US_detected_ca <- 0}}else{US_detected_ca <- 0}}else{US_detected_ca <- 0
      MRI_detected_ca <- 0}
  
  #Uses estimate based on increased cancer detection rate - check this then matches the number of interval cancers that are now found at screening
  
  result <- c(Screen_detected_ca,Mammo_detected_ca,MRI_detected_ca,US_detected_ca)
  return(result)
}
cmp_screening_result<-cmpfun(screening_result)

#Simulate survival by stage 
#Needs to know stage_cat to generate a survival time from current age (currently 10-yr with cancer survival is irrespective of age), mort_age and age are need for those suviving beyond 10 years.

#Further assumption to guard against (reverse)lead-time bias is that cancer-specific survival is calculated from the age the cancer would have been clinically detected. Assumes no mortality effect of treatment.

Ca_survival_time <- function(stage_cat, Mort_age,age, CD_age){
  
  if (stage_cat< 4){
    survival_time <- -(log(x = dqrunif(1,0,1))/gamma_stage[stage_cat]) #inverse of cdf when rate is gamma_stage[x]
    
    #adjust for additional mortality at ages above 65
    if (CD_age > 65){
      survival_time <- -(log(x = dqrunif(1,0,1))/((Incidence_Mortality[min((floor(CD_age)+1),100),5]/Incidence_Mortality[66,5])*gamma_stage[stage_cat]))
    }
    #data are for 10-year survival, after 10 years assume that pop mortality rates apply
    if(survival_time > 10){
      Mort_age <- qweibull(p = dqrunif(n = 1,min = pweibull(q = CD_age+10,shape = acmmortality_wb_a,scale = acmmortality_wb_b),max = 1),shape = acmmortality_wb_a, scale = acmmortality_wb_b)
      if(Mort_age > time_horizon){Mort_age <- time_horizon}
      survival_time <- Mort_age - age
    }
  }
  
  if (stage_cat == 4){
    if (age < 55){age_cat_M <- 1}
    if (age >=55 && age <75){age_cat_M <- 2}
    if (age >= 75){age_cat_M <- 3}
    survival_time <- -(log(dqrunif(1,0,1))/metastatic_survival[age_cat_M])
    #check lifetime does not exceed horizon and set to less than 100 if it does
    if (CD_age+survival_time >=100){survival_time <- time_horizon - CD_age}
  }
  if (stage_cat == 5){
    survival_time <- (Mort_age-CD_age) # no effect on mortality
  }
  
  if(CD_age+survival_time > time_horizon){survival_time <- time_horizon-CD_age}
  result <- CD_age+survival_time
  return(result)
}
cmp_ca_survival_time<-cmpfun(Ca_survival_time)