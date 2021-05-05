## ASSURE Breast Cancer Screening DES model##
#v2
#Gray et al, 2017, Evaluation of a stratified national breast screening program in the United Kingdom: an early model-based cost-effectiveness analysis, Value in Health 20 (8), 1100-1109
#Adaptations/corrections from version 1 used in the previous publication. To be documented in future publication.
#May be adapted for non-commercial use. Please cite the above study or more recent publications as the source.

#NOTES#
#necessary libraries - see below
#functions used within the model were stored as seperate R files - this is only for reading convenience
#Model uses individual patient sampling that can be acheived in different ways using the i, j and ii loops. 
#Assumes using parallel computing with doParallel package - computation time will still be substantial on a single computer
#One alternative is simulated at a time - controlled by the "Programme parameters" section inputs.
#Output is costs, QALYs, etc for that alternative as a dataframe and stored as a csv. These need to be further combined in another script or spreadsheet with results from other alternatives for cost-effectiveness analysis.

#Call libraries
library("doParallel")
library("MASS")
#set to directory containing other necessary files
setwd(dir="C:/Users/mdxassw4/Dropbox (The University of Manchester)/MERCADO/ASSURE model Export/ASSURE model Export")

#Creating clusters for parallel computations (available cores)
registerDoParallel(cores = 5)

#Start the clock
ptm <- proc.time()

#Call functions

source(file="incidence_sim_v2.R")

source(file = "screening_result_sim_USMRI2.R")

source(file = "NPI_by_size_sim_DCIS.R")

source(file = "survival_times_NPI_sim_v3.R")


  #Set the numbers in i (outer individual risk factor sampler) and j (inner random tumor etc variable draws) loops
  inum <- 1000000 #10000000 needed for the results
  jnum <- 1
  
  for (ii in 1:10){ #loop for repeating 1 mill simulation and saving
    
  ######################################Programme parameters######################################################
  
  #Supplental Screening on/off 1=on, 0=off
  supplemental_screeing <- 0
  #Risk stratification
  risk_strategy <- 4 #1=procas, 2=optimal, 3=risk tertiles, 4= 3-yearly, 5= 2-yearly, 6= 5-yearly, 7=2 rounds (50 and 60) other=no screening
  
  #Age of an individual at start of simulation
  start_age <- 38
  
  #Time horizon
  time_horizon <- 100 #sets time horizon
  
  #Discount rate
  discount_health <- 0.035 #the discount rate for health benefits
  discount_costs <- 0.035  #the discount rate for costs
  ######################################End Programme parameters###################################################
  
  #############################################################################################################
  
  #Parameters not varied in the PSA#
  #Proportion of cancers in current screening age range that were screen-detected
  prop_screen_detected <- 0.5 #proportion from screening programme reprots UPDATE
  #screen-detected size parameters for calculting tumour incidence/genesis time
  screen_detection_m <- 4.12
  screen_detection_sd <- 3.93
  
  #All-cause mortality
  mortality_wb_a <- 8.97 #weibull shape parameter for all cause mortality
  mortality_wb_b <- 86.74 #weibull scale parameter for all cause mortality
  
  #Post-BC all cause mortality
  gamma_survival_3 <- exp(-2.465) #exponential distribution scale parameter NPI 3
  gamma_survival_2 <- exp(-4.023)
  gamma_survival_1 <- exp(-5.413)
  
  #Incidence distribution
  Incidence_Mortality <- read.csv("Incidence_Mortality_ONS2.csv")
  
  #Risk factor and density variable correlations from PROCAS
  risk_factor_cov <- read.csv(file="risk_factor_covariance.csv",header=TRUE)[,c(2,3,4)]
  rf_means <- c(2.427836, 3.038518, 13.209247)
  
  #metastatic cancer parameters
  metastatic_prob <- read.csv(file = "metastatic_prob.csv",header = FALSE)
  
  #DCIS parameters
  DCIS_fraction <- 0.21 #UK screening programme statistics proporiton of SD cancers that are DCIS
  NPI_by_size_mat <- read.csv(file = "NPIbySize.csv")
  
  #Load age adjsuted utility weights
  utility_ages <- read.csv(file = "age specific utility weights.csv",header = FALSE)

  #Clinical detection parameters
  clin_detection_m <- 6.5 #normal mean parameter for tumour doubling at detection
  clin_detection_sd <- 0.535 #standard deviation for  doubling at detection
  
  #Screen detection size parameters
  screen_detection_m <- 6.12 #normal mean parameter for tumour doubling at detection
  screen_detection_sd <- 0.96 #standard deviation for  doubling at detection
  
  #Mammography with sensitivity conditional on tumour diameter parameters W-F
  beta1 <- 1.47 
  beta2 <- 6.51
  #Mammography sensitivity by VDG PREVENTICON
  Sen_VDG <- c(0.85,0.776,0.69,0.586)
  Sen_VDG_av <- 0.735
  
  #Supplemental screening sensitivity parameters CEPAC
  Mammo_cdr <- 4.2 #Cancer detection rate per 1000 high dense screens Mammo CEPAC
  MRI_cdr <- 5 #CDR for MRI in Mammo negative women (incremental)
  US_cdr <- 3 #CDR for US in Mammo negative women (incremental)
  
  
  #set an empirical distribution for breast cancer incidence
  #Use a discrete distribution, sampled using a uniform dist draw and then a linear approximation by using a draw from a triangular distribution to place within the interval. First two variables give age bands and risks
  age_bands <- seq(15,95,10)
  age_band_risk <- c(0,2,9,13,15,18,20,24,0)
  age_band_prob <- age_band_risk/100
  #Growth rate distribution parameters
  log_norm_mean <- 1.07
  log_norm_sd <- 1.31
  max_size <- 128 #mm diameter
  start_size <- 0.25 #starting size of tumours, diameter in mm
  Vc = (4/3)*pi*(start_size/2)^3 #Volume at start
  Vm = (4/3)*pi*(max_size/2)^3 #Max volume
  
  #Metatstatic survival parameters
  meta_survival_49 <- -0.527 #age <= 49
  meta_survival_69 <- -0.537 #age 50-69
  meta_survival_99 <- -0.849 # 70-99
  
  
  #Screening Times
  screen_startage <- 50
  screen_endage <- 70
  high_risk_screentimes <- seq(screen_startage,screen_endage,1)
  med_risk_screentimes <- seq(screen_startage,screen_endage,2)
  low_risk_screentimes <- seq(screen_startage, screen_endage,3)
  screen_count <- 1 #counter for number of screens so far
  #Screen Sensitivity Parameters
  #sensitivity_lognorm_mean <- 6
  #sensitivity_lognorm_sd <- 1
  sensitivity_max <- 0.95
  #Risk cut-offs
  risk_cutoffs_optimal <- c(1.6,2.4,2.7,4.1,5.2,7.6,15) #optimal cut-offs v1
  risk_cutoffs_procas <- c(2,3.5,5,8,100) #procas plan
  risk_cutoffs_tert <- c(2.328355,3.067665) #tertiles of risk
  #tertile cut-offs
  #Density cut-off
  density_cutoff <- 3 #VDG groups 3 and 4
  
  #######################Cost Data#########################################
  #From Madan, inflated to 2015 prices
  cost_screen <- 54
  
  cost_follow_up <- 95
  
  cost_biop <- 160
  
  cost_DCIS <- 8806
  
  cost_NPI_Good <- 11630
  
  cost_NPI_Moderate <- 12978
  
  cost_NPI_Poor <- 15405
  
  cost_metastatic <- 23449
  
  cost_US <- 80
  
  cost_MRI <- 220
  
  
  ####################False Positive and Overdiagnosis parameter################
  recall_rate <- 0.045 #approx UK recall rate
  biopsy_rate <- 0.024 #proporiton of referrals without cancer that have biopsy - Madan
  
  ##########################Utility Weigts#########################################
  # NPI-group specific utility decrements
  #Lidgren 0.696 (mean age 57, range(28-93)), metastatic 0.685 permanent
  #1st year
  utility_NPI1_y1 <- 0.696/0.822
  utility_NPI2_y1 <- 0.696/0.822
  utility_NPI3_y1 <- 0.696/0.822
  
  #Folowing Years
  #0.779
  utility_NPI1_follow <- 0.779/0.822
  utility_NPI2_follow <- 0.779/0.822
  utility_NPI3_follow <- 0.779/0.822
  
  #Metastatic cancer (permanent decrement)
  utility_metastatic <- 0.685/0.822
  
  #DCIS assumed no effect
  utility_DCIS <- 1
  
  utility_NPI_cat_y1 <- c(utility_NPI1_y1, utility_NPI2_y1,utility_NPI3_y1, utility_DCIS,utility_metastatic) #NPI1, NPI2, NPI3, DCIS, metastatic
  utility_NPI_cat_follow <- c(utility_NPI1_follow,utility_NPI2_follow,utility_NPI3_follow,utility_DCIS,utility_metastatic)
  
  ########################Model Diagnostics##############################################
  #cancer_counter <- 0

  #############Individual patient-level siumulation start################################
  
  
  #Start the clock
  #ptm <- proc.time()
  
  ############################Loop for PSA results######################
  
  ###########Outer Individual sampling loop##############################  
  #results <- c()
    total_screens <- 0
    total_cancers_detected <- 0
    total_costs <- 0
    total_US_costs <- 0
    total_MRI_costs <- 0
    total_life_years <- 0
    total_US <- 0
    total_MRI <- 0
    total_QALYs <- 0
    total_costs_follow_up <- 0
    
    #Start of i loop
    results <- foreach(i=1:inum,.combine = 'rbind',.packages = 'MASS') %dopar% { #parellised foreacg loop
      cancer_diagnostic <- rep(0,10) #stores age, size, mode of detection of each detected cancer
     
      #for (i in 1:inum){  
      risk_data <- mvrnorm(n=1,mu = rf_means,Sigma = risk_factor_cov)
      liferisk <- risk_data[3]/100 #as probability not percentage
      
      #ten_year_risk <- risk_factors[sample_vector[i],6]
      ten_year_risk <- risk_data[2] #As a percentage
      
      #Set risk group
      if(risk_strategy == 1){
        if (ten_year_risk<risk_cutoffs_procas[1]) {risk_group <- 1} 
        if(ten_year_risk>=risk_cutoffs_procas[1] & ten_year_risk <risk_cutoffs_procas[2]) {risk_group <- 2}
        if(ten_year_risk>=risk_cutoffs_procas[2] & ten_year_risk <risk_cutoffs_procas[3]) {risk_group <- 3}
        if (ten_year_risk>=risk_cutoffs_procas[3] & ten_year_risk <risk_cutoffs_procas[4]) {risk_group <- 4} 
        if (ten_year_risk >= risk_cutoffs_procas[4]){risk_group <- 5}
      }
      if(risk_strategy == 2){
        if (ten_year_risk<risk_cutoffs_optimal[1]) {risk_group <- 1} 
        if(ten_year_risk>=risk_cutoffs_optimal[1] & ten_year_risk <risk_cutoffs_optimal[2]) {risk_group <- 2}
        if(ten_year_risk>=risk_cutoffs_optimal[2] & ten_year_risk <risk_cutoffs_optimal[3]) {risk_group <- 3}
        if(ten_year_risk>=risk_cutoffs_optimal[3] & ten_year_risk <risk_cutoffs_optimal[4]) {risk_group <- 4}
        if(ten_year_risk>=risk_cutoffs_optimal[4] & ten_year_risk <risk_cutoffs_optimal[5]) {risk_group <- 5}
        if(ten_year_risk>=risk_cutoffs_optimal[5] & ten_year_risk <risk_cutoffs_optimal[6]) {risk_group <- 6}
        if(ten_year_risk>=risk_cutoffs_optimal[7]) {risk_group <- 7}
      }
      if(risk_strategy ==3){
        if (ten_year_risk<risk_cutoffs_tert[1]) {risk_group <- 1} 
        if(ten_year_risk>=risk_cutoffs_tert[1] & ten_year_risk <risk_cutoffs_tert[2]) {risk_group <- 2}
        if(ten_year_risk>=risk_cutoffs_tert[2]) {risk_group <- 3}
      }
      
      #Set density group
      
      VDG <- round(risk_data[1],0) #temporary, need to check VDG cutpoints and simulate density then grade
      if(VDG >4){VDG <- 4}
      if(VDG<1){VDG <- 1}
      
      #Select the MRI and US screening based on simple rules
      if(supplemental_screeing == 0){ 
        MRI_screening <- 0
        US_screening <- 0
      }else{
        if(VDG>=density_cutoff){
          if(ten_year_risk >= 8){MRI_screening <- 1}else{US_screening <- 1}
        }else{
          MRI_screening <- 0
          US_screening <- 0
        }}
      # Indicators for supplemental screening set
      
      ###############Screen times###############################
      screen_times <- c(999)
      if(risk_strategy ==1){
        if(risk_group ==1 | risk_group ==2){screen_times = low_risk_screentimes}
        if(risk_group ==3 | risk_group ==4){screen_times = med_risk_screentimes}
        if(risk_group ==5){screen_times = high_risk_screentimes}
      }
      if(risk_strategy ==2){
        if(risk_group ==1){screen_times = c(47,51,56,59,63,69)}
        if(risk_group ==2){screen_times = c(47, 51, 56, 59, 62, 64, 67, 69)}
        if(risk_group ==3){screen_times = c(47, 51, 54, 56, 59, 62, 64, 67, 69)}
        if(risk_group ==4){screen_times = c(47, 49, 51, 54, 56, 59, 62, 64, 66, 67, 69)}   
        if(risk_group ==5){screen_times = c(47, 49, 51, 54, 56, 59, 62, 63, 65, 66, 67, 69)}
        if(risk_group ==6){screen_times = c(47, 49, 51, 54, 56, 59, 60, 61, 62, 63, 65, 66, 67, 69)}
        if(risk_group ==7){screen_times = c(47, 48, 49, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, 63, 65, 66, 67, 69)}
      }
      if(risk_strategy == 3){
        if(risk_group ==1){screen_times = low_risk_screentimes}
        if(risk_group ==2){screen_times = med_risk_screentimes}
        if(risk_group ==3){screen_times = high_risk_screentimes}
      }
      if(risk_strategy ==4){ 
        #Changed to vary screening times by individual to produce more realistic incidence data
        start_screening <- 47+runif(n = 1,min = 0,max = 4)
        screen_times = seq(start_screening, start_screening+(3*8),3) #eight invitations offerred
        #if(runif(1,0,1)>0.72){screen_times<- c(999)} #72% uptake of screening
      }
      if(risk_strategy ==5){
        screen_times = med_risk_screentimes
      }
      if(risk_strategy == 6){
        start_screening <- 47+runif(n = 1,min = 0,max = 4)
        screen_times <- seq(start_screening, start_screening+(5*4),5) #assume 5 rounds
      }
      if(risk_strategy == 7){
        start_screening <- 47+runif(n = 1,min = 0,max = 4)
        screen_times <- seq(start_screening, start_screening+10,10)
      }
      ##########Counters i loop level######################
      #screen-detected cancer counts
      screen_detected_count <- 0
      sdfirst_counter <- 0
      sdlast_counter <-0
      #count of screens
      screen_counter <- 0
      lastscreen_counter <-0
      US_counter <- 0
      MRI_counter <- 0
      #recall count
      recall_counter <- 0
      #total cost
      cost_counter <- 0
      #total life years
      LY_counter <- 0
      #total QALYs
      QALY_counter <- 0
      #NPI group cancers counters
      NPI1_counter <- 0
      NPI2_counter <- 0
      NPI3_counter <- 0
      NPI4_counter <- 0
      NPI5_counter <- 0
      ###Preload incidence, mortality and clinical detection times for j cases
      
      mort_sample<- rweibull(n = jnum,shape = mortality_wb_a, scale = mortality_wb_b)
      clin_detect_sample <- rnorm(n = jnum,mean = clin_detection_m,sd = clin_detection_sd)
      clin_detect_sample[clin_detect_sample < 4] <- 4 # to prevent unrealistic left tail
      clin_detect_sample[clin_detect_sample >= 9] <- 8.99 # to prevent unrealistic right tail
      
      #J histories
      for (j in 1:jnum){
        #reset screen count
        #J level counters
        screen_count <- 0
        recall_count <- 0
        sdlast_cancer <-0
        lastscreen_count <- 0
        sdfirst_cancer <- 0
        NPI_cat <- 0
        MRI_count <- 0
        US_count <- 0
        incidence_age_record <- 0
        costs <- 0
        US_costs <- 0
        MRI_costs <- 0
        costs_follow_up <- 0

        #Lifetime cancer incidence#
        #determine if cancer occus, then set age of incidence 
        if (runif(1,0,1)<liferisk){
          ca_case <- 1
        #Individual growth rate
        grow_rate_i <- qlnorm(runif(1,0,1),meanlog = log_norm_mean,sdlog = sqrt(log_norm_sd))
        #Incidence age (under current programme)
        ca_incidence_i <- Incidence_function()
        ca_incidence_age <- ca_incidence_i[1]
        print(ca_incidence_i)
        #Clinical detection age
        CD_size <-  ca_incidence_i[4]#tumour diameter at CD

        if(ca_incidence_i[2] ==1){CD_age <- ca_incidence_i[1]} #endif
        if(ca_incidence_i[2] ==0){
          #Age at clinical detection of a currently screen detected tumour - incidence age + time to grow to clincial detection dize - time to grow to screen detection size
          CD_age <- ca_incidence_i[1] + ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[4]/2)^3))^0.25-1))/(0.25*grow_rate_i)) - 
            ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[3]/2)^3))^0.25-1))/(0.25*grow_rate_i))
          } #endif

        cancer_diagnostic[8] <- c(CD_age)
        t_gen <- ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(CD_size/2)^3))^0.25-1))/(0.25*grow_rate_i)) #time from genesis
        gen_age <- CD_age - t_gen # genesis age
        
        }else{
          ca_case <- 0
          ca_incidence_age <- 999 #redundent but ensures after end of simulation if called
          CD_age <- 999 #redundent but ensures after end of simulation if called
        } #endifelse
        
        #All cause moratlity
        #Get a mortality age
        
        Mort_age <- mort_sample[j]
        if(Mort_age <= start_age){Mort_age <-qweibull(p = runif(n = 1,min = pweibull(q = start_age,shape = mortality_wb_a,scale = mortality_wb_b), max = 1),shape = mortality_wb_a, scale = mortality_wb_b)}
        #Ca incidnece ('original' incidence time) trumps mortality because it probability conditional on survival
        if(ca_case == 1 & Mort_age <= ca_incidence_age){Mort_age <-qweibull(p = runif(n = 1,min = pweibull(q = CD_age,shape = mortality_wb_a,scale = mortality_wb_b), max = 1),shape = mortality_wb_a, scale = mortality_wb_b)}
        if(Mort_age >= time_horizon){Mort_age <- 99.99}
        cancer_diagnostic[7] <- c(Mort_age)
        
        
        #Other individual variables
        age <- start_age
        interval_ca <- 0
        screen_detected_ca <- 0

        #DES#
        #Time_to* variables
        Time_to_screen <- screen_times[1] - age #select the current next screen age and subtract age
        Time_to_death <- Mort_age - age #time to death from current age
        Time_to_CD <- CD_age - age  #Time to clinical detection
        
        #triple While loop condition check if abosrbing death, screen_detected or interval ca event has occured
        #update age at the end of each iteration
        
        while ((age < Mort_age) && (interval_ca == 0) && (screen_detected_ca == 0)){
          #events pre-diagnosis
          Event_list <- c(Time_to_screen,Time_to_death,Time_to_CD)
          Event_place <- which.min(Event_list) # pick the nearest event in time
          Next_event_time <- Event_list[Event_place] # the time to nearest event
          #Screening event
          if(Event_place == 1){
            screen_count <- screen_count + 1 #update screen count if screening event
            costs <- costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_screen
            if(screen_count == length(screen_times)){lastscreen_count <- 1} #count number recieving last screen
            if(US_screening == 1){US_count <- US_count + 1
            costs <- costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_US
            US_costs <- US_costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_US}
            if(MRI_screening == 1){MRI_count <- MRI_count + 1
            costs <- costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_MRI
            MRI_costs <- MRI_costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_MRI}
            
            
            #Ca size to screening function return result
            if (Event_place == 1 && ca_case ==1){
              #determine if tumour is present at this time
              t <- (age+Next_event_time) - gen_age # time from cancer genesis
              if(t > 0){
                
                Ca_size <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*grow_rate_i*t))^4 #tumour volume at time t
                Ca_size <- 2*(Ca_size/(4/3*pi))^(1/3)
                
              screen_result <- screening_result(Ca_size,VDG,MRI_screening,US_screening)
              
              if(screen_result[1] == 1){
              screen_detected_ca <-1
              
             # cancer_counter <- cancer_counter + 1
              cancer_diagnostic[1] <- c((age+Time_to_screen))
              cancer_diagnostic[3] <- c(Ca_size)
              cancer_diagnostic[4] <- c(1)
              cancer_diagnostic[5] <- c(screen_result[4])
              cancer_diagnostic[6] <- c(screen_result[3])
              cancer_diagnostic[10] <- c(screen_count)
              
              incidence_age_record = age+Time_to_screen
              costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_follow_up #follow-up costs of + test
              costs_follow_up = costs_follow_up + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_follow_up
              }
              if(screen_result[1] == 1 && screen_count == 1){sdfirst_cancer <-1} #ca detected in first screen
              if(screen_result[1] == 1 && screen_count == length(screen_times)){sdlast_cancer <-1} #ca detected on last screen we are interested in for matrix
              }else{screen_detected_ca <- 0}
              
              }else{screen_detected_ca <- 0}#end if screening and cancer case
            if(Event_place == 1 && screen_detected_ca == 0 && runif(1,0,1)<recall_rate){recall_count <- recall_count+1} # check if there is a false-positive screen
          } #end screening event if
          
          #CD Event
          if(Event_place == 3){
          interval_ca <-1
          incidence_age_record = age+Time_to_CD
          costs <- costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_follow_up
          
          #cancer_counter <- cancer_counter + 1
          cancer_diagnostic[1] <- c((age+Time_to_CD))
          cancer_diagnostic[3] <- c(CD_size)
          
          } 
          
          
          #Section if cancer is detected, clinically or by screening
          if(screen_detected_ca == 1 || interval_ca == 1){
            age <- age + Next_event_time #update age if cancer
            #Assign an NPI category based on tumour size
            if(interval_ca == 1){Ca_size <- CD_size}
            #If screen detected then already set in Ca_size
            NPI_cat <- NPI_by_size(Ca_size, screen_detected_ca)
            #update the NPI counters
            if(NPI_cat == 1){NPI1_counter = NPI1_counter+1
            costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_NPI_Good}
            if(NPI_cat == 2){NPI2_counter = NPI2_counter+1
            costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_NPI_Moderate}
            if(NPI_cat == 3){NPI3_counter = NPI3_counter+1
            costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_NPI_Poor}
            if(NPI_cat == 4){NPI4_counter = NPI4_counter+1
            costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_DCIS}
            if(NPI_cat == 5){NPI5_counter = NPI5_counter+1
            costs = costs + (1/(1+discount_costs)^(Next_event_time+age-start_age))*cost_metastatic}
            #Generate a cancer specific survival time, acconting for competing risks
            Ca_mort_age <- Ca_survival_time(NPI_cat,Mort_age,age, CD_age)
            Mort_age <- Ca_mort_age
            
            if(Mort_age >= time_horizon){Mort_age <- 99.99}
            cancer_diagnostic[9] <- c(Mort_age)
            cancer_diagnostic[2] <- c(NPI_cat) 
            
            
          }else{age <- age + Next_event_time #update age if no cancer
          }
          
          #update times for next event
          if(screen_count < length(screen_times)){Time_to_screen <- screen_times[screen_count+1] - age}else{Time_to_screen <- 101} #when screen times runs out set time to age 101
          Time_to_death <- Mort_age - age 
          Time_to_CD <- CD_age - age
          
        } #while1 end
        if((screen_detected_ca+interval_ca) == 0){cancer_diagnostic[1] <- Mort_age} # recorded age is age of death or cancer incidence
        
        #all ca/screen counters
        screen_detected_count <- screen_detected_count + screen_detected_ca
        screen_counter <- screen_counter + screen_count
        US_counter <- US_counter + US_count
        MRI_counter <- MRI_counter + MRI_count
        #FP recalls
        recall_counter <- recall_counter + recall_count
        #first screen detected ca counter
        sdfirst_counter <- sdfirst_counter + sdfirst_cancer
        #last ca/screen counters
        sdlast_counter <- sdlast_counter + sdlast_cancer
        lastscreen_counter <- lastscreen_counter + lastscreen_count
        #Life-year counter
        LY_counter <- LY_counter + (Mort_age-40)
        #QALY counter
        QALY_length <- ceiling(Mort_age)-start_age
        if(QALY_length<1){QALY_length <-1}
        if(QALY_length>time_horizon-start_age){QALY_length <-time_horizon-start_age}
        QALY_vect <- rep(0,QALY_length)
        for (y in 1:length(QALY_vect)){
          QALY_vect[y] <- utility_ages[match((ceiling((start_age+y)/5)*5),utility_ages[,1]),2]
          QALY_vect[y] <- QALY_vect[y]*(1/(1+discount_health)^y) # apply discounting
        }
        if (incidence_age_record > 0){
          QALY_vect[floor(incidence_age_record)-start_age] <- utility_NPI_cat_y1[NPI_cat]*QALY_vect[floor(incidence_age_record)-start_age]}
        
        if(incidence_age_record > 0 && ceiling(if(Mort_age<100){Mort_age}else{100}) > incidence_age_record+1){
          for (y in (incidence_age_record+1):min((incidence_age_record+9),ceiling(if(Mort_age<100){Mort_age}else{100}))){
            QALY_vect[y-start_age] <- QALY_vect[y-start_age]*utility_NPI_cat_follow[NPI_cat]
            
          }
        }
        
        QALY_counter <- QALY_counter + sum(QALY_vect,na.rm = TRUE) 
        #end of j sample loop
      }
      
      ###########Record the number of screens and screen_detected cancers across all simulated histories
      ######Don't need these counters for parellel version########
      #total_screens <- total_screens + screen_counter
      #total_cancers_detected <- total_cancers_detected + screen_detected_count
      #total_costs <- total_costs + costs
      #total_US_costs <- total_US_costs + US_costs
      #total_MRI_costs <- total_MRI_costs + MRI_costs
      #total_costs_follow_up <- total_costs_follow_up + costs_follow_up
      #(total_screens[i]*cost_screen+ (recall_counter*cost_follow_up) + (recall_counter*cost_biop*biopsy_rate) + screen_detected_count*(cost_follow_up+cost_biop)) + (NPI1_counter*cost_NPI_Good) + (NPI2_counter*cost_NPI_Moderate) + (NPI3_counter*cost_NPI_Poor) + (NPI4_counter*cost_DCIS) + (NPI5_counter*cost_metastatic)
      #total_life_years <- total_life_years + LY_counter
      #total_US <- total_US  + US_counter
      #total_MRI <- total_MRI + MRI_counter
      #total_QALYs <- total_QALYs + QALY_counter
      #Returned result for foreach loop
      c(LY_counter, QALY_counter, costs, screen_counter, (screen_detected_ca+interval_ca), cancer_diagnostic)
      
      #End of individual risk factors sampling loop
    }
    #results <- rbind(results, c(total_life_years, total_QALYs, total_costs, total_screens))
    #results <- rbind(results, c(total_screens, total_cancers_detected, total_costs, total_US_costs, total_MRI_costs, total_costs_follow_up, total_life_years, total_US, total_MRI, total_QALYs))
results <- data.frame(results)
names(results)[1] <- 'LY'
names(results)[2] <- 'QALY'
names(results)[3] <- 'cost'
names(results)[4] <- 'screens'
names(results)[5] <- 'cancer'
names(results)[6] <- 'age'
names(results)[7] <- 'NPI'
names(results)[8] <- 'ca_size'
names(results)[9] <- 'screen_detected'
names(results)[10] <- 'US'
names(results)[11] <- 'MRI'
names(results)[12] <- 'Initial_mortality'
names(results)[13] <- 'CD_age'
names(results)[14] <- 'Postca_mortality'
names(results)[15] <- 'screening_round'

#directory to save 1 million sets of case histories and name of files
save(results,file = paste("",ii,".Rdata",sep = "")) 

  } #end 1 mill sim loop
  #results #see result if parellel version
  #save results
  #see results
  merged_result <- matrix(0,nrow = 10,ncol = 5)
  for (i in 1:10){
    #name of saved files needed
    load(paste("",i,".Rdata",sep = ""))
    merged_result[i,1] <- mean(results$QALY)
    merged_result[i,2] <- mean(results$cost)
    merged_result[i,3] <- mean(results$screens)
    merged_result[i,4] <- mean(results$cancer)
    merged_result[i,5] <- mean(results$screen_detected)
  }
  #store main outputs as csv
write.csv(merged_result,file = "")  
