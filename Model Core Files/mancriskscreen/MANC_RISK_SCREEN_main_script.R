controls<-list("strategies"=c(0,1,2,3,4,9), #A vector of strategies to evaluate
               "gensample"=TRUE, #Whether to generate a new sample to simulate
               "MISCLASS"=FALSE, #whether to include risk misclassification in analysis
               "PREVENTATIVE_DRUG"=FALSE,#whether to include chemoprevention in analysis
               "PSA"=FALSE, #whether to conduct a probabilistic sensitivity analysis
               "intervals"=FALSE, #whether to conduct a PSA with wide intervals for GAM estimations
               "desired_cases"=10000, #apprximate number of cancer cases required in simulation
               "chunks"=10, #number of chunks to divide analysis into
               "mcruns"=1, #number of monte carlo runs in PSA/intervals
               "numcores"=16,
               "install"=FALSE) #set number of cores for parallel processing
               
DO_INSTALL <- controls$install

if (DO_INSTALL){
  #Install required packages
  install.packages("doParallel")
  install.packages("MASS")
  install.packages("dqrng")
  install.packages("compiler")
  install.packages("tidyverse")
  install.packages("iterators")
}

MISCLASS <- controls$MISCLASS # Set to TRUE to include impact of errors in risk prediction in model
PREVENTATIVE_DRUG <- controls$PREVENTATIVE_DRUG # Set to TRUE to simulate preventative drugs

# Add specifiers for output files
det_output_path <- "Deterministic results/"
psa_output_path <- "PSA results/"
if (MISCLASS & PREVENTATIVE_DRUG){
  dir.create("Deterministic results/misclassification_and_preventative_drug",
             showWarnings = FALSE)
  det_output_path <- "Deterministic results/misclassification_and_preventative_drug/"
  dir.create("PSA results/misclassification_and_preventative_drug",
             showWarnings = FALSE)
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
}else{
  if (MISCLASS){
    dir.create("Deterministic results/misclassification",
               showWarnings = FALSE)
    det_output_path <- "Deterministic results/misclassification/"
    dir.create("PSA results/misclassification",
               showWarnings = FALSE)
    psa_output_path <- "PSA results/misclassification/"
  }
  if (PREVENTATIVE_DRUG){
    dir.create("Deterministic results/preventative_drug",
               showWarnings = FALSE)
    det_output_path <- "Deterministic results/preventative_drug/"
    dir.create("PSA results/misclassification_and_preventative_drug",
               showWarnings = FALSE)
    psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
  }
}

sample_fname <- "possample_"

#Run required packages
library("doParallel")
library("MASS")
library("dqrng")
library("compiler")
library("tidyverse")
library("iterators")
library("tictoc")

#####Choose screening programme and related parameters##########
tic()
#Set the screening strategy: 1=PROCAS, 2=Risk tertiles, 3=3 yearly, 4=2 yearly,
#5=5 yearly, 6=2 rounds at 50 and 60 (10 yearly), 7=Low risk (5 yearly),
#8=Low risk (6 yearly),#9=Fully stratified screening programmes
#Other num=no screening
screen_strategies<-unlist(controls$strategies)
for (r in 1:length(screen_strategies)){
screen_strategy<-screen_strategies[r]

#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Generate new sample? 1=YES, any other number NO
gensample<-ifelse((controls$gensample==TRUE & r==1),1,0)

#Deterministic (0) or Probabilistic Analysis (1)
PSA=ifelse(controls$PSA==TRUE,1,0)

#Standard (0) or wide (1) distributions for PSA
#Wide intervals recommended for generating data to predict GAM model
intervals=ifelse(controls$intervals==TRUE,1,0)

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Set loop numbers
chunks<-controls$chunks #Number of chunks to split inum into for faster running time
expected_prev <- .12
desired_cases <- controls$desired_cases
inum <- ceiling((desired_cases / expected_prev)) #Individual women to be sampled to give desired number of positive cancer cases
inum <- chunks * ceiling(inum / chunks) # Make sure number of women is divisible by number of chunks
mcruns<-controls$mcruns #Monte Carlo runs used if PSA switched on
seed<-set.seed(controls$seed) #Set seed for random draws

#Register number of cores for foreach loop
closeAllConnections()
numcores<-controls$numcores
registerDoParallel(cores=numcores)

#Load file containing required functions for the model
source(file="MANC_RISK_SCREEN_functions.R")
source(file="risksample function.R")
source(file="negsample function.R")

#################################Import baseline parameters####################

source(file="params.R")

#########################CREATE SAMPLE OF WOMEN FOR MODEL###################
if (MISCLASS){
  if(gensample==1){dir.create("Risksamplewithmisclass", showWarnings = FALSE)
  create_sample_with_misclass(PSA,intervals,seed,screen_strategy)}
}else{
  if(gensample==1){dir.create("Risksample", showWarnings = FALSE)
  cmp_create_sample(PSA,intervals,seed,screen_strategy)}
}
################Outer Individual sampling loop##############################

#Set loop to divide i loop into a number of sub-loops in case of simulation break
for (ii in 1:chunks) {
  start_time <- Sys.time()
  if (MISCLASS){
    load(paste("Risksamplewithmisclass/",sample_fname,ii,".Rdata",sep = ""))
  }else{
    load(paste("Risksample/",sample_fname,ii,".Rdata",sep = ""))
  }
  prefix<-paste("^","X",ii,".",sep="")
  names(splitsample)<-sub(prefix,"",names(splitsample))
  
  if(MISCLASS){
    
    #Assign women to risk groups based on 10yr risk if using risk-stratified approach  
    if(screen_strategy==1 | screen_strategy==9) {
      splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk_est,risk_cutoffs_procas)
    } else
      if(screen_strategy==2) {
        splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk_est,risk_cutoffs_tert)
      } else
        if(screen_strategy==7 | screen_strategy==8) {
          splitsample$risk_group<-ifelse(splitsample$tenyrrisk_est<low_risk_cut,1,2)
        } }else{
          if(screen_strategy==1 | screen_strategy==9) {
            splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk,risk_cutoffs_procas)
          } else
            if(screen_strategy==2) {
              splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk,risk_cutoffs_tert)
            } else
              if(screen_strategy==7 | screen_strategy==8) {
                splitsample$risk_group<-ifelse(splitsample$tenyrrisk<low_risk_cut,1,2) 
              }}

  
  if (PREVENTATIVE_DRUG){
    # # Add extra fields for drug:
    nsample <- nrow(splitsample)
    # Use 1, 2 coding for menopause status to match indexing for drug efficacy/uptake
    splitsample$starting_menses_status <- ifelse(dqrunif(nsample, 0, 1)
                                                 <prob_premen, 1, 2)
    splitsample$takes_drug <- logical(nsample)
    splitsample$time_taking_drug <- numeric(nsample)
    
    if (PSA==1){
      # Redraw effects for PSA, assuming Monte Carlo draws of parameters are the same for each individual
      drug_matrix_list <- redraw_drug_pars(splitsample[1,])
      risk_red <- drug_matrix_list[[1]]
      uptake <- drug_matrix_list[[2]]
      persistence <- drug_matrix_list[[3]]
    }
  }
  
  #Assign women to supplemental screening if switched on and criteria met 
  #if(supplemental_screening==1){
    #for (i in 1:length(splitsample$MRI_screen)) {
      #if(splitsample[i,"VDG"]>=density_cutoff & splitsample[i,"tenyearrisk"]>=8){splitsample[i,"MRI_screen"]<1}else
        #if(splitsample[i,"VDG"]>=density_cutoff & splitsample[i,"tenyearrisk"]<8){splitsample[i,"US_screen"]<-1}}}
  
  #Create iterator for the data.frame of women to pass to parallel processors  
  itx<-iter(splitsample,by="row")
  
  #Open i loop: Simulating individual women through the strategy
  results <- foreach(i=itx,.combine = 'rbind',.packages = c('MASS','dqrng','tidyverse')) %dopar% {
    
    #Set up record of age, size, mode of detection of each detected cancer
    cancer_diagnostic <- rep(0,10)
    
    #Select an individual woman from the data.frame
    risk_data<-as.data.frame(i)
    
    #If PSA switched on, replace base case parameter values with Monte Carlo draws
    if(PSA==1){
      
      beta1<-risk_data$PSA_beta_1
      beta2<-risk_data$PSA_beta_2
      
      log_norm_mean<-risk_data$PSA_log_norm_mean
      log_norm_sd<-risk_data$PSA_log_norm_sd
      
      gamma_survival_1<-exp(risk_data$PSA_gamma_survival_1) 
      gamma_survival_2<-exp(risk_data$PSA_gamma_survival_2) 
      gamma_survival_3<-exp(risk_data$PSA_gamma_survival_3) 
      gamma_stage <- c(gamma_survival_1,gamma_survival_2,gamma_survival_3)
      
      meta_survival_54 <- exp(risk_data$PSA_meta_survival_54) 
      meta_survival_74 <- exp(risk_data$PSA_meta_survival_74) 
      meta_survival_99 <- exp(risk_data$PSA_meta_survival_99) 
      metastatic_survival <- c(meta_survival_54, meta_survival_74, meta_survival_99)
      
      Sen_VDG<-c(risk_data$PSA_VDG1_sen,risk_data$PSA_VDG2_sen,
                 risk_data$PSA_VDG3_sen,risk_data$PSA_VDG4_sen)
      Sen_VDG_av<-mean(Sen_VDG)
      
      MRI_cdr<-risk_data$PSA_MRI_cdr
      US_cdr<-risk_data$PSA_US_cdr
      
      risk_data$growth_rate<-risk_data$cancer*qlnorm(dqrunif(1,0,1),
                                                     meanlog=log_norm_mean,
                                                     sdlog=sqrt(log_norm_sd))
      
      utility_stage_cat_y1 <- c("stage1"=risk_data$PSA_util_1to3/0.822, 
                                "stage2"=risk_data$PSA_util_1to3/0.822,
                                "stage3"=risk_data$PSA_util_1to3/0.822,
                                "Metastatic"=risk_data$PSA_util_4/0.822,
                                "DCIS"=utility_DCIS)
      
      utility_stage_cat_follow <- c("stage1"=risk_data$PSA_util_1to3/0.822, 
                                    "stage2"=risk_data$PSA_util_1to3/0.822,
                                    "stage3"=risk_data$PSA_util_1to3/0.822,
                                    "Metastatic"=risk_data$PSA_util_4/0.822,
                                    "DCIS"=utility_DCIS)
      
      cost_strat<-risk_data$PSA_cost_strat
      cost_DCIS<-cost_DCIS_base*(1+risk_data$PSA_costvar)
      cost_screen<-cost_screen_base*(1+risk_data$PSA_costscreen)
      cost_follow_up <- cost_follow_up_base*(1+risk_data$PSA_cost_follow_up)
      cost_biop <- cost_biop_base*(1+risk_data$PSA_cost_biop)
      cost_US <- cost_US_base*(1+risk_data$PSA_cost_US)
      cost_MRI <-cost_MRI_base*(1+risk_data$PSA_cost_MRI)
      cost_drug <- cost_drug_base*(1+risk_data$PSA_cost_drug)
    }
    
    ############################## Set Screen times###############################
    
    #Assign screening intervals based on strategy and risk group    
    screen_times<-cmp_set_screen_times(risk_data,screen_strategy)
    
    ##########################Set counters at i loop level#########################
    
      screen_count <- 0 #Screens attended
      missed_screen<- 0 #Screens missed
      recall_count <- 0 #Number of recalls
      sdlast_cancer <-0 #Cancer detected as last screen
      lastscreen_count <- 0 #Attended last screen
      sdfirst_cancer <- 0 #Cancer detected at first screen
      stage_cat <- 0 #Stage of diagnosed cancer
      MRI_count <- 0 #Number of MRIS
      US_count <- 0 #Number of Ultrasounds
      incidence_age_record <- 0 #Age of cancer incidence
      costs <- 0 #Total costs
      US_costs <- 0 #Ultrasound costs
      MRI_costs <- 0 #MRI costs
      drug_costs <- 0 # Drug costs
      costs_follow_up <- 0 #Follow up costs
      
      #Get an all-cause mortality age and make sure this is greater than start
      #age and cancer incidence age
      Mort_age <- risk_data$life_expectancy
      
      #Other individual variables
      age <- start_age
      interval_ca <- 0 #Cancer clinically detected
      screen_detected_ca <- 0 #Cancer screen detected
      
      ##############################DES COMPONENT CANCER ###################################
        ca_case<-1
          
          #Determine cancer growth rate
          grow_rate_i<-risk_data$growth_rate
          
          # Do incididence time based on whether patient takes preventative drug
          if (PREVENTATIVE_DRUG & risk_data$risk_group!=0){
            #Determine when the cancer would be clinically diagnosed
            ca_incidence_i <- cmp_adj_incidence_function(risk_data,
                                                         uptake,
                                                         persistence,
                                                         risk_red)
            
            # Calculate cost of drug course based on time taking
            time_taking_drug <- ca_incidence_i[[5]]
            if (cost_in_full_courses){
              prop_drug_admin <- 1.
            }else{
              prop_drug_admin <- time_taking_drug / course_length[starting_menses_status]
            }
            costs <- costs + prop_drug_admin * cost_drug[risk_data$starting_menses_status]
            drug_costs <- drug_costs + prop_drug_admin * cost_drug[risk_data$starting_menses_status]
          }else{
            #Determine when the cancer would be clinically diagnosed
            ca_incidence_i <- cmp_incidence_function(risk_data)
          }
          ca_incidence_age <- ca_incidence_i[1]
          
          #Determine size at clinical detection age
          CD_size <- ca_incidence_i[4]#tumour diameter at CD
          
          #The detection age is either the age at clinical detection 
          #or a formula is applied to determine the age at screen 
          #detection
          CD_age <- ca_incidence_i[1]
          cancer_diagnostic[8]<-CD_age
          
          #Calculate tumour genesis age
          t_gen <- ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(CD_size/2)^3))^0.25-1))/(0.25*grow_rate_i)) #Calculate time to get to clinical detection size
          gen_age <- CD_age - t_gen
          
          # #If cancer occurs after age of death, re-draw age of death
          if(Mort_age <= CD_age){Mort_age <-qweibull(
            p = dqrunif(n = 1,min = pweibull(
              q = CD_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),max = 1),
            shape = acmmortality_wb_a, scale = acmmortality_wb_b)}
          if(Mort_age >= time_horizon){Mort_age <- 99.99}
          if(CD_age>=Mort_age){CD_age<-(Mort_age-0.01)}
          
          cancer_diagnostic[7] <- c(Mort_age)
      
      Time_to_screen <- screen_times[1] - age #Select the current next screen age and subtract age
      Time_to_death <- Mort_age - age #Time to death from current age
      Time_to_CD <- CD_age - age  #Time to clinical detection
      
      #Triple While loop condition check if absorbing death, screen_detected or
      #interval ca event has occurred. Update age at the end of each iteration
      
      while ((age < Mort_age) && (interval_ca == 0) && (screen_detected_ca == 0)){
        
        #Events pre-diagnosis
        Event_list <- c(Time_to_screen,Time_to_death,Time_to_CD)
        Event_place <- which.min(Event_list) #Pick the nearest event in time
        Next_event_time <- Event_list[Event_place] #The time to nearest event
        
        #Calculate current discount rate
        current_discount<-(1/((1+discount_cost)^((Next_event_time+age-screen_startage))))
        
        #Open screening event
        if(Event_place == 1){
                
                #Woman attends screen    
                screen_count<-screen_count+1
                
                #Add cost of a mammography and risk prediction if first screen for relevant strategies
                costs<-costs+(cost_screen*current_discount)
                if(screen_count==1 & screen_strategy<3 & risk_data$risk_predicted==1 |
                   screen_count==1 & screen_strategy==7 & risk_data$risk_predicted==1 |
                   screen_count==1 & screen_strategy==8 & risk_data$risk_predicted==1 |
                   screen_count==1 & screen_strategy==9 & risk_data$risk_predicted==1){costs<-costs+(cost_strat*current_discount)}
                if(screen_count == length(screen_times)){lastscreen_count <- 1}
                
                #Add costs of supplemental screening if relevant
                #if(risk_data$US_screen == 1){US_count <- US_count + 1
                #costs <- costs + (cost_US*current_discount)
                #US_costs<-US_costs+(cost_US*current_discount)}
                #if(risk_data$MRI_screen == 1){MRI_count <- MRI_count + 1
                #costs <- costs + (cost_MRI*current_discount)
                #MRI_costs <- MRI_costs + (cost_MRI*current_discount)}
                  
                  #Determine if tumour is present at screen
                  t <- (age+Next_event_time) - gen_age
                  if (t>0){
                    
                    #Determine size of tumour
                    Ca_size <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*grow_rate_i*t))^4 #tumour volume at time t
                    Ca_size <- 2*(Ca_size/(4/3*pi))^(1/3)
                    
                    #Determine if screening detects the cancer
                    screen_result <- cmp_screening_result(Ca_size,VDG=risk_data$VDG,MRI_screening = risk_data$MRI_screen,US_screening=risk_data$US_screen)
                    
                    #If a cancer is detected add a cancer and details to the counters
                    if(screen_result[1] == 1){
                      screen_detected_ca <-1
                      cancer_diagnostic[1] <- c((age+Time_to_screen))
                      cancer_diagnostic[3:6] <- c(Ca_size,
                                                  1,
                                                  screen_result[4],
                                                  screen_result[3])
                      cancer_diagnostic[10] <- c(screen_count)
                      incidence_age_record = age+Time_to_screen
                      
                      #Add cost of diagnostic follow up
                      costs = costs + (cost_follow_up*current_discount) 
                      costs_follow_up = costs_follow_up + (cost_follow_up*current_discount)
                    }
                    if(screen_result[1] == 1 && screen_count == 1){sdfirst_cancer <-1} #ca detected in first screen
                    if(screen_result[1] == 1 && screen_count == length(screen_times)){sdlast_cancer <-1} #ca detected on last screen 
                  } else{screen_detected_ca <- 0} 
                
                #If a cancer is not found does a false-positive occur?
                if(Event_place == 1 && screen_detected_ca == 0 && dqrunif(1,0,1)<recall_rate){
                  recall_count <- recall_count+1
                  
                  #Add costs of false-positive recall
                  costs=costs+(cost_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)
                  costs_follow_up=costs_follow_up+(costs_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)}
              } #End screening event
        
        #Clinical cancer diagnosis event
        if(Event_place == 3){
          interval_ca <-1
          incidence_age_record = age+Time_to_CD
          
          #Add costs of diagnosis for clinical diagnosis
          costs <- costs + (cost_follow_up*current_discount)
          
          cancer_diagnostic[1] <- c((age+Time_to_CD))
          cancer_diagnostic[3] <- c(CD_size)
        } 
        
        #If a cancer detected clinically or by screening
        if(screen_detected_ca == 1 || interval_ca == 1){
          age <- age + Next_event_time
          if(interval_ca == 1){Ca_size <- CD_size}
          
          #Assign a stage based on tumour size
          stage_cat <- cmp_stage_by_size(Ca_size)
          
          #Add the cost of DCIS
          if(stage_cat == 5){
          costs = costs + (cost_DCIS*current_discount)}
          
          #Generate a cancer specific survival time, accounting for competing risks
          Ca_mort_age <- cmp_ca_survival_time(stage_cat,Mort_age,age,CD_age)
          
          #Set up variables to look up treatment costs
          if(stage_cat<3){iStage<-"Early"} else {iStage<-"Late"}
          if(age<65){iAge<-"18.64"} else {iAge<-"65plus"}
          
          #If deterministic analysis then look up a treatment cost for the cancer
          if(PSA==0){
            if(stage_cat <5){costs<-costs+(as.numeric(fnLookupBase(iStage,iAge,min(c(round(Ca_mort_age-age),9))))*current_discount)}
          } else {
            #If PSA analysis then look up treatment cost and apply cost variation
            if(stage_cat <5){costs<-costs+((1+risk_data$PSA_costvar)*as.numeric(fnLookupBase(iStage,iAge,min(c(round(Ca_mort_age-age),9))))*current_discount)}
          }         
          
          #Record age of death and stage of cancer
          cancer_diagnostic[9] <- min(c(Ca_mort_age), c(Mort_age))
          cancer_diagnostic[2] <- c(stage_cat) 
          
        }else{age <- age + Next_event_time #Update age if no cancer
        }
        
        #Update times for next event
        if(screen_count+missed_screen < length(screen_times)){Time_to_screen <- screen_times[screen_count+1] - age}else{Time_to_screen <- 101} #when screen times runs out set time to age 101
        Time_to_death <- Mort_age - age 
        Time_to_CD <- CD_age - age
        
      } #End first while loop
      if((screen_detected_ca+interval_ca) == 0){cancer_diagnostic[1] <- Mort_age} # Recorded age is age of death or cancer incidence
      
      #Update Life-year counter
      LY_counter <- Ca_mort_age-start_age
      
      #Record total QALYs for J loop
      QALY_counter <- sum(cmp_QALY_counter(Mort_age=Ca_mort_age,
                                           incidence_age_record,
                                           stage_cat),na.rm = TRUE)
      
        if (PREVENTATIVE_DRUG & risk_data$risk_group!=0){ # Don't model impact of drug for strategies without risk stratification
          # Decide if individual takes drugs and add cost if so
          if (dqrunif(1,0,1) < uptake[risk_data$risk_group, risk_data$starting_menses_status]){
            time_taking_drug <- min(rexp(1,
                                         rate = persistence[risk_data$risk_group,
                                                            risk_data$starting_menses_status]),
                                    course_length)
          }
          else{
            time_taking_drug <- 0
          }
          
          # Calculate cost of drug course based on time taking
          if (cost_in_full_courses){
            prop_drug_admin <- 1.
          }else{
            prop_drug_admin <- time_taking_drug / course_length[starting_menses_status]
          }
          costs <- costs + prop_drug_admin * cost_drug[risk_data$starting_menses_status]
          drug_costs <- drug_costs + prop_drug_admin * cost_drug[risk_data$starting_menses_status]
        }
      
    #If deterministic analysis then record outputs
    if(PSA==0){
      return(c(QALY_counter,
               costs,
               screen_count,
               cancer_diagnostic[8],
               (screen_detected_ca+interval_ca),
               screen_detected_ca,
               screen_strategy,
               risk_data$growth_rate,
               LY_counter,
               cancer_diagnostic[2:3],
               min(c(Ca_mort_age), c(Mort_age)),cancer_diagnostic[10]))}else{
        #If PSA then record outputs + monte carlo draws
        return(as.numeric(c(QALY_counter,
                            costs,
                            screen_count,
                            cancer_diagnostic[8],
                            (screen_detected_ca+interval_ca),
                            screen_detected_ca,
                            screen_strategy,
                            risk_data$growth_rate,
                            LY_counter,
                            c(risk_data[15:41]))))
      }
  }
  
  #Create a results data.frame
  results <- data.frame(results)
  names(results) <- c('QALY',
                      'Cost',
                      'Screens',
                      "Cancer Diagnosed Age",
                      "Cancer",
                      "screen detected",
                      "alternative",
                      "Growth rate",
                      "Life Years",
                      "Stage",
                      "Cancer Size",
                      "Death Age",
                      "Cancer Screen Number")
  
  #If PSA add additional columns for Monte Carlo draws
  if(PSA==1){
    names(results)[10:36]<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                             "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                             "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                             'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                             "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                             "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                             "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","PSA_cost_drug","mcid")
  }
  
  #Save results from this chunk as an Rdata file
  if(PSA==0){
    save(results,file = paste(det_output_path,
                              "Determ_",
                              screen_strategy,
                              "_",
                              ii,
                              ".Rdata",
                              sep = ""))}else{
      save(results,file = paste(psa_output_path, "PSA_",
                                screen_strategy,
                                "_",
                                ii,
                                ".Rdata",
                                sep = "")) 
    }
  
  #Print simulation progress
  print(paste(100*ii/chunks,"%"))
  time.now <- Sys.time()
  elapsed <- as.numeric(difftime(time.now, start_time, units = "secs"))
  cat("Chunk", ii, "took", elapsed, "seconds.\n")
} #End i loop

negsamplefn(screen_strategy,MISCLASS,PSA)

#Create summarised results
merged_result <- matrix(0,nrow = chunks,ncol = 7)
if(PSA==0){
  for (i in 1:chunks){
    #Record average outputs for each chunk and save in an excel file
    load(paste(det_output_path,
               "Determ_",
               screen_strategy,
               "_",
               i,
               ".Rdata",
               sep = ""))
    results<-results %>% filter(results[,4]>50 | results[,4]==0)
    merged_result[i,1] <- mean(results[,1])
    merged_result[i,2] <- mean(results[,2])
    merged_result[i,3] <- mean(results[,3]) 
    merged_result[i,4] <- mean(results[,5])
    merged_result[i,5] <- mean(results[,6])
    merged_result[i,6] <- mean(results[,7])
    merged_result[i,7] <- mean(results[,9])
  }
  write.csv(merged_result,file = paste("Analysis/Summary results_",
                                       "Detresults_strat_",
                                       screen_strategy,
                                       ".csv",
                                       sep=""))}else{
    for (i in 1:chunks){
      #Record average outputs for each chunk and save in an excel file
      load(paste(psa_output_path, "PSA_",
                 screen_strategy,
                 "_",
                 i,
                 ".Rdata",
                 sep = ""))
      results<-results %>% filter(results[,4]>50 | results[,4]==0)
      merged_result[i,1] <- mean(results[,1])
      merged_result[i,2] <- mean(results[,2])
      merged_result[i,3] <- mean(results[,3]) 
      merged_result[i,4] <- mean(results[,5])
      merged_result[i,5] <- mean(results[,6])
      merged_result[i,6] <- mean(results[,7])
      merged_result[i,7] <- mean(results[,9])
    } 
    write.csv(merged_result,file = paste("Analysis/Summary results_",
                                         "PSAresults_strat_",
                                         screen_strategy,
                                         ".csv",
                                         sep=""))
  }

toc()
}


