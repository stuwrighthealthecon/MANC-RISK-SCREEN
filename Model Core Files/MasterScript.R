#Install required packages
install.packages("doParallel")
install.packages("MASS")
install.packages("dqrng")
install.packages("compiler")
install.packages("tidyverse")
install.packages("iterators")

#Run required packages
library("doParallel")
library("MASS")
library("dqrng")
library("compiler")
library("tidyverse")
library("iterators")

#####Choose screening programme and related parameters##########

#Set the screening strategy: 1=PROCAS, 2=Risk tertiles, 3=3 yearly, 4=2 yearly,
#5=5 yearly, 6=2 rounds at 50 and 60 (10 yearly), 7=Low risk (5 yearly),
#8=Low risk (6 yearly),#9=Fully stratified screening programmes
#Other num=no screening
screen_strategy<-0

#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Generate new sample? 1=YES, any other number NO
gensample<-0

#Deterministic (0) or Probabilistic Analysis (1)
PSA=0

#Standard (0) or wide (1) distributions for PSA
#Wide intervals recommended for generating data to predict GAM model
intervals=0

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Set loop numbers
inum<-3000000 #Individual women to be sampled
jnum<-1 #Lifetimes to be simulated per woman
mcruns<-1 #Monte Carlo runs used if PSA switched on
chunks<-10 #Number of chunks to split inum into for faster running time
seed<-set.seed(1) #Set seed for random draws

#Register number of cores for foreach loop
numcores<-16
registerDoParallel(cores=numcores)

#Load file containing required functions for the model
source(file="MANC_RISK_SCREEN_functions Version 1.1.R")
source(file="risksample function.R")

#################################Define baseline parameters####################

#Screening uptake
uptakefirstscreen<-0.605 #Uptake for the first screen
uptakeotherscreen<-0.852 #Uptake if woman has attended >=1 screen
uptakenoscreen<-0.191 #Uptake if woman has not previously attended screening

#Uptake for risk stratification
risk_uptake<-1 #Uptake for risk prediction
risk_feedback<-1 #Uptake for attending risk feedback consultation
screen_change<-1 #Uptake for changing screening intervals based on risk

#Age of an individual at start of simulation
start_age<-38 #Default is 38 to give time for tumours to develop pre-screening

#Set time horizon
time_horizon<-100

#Set health and cost discount rates
discount_health<-0.035
discount_cost<-0.035

##############Set clinical parameters for screening########################

#Set proportion of cancers in screening age range detected by screen
prop_screen_detected<-0.431

#Set parameters of a Weibull survival curve to represent all cause mortality
acmmortality_wb_a<-7.937
acmmortality_wb_b<-86.788

#Set parameters for all cause mortality following breast cancer
gamma_survival_1<-exp(-5.462) #Exponential distribution scale parameter stage 1
gamma_survival_2<-exp(-3.814) #Exponential distribution scale parameter stage 2
gamma_survival_3<-exp(-2.723) #Exponential distribution scale parameter stage 3
gamma_stage <- c(gamma_survival_1,gamma_survival_2,gamma_survival_3)

#Read in distribution of cancer incidence by age
Incidence_Mortality<-read.csv("Incidence_Mortality_ONS2.csv")

#Set metastatic cancer probabilities by age
metastatic_prob <- data.frame(c(25,35,45,55,65,75,85),
                              c(0.046218154,0.086659039,0.109768116,0.127099924,
                                0.142505975,0.159837783,1.73E-01))

#Create matrix of probability of cancer stage by cancer size
stage_by_size_mat<-data.frame("v1"=c(0.383,0.567,0.611,0.557,0,0),
                              "v2"=c(0.033,0.111,0.180,0.208,0.723,0.563),
                              "v3"=c(0.058,0.057,0.089,0.147,0.206,0.351),
                              "v5"=c(0.525,0.265,0.120,0.088,0.071,0.086))

#Set mean and sd of tumour doublings at clinical detection
clin_detection_m <- 6.5 
clin_detection_sd <- 0.535

#Set mean and sd of tumour doublings at screen detection
screen_detection_m <- 6.12 
screen_detection_sd <- 0.96 

#Mammography with sensitivity conditional on tumour diameter parameters
beta1 <- 1.47 
beta2 <- 6.51

#Mammography sensitivity by volpara density grade
VDG_interval<-c(4.5,7.5,15.5)
Sen_VDG <- c(0.85,0.776,0.695,0.61)
Sen_VDG_av <- 0.757

#Supplemental screening sensitivity parameters
Mammo_cdr <- 4.2 #Cancer detection rate per 1000 high dense screens
MRI_cdr <- 5 #CDR for MRI in Mammogram negative women (incremental)
US_cdr <- 3 #CDR for US in Mammogram negative women (incremental)

#Set tumour growth rate parameters
log_norm_mean <- 1.07
log_norm_sd <- 1.31
max_size <- 128 #Maximum size of tumours, mm diameter
start_size <- 0.25 #Starting size of tumours, diameter in mm
Vc = (4/3)*pi*(start_size/2)^3 #Volume at start
Vm = (4/3)*pi*(max_size/2)^3 #Max volume

#Metatstatic survival parameters
meta_survival_54 <- exp(-1.787) #Age <= 54
meta_survival_74 <- exp(-1.388) #Age 55-74
meta_survival_99 <- exp(-1.011) #Age 75+
metastatic_survival <- c(meta_survival_54, meta_survival_74, meta_survival_99)

#Set screening ages
screen_startage <- 50 #Starting age
screen_endage <- 70 #Ending age
high_risk_screentimes <- seq(screen_startage,screen_endage,1) #Yearly screening
med_risk_screentimes <- seq(screen_startage,screen_endage,2) #Two yearly screening
low_risk_screentimes <- seq(screen_startage,screen_endage,3) #Three yearly screening

#Set maximum screening sensitivity
sensitivity_max <- 0.95

#Risk cut-offs for different screening approaches
risk_cutoffs_procas <- c(1.5,3.5,5,8,100) #PROCAS and PROCAS Full
risk_cutoffs_tert <- c(1.946527,2.942792) #Tertiles of risk
low_risk_cut<-1.5 #Cut off in low risk only strategies

#Cancer size cut-points
ca_size_cut <- c(0.025, 5, 10, 15, 20, 30, 128) #Size cut-points for deciding stage

#######################Cost Data#########################################

cost_strat<-8.69 #Cost of risk prediction
cost_screen_base <- 62.21 #Cost of Mammography
cost_follow_up_base <- 109.4 #Cost of follow-up
cost_biop_base <- 297.98 #Cost of biopsyy
cost_DCIS_base <- 10107.80 #Cost of treating DCIS
cost_US_base <- 53.41 #Cost of ultrasound
cost_MRI_base <-117.10 #Cost of MRI

#If deterministic analysis then set costs as base costs
if(PSA==0){
  cost_DCIS<-cost_DCIS_base
  cost_screen<-cost_screen_base
  cost_follow_up <- cost_follow_up_base
  cost_biop <- cost_biop_base
  cost_US <- cost_US_base
  cost_MRI <-cost_MRI_base
}

#Set up look-up table for treatment costs
tbl <- tribble(~Yr, ~Early_18.64, ~Late_18.64, ~Diff1, ~Early_65plus, ~Late_65plus, ~Diff2,
               0, 464, 607, 143, 1086, 1324, 238,
               1, 10746, 13315, 2569, 7597, 8804, 1207,
               2, 3357, 5785, 2429, 2529, 3650, 1121,
               3, 1953, 3782, 1829, 2156, 3170, 1014,
               4, 1627, 2932, 1305, 2230, 2924, 693,
               5, 1617, 2841, 1225, 2077, 2957, 880,
               6, 1547, 2645, 1099, 2174, 2783, 609,
               7, 1394, 2618, 1225, 2063, 2903, 840,
               8, 1376, 2559, 1183, 2134, 2454, 320,
               9, 1279, 1848, 569, 2204, 2932, 728) %>%
  dplyr::select(-Diff1, -Diff2) %>%
  pivot_longer(cols      = contains("6"),
               names_to  = c("Stage", "Age"),
               names_sep = "_",
               values_to = "Cost") %>%
  group_by(Stage, Age) %>%
  mutate(DCost      = Cost - first(Cost),
         DCost.i    = DCost * 1.2524778811488, # NHSCII inflator for 2010/11-->2021/2022
         disc       = 1/1.035^(Yr-0.5),
         DCost.i.d  = DCost.i * disc,
         CDCost.i.d = cumsum(DCost.i.d),
         Yr1        = as.factor(Yr==1),
         Yr2        = as.factor(Yr==2),
         Yr3        = as.factor(Yr==3)) %>%
  filter(Yr > 0) %>%
  arrange(Stage, Age, Yr)

#Create predictive model of cost by age of diagnosis, stage, life expectancy
mod <- lm(data = tbl,
          formula = log(DCost) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age)

#Prediction matrix
tblNewDat <- crossing(Yr=1:50, Stage=c("Early", "Late"), Age=c("18.64", "65plus")) %>%
  mutate(Yr1 = as.factor(Yr==1),
         Yr2 = as.factor(Yr==2),
         Yr3 = as.factor(Yr==3))

#Generate cost predictions
tblNewDat %>%
  bind_cols(pred = mod %>% predict(newdata = tblNewDat)) %>%
  mutate(DCost.p = exp(pred)) -> tblPred

#Make lookup table for costs
tblLookup <- tblPred %>%
  filter(Yr==1) %>%
  mutate(across(c(Yr, pred, DCost.p), ~0)) %>%
  bind_rows(tblPred) %>%
  group_by(Stage, Age) %>%
  mutate(DCost.p.i    = DCost.p * 1.219312579, # NHSCII inflator for 2010/11-->2020/21
         disc         = 1/1.035^(Yr-0.5),
         DCost.p.i.d  = DCost.p.i * disc,
         CDCost.p.i.d = cumsum(DCost.p.i.d),
         StageEarly   = Stage=="Early",
         AgeYoung     = Age=="18.64") %>%
  arrange(Stage, Age, Yr) %>%
  ungroup()

##########False Positive and Overdiagnosis parameters################
recall_rate <- 0.0456 #UK recall rate
biopsy_rate <- 0.024 #Proporiton of referrals without cancer that have biopsy

#######################Utility Weights#########################################

#Set age adjusted utility values
utility_ages<-data.frame(c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
                         c(0.9383,0.9145,0.9069,0.8824,0.8639,0.8344,0.8222,0.8072,0.8041,0.779,0.7533,0.6985,0.6497,0.6497,0.6497))

#Set time independent utility decrements
#Utility of DCIS
utility_DCIS <- 1 #assumes no effect

#Set first year cancer utilities: 
utility_stage_cat_y1 <- c("stage1"=0.82/0.822, 
                          "stage2"=0.82/0.822,
                          "stage3"=0.75/0.822,
                          "Metastatic"=0.75/0.822,
                          "DCIS"=utility_DCIS)

#Set following year cancer utilities:
utility_stage_cat_follow <- c("stage1"=0.82/0.822,
                              "stage2"=0.82/0.822,
                              "stage3"=0.75/0.822,
                              "Metastatic"=0.75/0.822,
                              "DCIS"=utility_DCIS)

#########################CREATE SAMPLE OF WOMEN FOR MODEL###################
if(gensample==1){create_sample(PSA,intervals,seed)}

################Outer Individual sampling loop##############################

#Set loop to divide i loop into a number of sub-loops in case of simulation break
for (ii in 1:chunks) {
  load(paste("Risksample/risksample_",ii,".Rdata",sep = ""))
  prefix<-paste("^","X",ii,".",sep="")
  names(splitsample)<-sub(prefix,"",names(splitsample))
  
  #Assign women to risk groups based on 10yr risk if using risk-stratified approach  
  if(screen_strategy==1 | screen_strategy==9) {
    splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk,risk_cutoffs_procas)
  } else
    if(screen_strategy==2) {
      splitsample$risk_group<-1+findInterval(splitsample$tenyrrisk,risk_cutoffs_tert)
    } else
      if(screen_strategy==7 | screen_strategy==8) {
        splitsample$risk_group<-ifelse(splitsample$tenyrrisk<low_risk_cut,1,2)
      }  
  
  #Assign women to supplemental screening if switched on and criteria met 
  if(supplemental_screening==1){
    for (i in 1:length(splitsample$MRI_screen)) {
      if(splitsample[i,"VDG"]>=density_cutoff & splitsample[i,"tenyearrisk"]>=8){splitsample[i,"MRI_screen"]<1}else
        if(splitsample[i,"VDG"]>=density_cutoff & splitsample[i,"tenyearrisk"]<8){splitsample[i,"US_screen"]<-1}}}
  
  #If risk-stratified screening used then determine if each woman chooses to have
  #risk predicted, attends risk consultation, and changes interval
  if(screen_strategy==1 | screen_strategy==2 | (screen_strategy>6 & screen_strategy<10)){
    splitsample$risk_predicted<-ifelse(dqrunif(
      length(splitsample$risk_predicted),0,1)<
        c(rep(risk_uptake,length(splitsample$risk_predicted))),1,0)
    splitsample$feedback<-ifelse(splitsample$risk_predicted==1 & 
                                   dqrunif(length(splitsample$feedback),0,1)<
                                   c(rep(risk_feedback)),1,0)
    splitsample$interval_change<-ifelse(splitsample$feedback==1 & 
                                          dqrunif(length(splitsample$interval_change),0,1)<
                                          c(rep(screen_change)),1,0)
  }
  
  #Create iterator for the data.frame of women to pass to parallel processors  
  itx<-iter(splitsample,by="row")
  
  #Set counters for individual sampling loop
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
    }
    
    ############################## Set Screen times###############################
    
    #Assign screening intervals based on strategy and risk group    
    screen_times <- c(999)
    if (screen_strategy==1 & risk_data$interval_change==1) {
      if (risk_data$risk_group<4) {screen_times<-low_risk_screentimes} else
        if (risk_data$risk_group>3 & risk_data$risk_group<5) {screen_times<-med_risk_screentimes} else
          if (risk_data$risk_group>4) {screen_times<-high_risk_screentimes}
    } else if(screen_strategy==1 & risk_data$interval_change==0) {screen_times<-low_risk_screentimes}
    if(screen_strategy==2 & risk_data$interval_change==1){
      if(risk_data$risk_group==1){screen_times<-low_risk_screentimes} else
        if(risk_data$risk_group==2){screen_times<-med_risk_screentimes} else
          if(risk_data$risk_group==3){screen_times<-high_risk_screentimes}
    } else if(screen_strategy==1 & risk_data$interval_change==0) {screen_times<-low_risk_screentimes}
    if(screen_strategy==3){
      screen_times <- low_risk_screentimes
    }
    if(screen_strategy==4){
      screen_times <- med_risk_screentimes
    }
    if(screen_strategy==5){
      screen_times <- seq(screen_startage, screen_startage+(5*4),5)
    }
    if(screen_strategy==6){
      screen_times <- seq(screen_startage, screen_startage+10,10)
    }
    if(screen_strategy==7 & risk_data$interval_change==1){
      if(risk_data$risk_group==1){screen_times<-seq(screen_startage, screen_startage+(5*4),5)}
      if(risk_data$risk_group==2){screen_times<-low_risk_screentimes}
    } else if(screen_strategy==7 & risk_data$interval_change==0) {screen_times<-low_risk_screentimes}
    if(screen_strategy==8 & risk_data$interval_change==1){
      if(risk_data$risk_group==1){screen_times<-seq(screen_startage,screen_startage+(6*3),6)}
      if(risk_data$risk_group==2){screen_times<-low_risk_screentimes}
    } else if (screen_strategy==8 & risk_data$interval_change==0) {screen_times<-low_risk_screentimes}
    if(screen_strategy==9 & risk_data$interval_change==1){
      if (risk_data$risk_group==1) {screen_times<-seq(screen_startage, screen_startage+(5*4),5)} else
        if (risk_data$risk_group==2 | risk_data$risk_group==3) {screen_times<-low_risk_screentimes} else
          if (risk_data$risk_group==4) {screen_times<-med_risk_screentimes} else
            if (risk_data$risk_group==5) {screen_times<-high_risk_screentimes}
    } else if(screen_strategy==9 & risk_data$interval_change==0) {screen_times<-low_risk_screentimes}
    
    ##########################Set counters at i loop level#########################
    
    #screen-detected cancer counts
    screen_detected_count <- 0 #Cancer detected by screening
    sdfirst_counter <- 0 #Cancer found at first screen
    sdlast_counter <-0 #Cancer found at last screen
    #Count of screens
    screen_counter <- 0 #Number of Screens
    lastscreen_counter <-0 #Last screen attended
    US_counter <- 0 #Number of ultrasounds
    MRI_counter <- 0 #Number of MRIs
    #Recall count
    recall_counter <- 0 #Number of recalls
    #Total cost
    cost_counter <- 0 #Total costs
    #Total life years
    LY_counter <- 0 #Total life years
    #Total QALYs
    QALY_counter <- 0 #Total QALYs
    #Cancer stage counters
    stage1_counter <- 0 #Stage 1 cancer found
    stage2_counter <- 0 #Stage 2 cancer found
    stage3_counter <- 0 #Stage 3 cancer found
    stage4_counter <- 0 #Stage 4 cancer found
    DCIS_counter <- 0 #DCIS found
    
    #######J loop for individual experience of breast cancer screening##########
    for (j in jnum){
      
      #Set J level counters
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
      costs_follow_up <- 0 #Follow up costs
      
      #Lifetime cancer incidence
      #Determines if a cancer occurs and at what age
      if (risk_data$cancer==1){
        ca_case<-1
        
        #Determine cancer growth rate
        grow_rate_i<-risk_data$growth_rate
        
        #Determine when the cancer would be clinically diagnosed
        ca_incidence_i <- cmp_incidence_function()
        ca_incidence_age <- ca_incidence_i[1]
        
        #Determine size at clinical detection age
        CD_size <- ca_incidence_i[4]#tumour diameter at CD
        
        #The detection age is either the age at clinical detection 
        #or a formula is applied to determine the age at screen 
        #detection
        if(ca_incidence_i[2] ==1){CD_age <- ca_incidence_i[1]} else
          CD_age <- ca_incidence_i[1] + ((log((Vm/Vc)^0.25-1)-
          log((Vm/((4/3)*pi*(ca_incidence_i[4]/2)^3))^0.25-1))/(0.25*grow_rate_i)) - 
          ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[3]/2)^3))^0.25-1))/(0.25*grow_rate_i))
        cancer_diagnostic[8] <- c(CD_age)
        
        #Calculate tumour genesis age
        t_gen <- ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(CD_size/2)^3))^0.25-1))/(0.25*grow_rate_i)) #Calculate time to get to clinical detection size
        gen_age <- CD_age - t_gen
      } else {
        ca_case <- 0
        ca_incidence_age <- 999 #Redundant but ensures after end of simulation if called
        CD_age <- 999 #Redundant but ensures after end of simulation if called
      }
      
      #Get an all-cause mortality age and make sure this is greater than start
      #age and cancer incidence age
      Mort_age <- risk_data$life_expectancy
      
      #If cancer occurs after age of death, re-draw age of death
      if(ca_case == 1 & Mort_age <= ca_incidence_age){Mort_age <-qweibull(
        p = dqrunif(n = 1,min = pweibull(
          q = CD_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),max = 1),
        shape = acmmortality_wb_a, scale = acmmortality_wb_b)}
      if(Mort_age >= time_horizon){Mort_age <- 99.99}
      cancer_diagnostic[7] <- c(Mort_age)
      
      #Other individual variables
      age <- start_age
      interval_ca <- 0 #Cancer clinically detected
      screen_detected_ca <- 0 #Cancer screen detected
      
      ##############################DES COMPONENT ###################################
      
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
        current_discount<-(1/((1+discount_cost)^(Next_event_time+age-screen_startage)))
        
        #Open screening event
        if(Event_place == 1){
          
          #Check if woman attends screen
          if (screen_count==0 & missed_screen==0 & dqrunif(1,0,1)>uptakefirstscreen |
              screen_count==0 & missed_screen>0 & dqrunif(1,0,1)>uptakenoscreen|
              screen_count>0 & dqrunif(1,0,1)>uptakeotherscreen) {missed_screen<-missed_screen+1}else{
                
                #Woman attends screen    
                screen_count<-screen_count+1
                
                #Add cost of a mammography and risk prediction if first screen for relevant strategies
                costs<-costs+(cost_screen*current_discount)
                if(screen_count==1 & screen_strategy<3 & risk_data$risk_predicted==1){costs<-costs+(cost_strat*current_discount)}
                if(screen_count==1 & screen_strategy==7 & risk_data$risk_predicted==1){costs<-costs+(cost_strat*current_discount)}
                if(screen_count==1 & screen_strategy==8 & risk_data$risk_predicted==1){costs<-costs+(cost_strat*current_discount)}
                if(screen_count==1 & screen_strategy==9 & risk_data$risk_predicted==1){costs<-costs+(cost_strat*current_discount)}
                if(screen_count == length(screen_times)){lastscreen_count <- 1}
                
                #Add costs of supplemental screening if relevant
                if(risk_data$US_screen == 1){US_count <- US_count + 1
                costs <- costs + (cost_US*current_discount)
                US_costs<-US_costs+(cost_US*current_discount)}
                if(risk_data$MRI_screen == 1){MRI_count <- MRI_count + 1
                costs <- costs + (cost_MRI*current_discount)
                MRI_costs <- MRI_costs + (cost_MRI*current_discount)}
                
                #If the next event is a screen and a cancer is present:
                if (Event_place == 1 && ca_case ==1){
                  
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
                      cancer_diagnostic[3] <- c(Ca_size)
                      cancer_diagnostic[4] <- c(1)
                      cancer_diagnostic[5] <- c(screen_result[4])
                      cancer_diagnostic[6] <- c(screen_result[3])
                      cancer_diagnostic[10] <- c(screen_count)
                      incidence_age_record = age+Time_to_screen
                      
                      #Add cost of diagnostic follow up
                      costs = costs + (cost_follow_up*current_discount) 
                      costs_follow_up = costs_follow_up + (cost_follow_up*current_discount)
                    }
                    if(screen_result[1] == 1 && screen_count == 1){sdfirst_cancer <-1} #ca detected in first screen
                    if(screen_result[1] == 1 && screen_count == length(screen_times)){sdlast_cancer <-1} #ca detected on last screen 
                  } else{screen_detected_ca <- 0} 
                } else{screen_detected_ca <- 0} 
                
                #If a cancer is not found does a false-positive occur?
                if(Event_place == 1 && screen_detected_ca == 0 && dqrunif(1,0,1)<recall_rate){
                  recall_count <- recall_count+1
                  
                  #Add costs of false-positive recall
                  costs=costs+(cost_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)
                  costs_follow_up=costs_follow_up+(costs_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)}
              }} #End screening event
        
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
          
          #Record the stage
          if(stage_cat == 1){stage1_counter = stage1_counter+1}
          if(stage_cat == 2){stage2_counter = stage2_counter+1}
          if(stage_cat == 3){stage3_counter = stage3_counter+1}
          if(stage_cat == 4){stage4_counter = stage4_counter+1}
          if(stage_cat == 5){DCIS_counter = DCIS_counter+1
          
          #Add the cost of DCIS
          costs = costs + (cost_DCIS*current_discount)}
          
          #Generate a cancer specific survival time, accounting for competing risks
          Ca_mort_age <- cmp_ca_survival_time(stage_cat,Mort_age,age,ca_incidence_age)
          
          #Reduce age of death if cancer causes woman to die earlier
          if(Ca_mort_age<Mort_age){Mort_age<-Ca_mort_age}
          
          #Set up variables to look up treatment costs
          if(stage_cat<3){iStage<-"Early"} else {iStage<-"Late"}
          if(age<65){iAge<-"18.64"} else {iAge<-"65plus"}
          
          #If deterministic analysis then look up a treatment cost for the cancer
          if(PSA==0){
            if(stage_cat <5){costs<-costs+(as.numeric(fnLookupBase(iStage,iAge,min(c(round(Mort_age-age),50)))*current_discount))}
          } else {
            #If PSA analysis then look up treatment cost and apply cost variation
            if(stage_cat <5){costs<-costs+((1+risk_data$PSA_costvar)*as.numeric(fnLookupBase(iStage,iAge,min(c(round(Mort_age-age),50)))*current_discount))}
          }         
          
          #Record age of death and stage of cancer
          cancer_diagnostic[9] <- c(Mort_age)
          cancer_diagnostic[2] <- c(stage_cat) 
          
        }else{age <- age + Next_event_time #Update age if no cancer
        }
        
        #Update times for next event
        if(screen_count+missed_screen < length(screen_times)){Time_to_screen <- screen_times[screen_count+1] - age}else{Time_to_screen <- 101} #when screen times runs out set time to age 101
        Time_to_death <- Mort_age - age 
        Time_to_CD <- CD_age - age
        
      } #End first while loop
      if((screen_detected_ca+interval_ca) == 0){cancer_diagnostic[1] <- Mort_age} # Recorded age is age of death or cancer incidence
      
      #Update all ca/screen counters
      screen_detected_count <- screen_detected_count + screen_detected_ca
      screen_counter <- screen_counter + screen_count
      US_counter <- US_counter + US_count
      MRI_counter <- MRI_counter + MRI_count
      
      #Update false-positive recalls
      recall_counter <- recall_counter + recall_count
      
      #Update first screen detected ca counter
      sdfirst_counter <- sdfirst_counter + sdfirst_cancer
      
      #Update last ca/screen counters
      sdlast_counter <- sdlast_counter + sdlast_cancer
      lastscreen_counter <- lastscreen_counter + lastscreen_count
      
      #Update Life-year counter
      LY_counter <- LY_counter + (Mort_age-start_age)
      
      #QALY counter
      #Set up a QALY vector of length equal to life years
      QALY_length <- ceiling(Mort_age)-(screen_startage-1)
      
      #If less than 1 life year lived, set length to 1
      if(QALY_length<1){QALY_length <-1}
      
      #Ensure people don't live past end of time horizon 
      if(QALY_length>time_horizon-screen_startage){QALY_length <-time_horizon-screen_startage}
      
      #Fill QALY vector with 0's
      QALY_vect <- rep(0,QALY_length)
      
      #Fill QALY vector with discounted age related utility values
      for (y in 1:length(QALY_vect)){
        QALY_vect[y] <- (utility_ages[match((ceiling(((screen_startage-1)+y)/5)*5),utility_ages[,1]),2])*(1/(1+discount_health)^y)
        QALY_vect[QALY_length]<-QALY_vect[QALY_length]*(1-(ceiling(Mort_age)-Mort_age))
      }
      #If cancer occurs then fill QALY vector with discounted cancer utilities from incidence age
      #NB this code accounts for partial years spent in different health states
      if (incidence_age_record > 0){
        QALY_vect[floor(incidence_age_record)-screen_startage] <- utility_stage_cat_y1[stage_cat]*QALY_vect[floor(incidence_age_record)-screen_startage]*(1-(incidence_age_record-floor(incidence_age_record)))}
      if(incidence_age_record>0 & Mort_age-incidence_age_record>1){
        QALY_vect[(floor(incidence_age_record)-screen_startage)+1]<-(utility_stage_cat_y1[stage_cat]*QALY_vect[(floor(incidence_age_record)-screen_startage)+1]*(incidence_age_record-floor(incidence_age_record)))+
          (utility_stage_cat_follow[stage_cat]*QALY_vect[(floor(incidence_age_record)-screen_startage)+1]*(1-(incidence_age_record-floor(incidence_age_record))))}
      if(incidence_age_record > 0 && ceiling(if(Mort_age<100){Mort_age}else{100}) > incidence_age_record+2){
        for (y in (incidence_age_record+2):min((incidence_age_record+8),ceiling(if(Mort_age<100){Mort_age}else{100}))){
          QALY_vect[y-screen_startage] <- QALY_vect[y-screen_startage]*utility_stage_cat_follow[stage_cat]
        }
      }
      
      #Record total QALYs for J loop
      QALY_counter <- QALY_counter + sum(QALY_vect,na.rm = TRUE)
    } #end j loop
    
    #If deterministic analysis then record outputs
    if(PSA==0){
      c(QALY_counter, costs, screen_counter,cancer_diagnostic[8],(screen_detected_ca+interval_ca),screen_detected_ca, screen_strategy,risk_data$growth_rate,LY_counter-(screen_startage-start_age),cancer_diagnostic)}else{
        #If PSA then record outputs + monte carlo draws
        as.numeric(c(QALY_counter, costs, screen_counter,cancer_diagnostic[8],(screen_detected_ca+interval_ca),screen_detected_ca,screen_strategy,risk_data$growth_rate,LY_counter-(screen_startage-start_age), c(risk_data[15:40])))
      }
  }
  
  #Create a results data.frame
  results <- data.frame(results)
  names(results)[1] <- 'QALY'
  names(results)[2] <- 'Cost'
  names(results)[3] <- 'Screens'
  names(results)[4] <- "Cancer Diagnosed Age"
  names(results)[5] <- "Cancer"
  names(results)[6] <- "screen detected"
  names(results)[7] <-"alternative"
  names(results)[8] <- "Growth rate"
  names(results)[9] <- "Life Years"
  
  #If PSA add additional columns for Monte Carlo draws
  if(PSA==1){
    names(results)[10:35]<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                             "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                             "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                             'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                             "PSA_log_norm_mean","PSA_log_norm_sd","PSA_cost_strat","PSA_costvar",
                             "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                             "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","mcid")
  }
  
  #Save results from this chunk as an Rdata file
  if(PSA==0){
    save(results,file = paste("Deterministic results/Determ_",screen_strategy,"_",ii,".Rdata",sep = ""))}else{
      save(results,file = paste("PSA results/PSA_",screen_strategy,"_",ii,".Rdata",sep = "")) 
    }
  
  #Print simulation progress
  print(paste(ii*10,"%"))
} #End i loop

#Create summarised results
merged_result <- matrix(0,nrow = chunks,ncol = 7)
if(PSA==0){
  for (i in 1:chunks){
    #Record average outputs for each chunk and save in an excel file
    load(paste("Deterministic results/Determ_",screen_strategy,"_",i,".Rdata",sep = ""))
    results<-results %>% filter(results[,4]>50 | results[,4]==0)
    merged_result[i,1] <- mean(results[,1])
    merged_result[i,2] <- mean(results[,2])
    merged_result[i,3] <- mean(results[,3]) 
    merged_result[i,4] <- mean(results[,5])
    merged_result[i,5] <- mean(results[,6])
    merged_result[i,6] <- mean(results[,7])
    merged_result[i,7] <- mean(results[,9])
  }
  write.csv(merged_result,file = paste("Detresults_strat_",screen_strategy,".csv"))}else{
    for (i in 1:chunks){
      #Record average outputs for each chunk and save in an excel file
      load(paste("PSA results/PSA_",screen_strategy,"_",i,".Rdata",sep = ""))
      results<-results %>% filter(results[,4]>50 | results[,4]==0)
      merged_result[i,1] <- mean(results[,1])
      merged_result[i,2] <- mean(results[,2])
      merged_result[i,3] <- mean(results[,3]) 
      merged_result[i,4] <- mean(results[,5])
      merged_result[i,5] <- mean(results[,6])
      merged_result[i,6] <- mean(results[,7])
      merged_result[i,7] <- mean(results[,9])
    } 
    write.csv(merged_result,file = paste("PSAresults_strat_",screen_strategy,".csv"))
  }


