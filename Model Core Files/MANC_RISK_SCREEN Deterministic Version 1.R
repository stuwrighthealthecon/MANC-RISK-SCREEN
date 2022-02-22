#########  MANC-RISK-SCREEN #########

#This model was developed by Ewan Gray and Katherine Payne
#and validated by Stuart Wright, Gabriel Rogers, Katherine
#Payne, and Rob Hainsworth at the Manchester Centre for 
#Health Economics. The model is currently maintained by
#Stuart Wright (stuart.j.wright@manchester.ac.uk)

#This discrete event simulation model estimates the QALYs,
#costs, and clinical outcomes of women attending breast 
#cancer screening delivered using one of 6 strategies or for
#no screening

#For further details see readme file and text algorithm 

#Install required packages
install.packages("doParallel")
install.packages("MASS")
install.packages("dqrng")
install.packages("compiler")
install.packages("tidyverse")

#Run required packages
library("doParallel")
library("MASS")
library("dqrng")
library("compiler")
library("tidyverse")

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Register number of cores for foreach loop
registerDoParallel(cores=8)

#Set timer to record duration of simulation
ptm <- proc.time()

#Load file containing required functions for the model
source(file="MANC_RISK_SCREEN_functions Version 1.R")

#Set loop numbers
#To attain stable results it is recommended that inum is set
#to 10,000,000. However, this will significantly slow the 
#model
inum<-10000
jnum<-1

#####Choose screening programme and related parameters##########

#Set the screening strategy: 1=PROCAS, 2=Risk tertiles, 
#3=3 yearly, 4=2 yearly, 5=5 yearly, 
#6=2 rounds at 50 and 60 (10 yearly), 
#7=Low risk (5 yearly), 8=Low risk (6 yearly)
#Other num=no screening
screen_strategy<-1

#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Age of an individual at start of simulation
start_age<-38

#Set time horizon
time_horizon<-100

#Set health and cost discount rates
discount_health<-0.035
discount_cost<-0.035

#########################################################################

##############Set clinical parameters for screening########################

#Set proportion of cancers in screening age range detected by screen
prop_screen_detected<-0.431

#Set the mean and standard deviation of the doubling rate for tumours
screen_detection_m<-4.12
screen_detection_sd<-3.93

#Set parameters of a Weibull survival curve to represent all cause mortality
acmmortality_wb_a<-8.97
acmmortality_wb_b<-86.74

#Set parameters for all cause mortality following breast cancer
gamma_survival_3<-exp(-2.723) #exponential distribution scale parameter stage 3
gamma_survival_2<-exp(-3.814) #exponential distribution scale parameter stage 2
gamma_survival_1<-exp(-5.462) #exponential distribution scale parameter stage 1
gamma_stage <- c(gamma_survival_1,gamma_survival_2,gamma_survival_3)

#Set incidence disribution
Incidence_Mortality<-read.csv("Incidence_Mortality_ONS2.csv")

#Import synthetic dataset of breast density, 
#10 year, and lifetime breast cancer risk derived from 
#PROCAS2 study
risk_mat<-read.csv("synthetic_risk_data.csv")[,2:4]

#Set metastatic cancer probabilities by age
metastatic_prob <- data.frame(c(25,35,45,55,65,75,85),
                              c(0.046218154,0.086659039,0.109768116,0.127099924,0.142505975,0.159837783,1.73E-01))

#Set proportion of ductal carcinoma in situ (DCIS)
#detected in screening
DCIS_fraction<-0.211

#Create matrix of Nottingham Prognostic Indicator by cancer size
stage_by_size_mat<-data.frame("v1"=c(0.76,0.7,0.55,0.4,0.07,0.06),
                            "v2"=c(0.22,0.27,0.43,0.55,0.64,0.5),
                            "v3"=c(0.02,0.02,0.02,0.05,0.29,0.44))

#Set mean and sd of tumour doublings at clinical detection
clin_detection_m <- 6.5 
clin_detection_sd <- 0.535

#Set mean and sd of tumour doublings at screen detection
screen_detection_m <- 6.12 
screen_detection_sd <- 0.96 

#Mammography with sensitivity conditional on tumour diameter parameters W-F
beta1 <- 1.47 
beta2 <- 6.51

#Mammography sensitivity by volpara density grade from PREVENTICON
Sen_VDG <- c(0.85,0.776,0.695,0.61)
Sen_VDG_av <- 0.757

#Supplemental screening sensitivity parameters from CEPAC
Mammo_cdr <- 4.2 #Cancer detection rate per 1000 high dense screens Mammo CEPAC
MRI_cdr <- 5 #CDR for MRI in Mammo negative women (incremental)
US_cdr <- 3 #CDR for US in Mammo negative women (incremental)

#Set tumour growth rate parameters
log_norm_mean <- 1.07
log_norm_sd <- 1.31
max_size <- 128 #mm diameter
start_size <- 0.25 #starting size of tumours, diameter in mm
Vc = (4/3)*pi*(start_size/2)^3 #Volume at start
Vm = (4/3)*pi*(max_size/2)^3 #Max volume

#Metatstatic survival parameters
meta_survival_54 <- exp(-1.787) #age <= 49
meta_survival_74 <- exp(-1.388) #age 50-69
meta_survival_99 <- exp(-1.011) # 70-99

metastatic_survival <- c(meta_survival_54, meta_survival_74, meta_survival_99)

#Set screening ages
screen_startage <- 50
screen_endage <- 70
high_risk_screentimes <- seq(screen_startage,screen_endage,1) #Yearly
med_risk_screentimes <- seq(screen_startage,screen_endage,2) #Two yearly
low_risk_screentimes <- seq(screen_startage,screen_endage,3) #Three yearly

#Set maximum screening sensitivity
sensitivity_max <- 0.95

#Risk cut-offs for different screening approaches
risk_cutoffs_procas <- c(2,3.5,5,8,100) #procas plan
risk_cutoffs_tert <- c(2.328355,3.067665) #tertiles of risk

#Breast density cut-offs for supplemental sreening
density_cutoff <- 3 #VDG groups 3 and 4

#######################Cost Data#########################################

cost_strat<-8.17
cost_screen <- 54
cost_follow_up <- 95
cost_biop <- 290
cost_DCIS <- 8806
cost_US <- 52
cost_MRI <-114

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
  select(-Diff1, -Diff2) %>%
  pivot_longer(cols      = contains("6"),
               names_to  = c("Stage", "Age"),
               names_sep = "_",
               values_to = "Cost") %>%
  group_by(Stage, Age) %>%
  mutate(DCost      = Cost - first(Cost),
         DCost.i    = DCost * 1.219312579, # NHSCII inflator for 2010/11-->2020/21
         disc       = 1/(1+discount_health)^(Yr-0.5),
         DCost.i.d  = DCost.i * disc,
         CDCost.i.d = cumsum(DCost.i.d),
         Yr1        = as.factor(Yr==1),
         Yr2        = as.factor(Yr==2),
         Yr3        = as.factor(Yr==3)) %>%
  filter(Yr > 0) %>%
  arrange(Stage, Age, Yr)

modC <- lm(data = tbl,
           formula = (CDCost.i.d) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age)
  
##########False Positive and Overdiagnosis parameters################
recall_rate <- 0.045 #approx UK recall rate
biopsy_rate <- 0.024 #proporiton of referrals without cancer that have biopsy - Madan

#######################Utility Weights#########################################

#Set age adjusted utility values
utility_ages<-data.frame(c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
                         c(0.9383,0.9145,0.9069,0.8824,0.8639,0.8344,0.8222,0.8072,0.8041,0.779,0.7533,0.6985,0.6497,0.6497,0.6497))

#Set time independent utility decrements
#Metastatic cancer 
utility_metastatic <- 0.685/0.822
utility_DCIS <- 1 #assumes no effect

#Set first year utilities: 
#Lidgren 0.696 (mean age 57, range(28-93)), metastatic 0.685 permanent
utility_stage_cat_y1 <- c("Stage1"=0.696/0.822, 
                        "Stage2"=0.696/0.822,
                        "Stage3"=0.696/0.822,
                        "Metastatic"=utility_metastatic,
                        "DCIS"=utility_DCIS) 

#Set following year utilities:
#0.779
utility_stage_cat_follow <- c("Stage1"=0.779/0.822,
                            "Stage2"=0.779/0.822,
                            "Stage3"=0.779/0.822,
                            "Metastatic"=utility_metastatic,
                            "DCIS"=utility_DCIS)

################Outer Individual sampling loop##############################

#Set loop to divide i loop into 10 sub-loops in case of simulation break
for (ii in 1:10) {
  
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

#Open loop
results <- foreach(i=1:(inum/10),.combine = 'rbind',.packages = c('MASS','dqrng','tidyverse')) %dopar% {

#Set up record of age, size, mode of detection of each detected cancer
cancer_diagnostic <- rep(0,10)

#Draw a breast density, 10 year, and lifetime risk of cancer for the individual
risk_data<-risk_mat[sample(nrow(risk_mat),1),]
ten_year_risk<-risk_data[2]

#If risk based screening is being used then place 
#individual into a risk group
if(screen_strategy==1) {
  risk_group<-1+findInterval(ten_year_risk,risk_cutoffs_procas)
} else
if(screen_strategy==2) {
  risk_group<-1+findInterval(ten_year_risk,risk_cutoffs_tert)
}
if(screen_strategy==7 | screen_strategy==8) {
  if(ten_year_risk<1.5){risk_group<-1}
  if(ten_year_risk>=1.5){risk_group<-2}
}

#Set VDG based on breast density
if(risk_data[1]<4.5){VDG<-1} else
  if(risk_data[1]>=4.5 & risk_data[1]<7.5){VDG<-2} else
  if(risk_data[1]>=7.5 & risk_data[1]<15.5){VDG<-3} else
  if(risk_data[1]>=15.5){VDG<-4}


#Set level of supplemental screening
if(supplemental_screening==0){
  MRI_screening<-0
  US_screening<-0} else {
    if (VDG>density_cutoff){
      MRI_screening<-0
      US_screening<-0
      if(ten_year_risk>8){MRI_screening<-1} else{US_screening<-1}
    } else {
       MRI_screening<-1
       US_screening<-1
    }
  }

###############Screen times###############################

screen_times <- c(999)
if (screen_strategy==1) {
  if (risk_group<3) {screen_times<-low_risk_screentimes} else
  if (risk_group>2 & risk_group<5) {screen_times<-med_risk_screentimes} else
  if (risk_group>4) {screen_times<-high_risk_screentimes}
}
  if(screen_strategy==2){
  if(risk_group==1){screen_times<-low_risk_screentimes} else
  if(risk_group==2){screen_times<-med_risk_screentimes} else
  if(risk_group==3){screen_times<-high_risk_screentimes}
  }
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
  if(screen_strategy==7){
    if(risk_group==1){screen_times<-seq(screen_startage, screen_startage+(5*4),5)}
    if(risk_group==2){screen_times<-low_risk_screentimes}
  }
  if(screen_strategy==8){
    if(risk_group==1){screen_times<-seq(screen_startage,screen_startage+(6*3),6)}
    if(risk_group==2){screen_times<-low_risk_screentimes}
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
#Cancer stage counters
stage1_counter <- 0
stage2_counter <- 0
stage3_counter <- 0
stage4_counter <- 0
DCIS_counter <- 0

###Preload incidence, mortality and clinical detection times for j cases
mort_sample<- rweibull(n = jnum,shape = acmmortality_wb_a, scale = acmmortality_wb_b)
clin_detect_sample <- dqrnorm(n = jnum,mean = clin_detection_m,sd = clin_detection_sd)
clin_detect_sample[clin_detect_sample < 4] <- 4 #Prevent unrealistic left tail
clin_detect_sample[clin_detect_sample >= 9] <- 8.99 #Prevent unrealistic right tail

#######J loop for individual experience of breast cancer screening)
for (j in jnum){
  
#Set J level counters
screen_count <- 0
recall_count <- 0
sdlast_cancer <-0
lastscreen_count <- 0
sdfirst_cancer <- 0
stage_cat <- 0
MRI_count <- 0
US_count <- 0
incidence_age_record <- 0
costs <- 0
US_costs <- 0
MRI_costs <- 0
costs_follow_up <- 0
 
#Lifetime cancer incidence
#Determines if a cancer occurs and at what age
if(dqrunif(1,0,1)<(risk_data[3]/100)){
  ca_case<-1

#Determine cancer growth rate
grow_rate_i<-qlnorm(dqrunif(1,0,1),meanlog=log_norm_mean,sdlog=sqrt(log_norm_sd))

#Incidence age (under current programme)
ca_incidence_i <- cmp_incidence_function()
ca_incidence_age <- ca_incidence_i[1]

#Clinical detection age
CD_size <- ca_incidence_i[4]#tumour diameter at CD
  
#The detection age is either the age at clinical detection 
#or a formula is applied to determine the age at screen 
#detection
if(ca_incidence_i[2] ==1){CD_age <- ca_incidence_i[1]} else
  CD_age <- ca_incidence_i[1] + ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[4]/2)^3))^0.25-1))/(0.25*grow_rate_i)) - 
    ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[3]/2)^3))^0.25-1))/(0.25*grow_rate_i))
  cancer_diagnostic[8] <- c(CD_age)
  
#Calculate tumour genesis age
t_gen <- ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(CD_size/2)^3))^0.25-1))/(0.25*grow_rate_i)) #Calculate time to get to clinical detection size
gen_age <- CD_age - t_gen
} else {
  ca_case <- 0
  ca_incidence_age <- 999 #redundent but ensures after end of simulation if called
  CD_age <- 999 #redundent but ensures after end of simulation if called
}

#All cause moratlity
#Get a mortality age and make sure this is greater than start age and cancer incidence age
Mort_age <- mort_sample[j]
if(Mort_age <= start_age){Mort_age <-qweibull(p = dqrunif(n = 1,min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b), max = 1),shape = acmmortality_wb_a, scale = acmmortality_wb_b)}

#Ca incidence ('original' incidence time) trumps mortality
#because it is probability conditional on survival
if(ca_case == 1 & Mort_age <= ca_incidence_age){Mort_age <-qweibull(p = dqrunif(n = 1,min = pweibull(q = CD_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b), max = 1),shape = acmmortality_wb_a, scale = acmmortality_wb_b)}
if(Mort_age >= time_horizon){Mort_age <- 99.99}
cancer_diagnostic[7] <- c(Mort_age)
  
#Other individual variables
age <- start_age
interval_ca <- 0
screen_detected_ca <- 0
  
  
#####################DES COMPONENT #######################
  
Time_to_screen <- screen_times[1] - age #select the current next screen age and subtract age
Time_to_death <- Mort_age - age #time to death from current age
Time_to_CD <- CD_age - age  #Time to clinical detection
  
#triple While loop condition check if abosrbing death, 
#screen_detected or interval ca event has occured
#update age at the end of each iteration
  
while ((age < Mort_age) && (interval_ca == 0) && (screen_detected_ca == 0)){
    
  #events pre-diagnosis
  Event_list <- c(Time_to_screen,Time_to_death,Time_to_CD)
  Event_place <- which.min(Event_list) # pick the nearest event in time
  Next_event_time <- Event_list[Event_place] # the time to nearest event
  current_discount<-(1/((1+discount_cost)^(Next_event_time+age-start_age)))
    
  #Open screening event
  if(Event_place == 1){
    screen_count<-screen_count+1
    costs<-costs+(cost_screen*current_discount)
    if(screen_count==1 & screen_strategy<3){costs<-costs+(cost_strat*current_discount)}
    if(screen_count==1 & screen_strategy==7){costs<-costs+(cost_strat*current_discount)}
    if(screen_count==1 & screen_strategy==8){costs<-costs+(cost_strat*current_discount)}
    if(screen_count == length(screen_times)){lastscreen_count <- 1}
    if(US_screening == 1){US_count <- US_count + 1
    costs <- costs + (cost_US*current_discount)
    US_costs<-US_costs+(cost_US*current_discount)}
    if(MRI_screening == 1){MRI_count <- MRI_count + 1
    costs <- costs + (cost_MRI*current_discount)
    MRI_costs <- MRI_costs + (cost_MRI*current_discount)}
      
    #If the next event is a screen:
    if (Event_place == 1 && ca_case ==1){
      
    #Determine if tumour is present
      t <- (age+Next_event_time) - gen_age
      if (t>0){
        
    #Determine size of tumour
      Ca_size <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*grow_rate_i*t))^4 #tumour volume at time t
      Ca_size <- 2*(Ca_size/(4/3*pi))^(1/3)
      
    #Determine if screening detects the cancer
      screen_result <- cmp_screening_result(Ca_size,VDG,MRI_screening,US_screening)
     
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
      costs = costs + (cost_follow_up*current_discount) 
      costs_follow_up = costs_follow_up + (cost_follow_up*current_discount)
      }
      if(screen_result[1] == 1 && screen_count == 1){sdfirst_cancer <-1} #ca detected in first screen
      if(screen_result[1] == 1 && screen_count == length(screen_times)){sdlast_cancer <-1} #ca detected on last screen 
      } else{screen_detected_ca <- 0} 
      } else{screen_detected_ca <- 0} 
      
    #Does a false-positive occur?
    if(Event_place == 1 && screen_detected_ca == 0 && dqrunif(1,0,1)<recall_rate){
      recall_count <- recall_count+1
      costs=costs+(cost_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)
      costs_follow_up=costs_follow_up+(costs_follow_up*current_discount)+(biopsy_rate*cost_biop*current_discount)}
      } #End screening event

    #Clinical cancer diagnosis event
    if(Event_place == 3){
      interval_ca <-1
      incidence_age_record = age+Time_to_CD
      costs <- costs + (cost_follow_up*current_discount)
      
      cancer_diagnostic[1] <- c((age+Time_to_CD))
      cancer_diagnostic[3] <- c(CD_size)
    } 
    
#Cancer detected clinically or by screening
  if(screen_detected_ca == 1 || interval_ca == 1){
    age <- age + Next_event_time
    if(interval_ca == 1){Ca_size <- CD_size}
    
    #Assign a Stage based on tumour size
    stage_cat <- cmp_stage_by_size(Ca_size, screen_detected_ca)
    if(stage_cat == 1){stage1_counter = stage1_counter+1}
    if(stage_cat == 2){stage2_counter = stage2_counter+1}
    if(stage_cat == 3){stage3_counter = stage3_counter+1}
    if(stage_cat == 4){stage4_counter = stage4_counter+1}
    if(stage_cat == 5){DCIS_counter = DCIS_counter+1
    costs = costs + (cost_DCIS*current_discount)}

    #Generate a cancer specific survival time, accounting for competing risks
    Ca_mort_age <- cmp_ca_survival_time(stage_cat,Mort_age,age,CD_age)
    if(Ca_mort_age<Mort_age){Mort_age<-Ca_mort_age}
    
    if(stage_cat<3){iStage<-"Early"} else {iStage<-"Late"}
    if(age<65){iAge<-"18.64"} else {iAge<-"65plus"}
    if(stage_cat <5){costs=costs+(fnModPred(iStage,iAge,Mort_age-age)*current_discount)}
    
    cancer_diagnostic[9] <- c(Mort_age)
    cancer_diagnostic[2] <- c(stage_cat) 
    
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
  LY_counter <- LY_counter + (Mort_age-start_age)
  
  #QALY counter
  QALY_length <- ceiling(Mort_age)-start_age
  if(QALY_length<1){QALY_length <-1}
  if(QALY_length>time_horizon-start_age){QALY_length <-time_horizon-start_age}
  QALY_vect <- rep(0,QALY_length)
  for (y in 1:length(QALY_vect)){
    QALY_vect[y] <- (utility_ages[match((ceiling((start_age+y)/5)*5),utility_ages[,1]),2])*(1/(1+discount_health)^y)
    QALY_vect[QALY_length]<-QALY_vect[QALY_length]*(1-(ceiling(Mort_age)-Mort_age))
  }
  if (incidence_age_record > 0){
    QALY_vect[floor(incidence_age_record)-start_age] <- utility_stage_cat_y1[stage_cat]*QALY_vect[floor(incidence_age_record)-start_age]*(1-(incidence_age_record-floor(incidence_age_record)))}
  if(incidence_age_record>0 & Mort_age-incidence_age_record>1){
    QALY_vect[(floor(incidence_age_record)-start_age)+1]<-(utility_stage_cat_y1[stage_cat]*QALY_vect[(floor(incidence_age_record)-start_age)+1]*(incidence_age_record-floor(incidence_age_record)))+
                                                           (utility_stage_cat_follow[stage_cat]*QALY_vect[(floor(incidence_age_record)-start_age)+1]*(1-(incidence_age_record-floor(incidence_age_record))))}
  if(incidence_age_record > 0 && ceiling(if(Mort_age<100){Mort_age}else{100}) > incidence_age_record+2){
    for (y in (incidence_age_record+2):min((incidence_age_record+8),ceiling(if(Mort_age<100){Mort_age}else{100}))){
      QALY_vect[y-start_age] <- QALY_vect[y-start_age]*utility_stage_cat_follow[stage_cat]
    }
  }
  
  QALY_counter <- QALY_counter + sum(QALY_vect,na.rm = TRUE)
} #end j loop

c(LY_counter, QALY_counter, costs, screen_counter, (screen_detected_ca+interval_ca), cancer_diagnostic)
}
results <- data.frame(results)
names(results)[1] <- 'LY'
names(results)[2] <- 'QALY'
names(results)[3] <- 'cost'
names(results)[4] <- 'screens'
names(results)[5] <- 'cancer'
names(results)[6] <- 'age'
names(results)[7] <- 'stage'
names(results)[8] <- 'ca_size'
names(results)[9] <- 'screen_detected'
names(results)[10] <- 'US'
names(results)[11] <- 'MRI'
names(results)[12] <- 'Initial_mortality'
names(results)[13] <- 'CD_age'
names(results)[14] <- 'Postca_mortality'
names(results)[15] <- 'screening_round'

#directory to save inum/10 sets of case histories and name of files
save(results,file = paste("",ii,".Rdata",sep = "")) 

} #End 1 million simulation loop
#results #see result if parellel version
#save results
#see results
merged_result <- matrix(0,nrow = 10,ncol = 5)
for (i in 1:10){
  #name of saved files needed
  load(paste("",i,".Rdata",sep = ""))
  merged_result[i,1] <- mean(results[,2])
  merged_result[i,2] <- mean(results[,3])
  merged_result[i,3] <- mean(results[,4])
  merged_result[i,4] <- mean(results[,5])
  merged_result[i,5] <- mean(results[,9])
}
#store main outputs as csv
write.csv(merged_result,file = "results.csv")


