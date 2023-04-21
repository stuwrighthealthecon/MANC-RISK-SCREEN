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
install.packages("iterators")

#Run required packages
library("doParallel")
library("MASS")
library("dqrng")
library("compiler")
library("tidyverse")
library("iterators")
library("tictoc")

tic("100k:7 cores:PROCASFULL")

#Generate new sample?
#1=YES, any other number NO
gensample<-1

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Register number of cores for foreach loop
numcores<-7
registerDoParallel(cores=numcores)

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
#7=Low risk (5 yearly), 8=Low risk (6 yearly),
#9=Fully stratified screening programmes
#Other num=no screening
screen_strategy<-9

#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Screening uptake
uptakefirstscreen<-0.605
uptakeotherscreen<-0.852
uptakenoscreen<-0.191

#Uptake for risk stratification
risk_uptake<-1
  #0.6 #Proportion of women who want risk predicted
risk_feedback<-1
  #0.95 #Proportion of women who attend risk consultation
screen_change<-1
  #0.8 #Proportion of women with high/moderate/low risk who change screening interval

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
acmmortality_wb_a<-7.937
acmmortality_wb_b<-86.788

#Set parameters for all cause mortality following breast cancer
gamma_survival_1<-exp(-5.462) #exponential distribution scale parameter stage 1
gamma_survival_2<-exp(-3.814) #exponential distribution scale parameter stage 2
gamma_survival_3<-exp(-2.723) #exponential distribution scale parameter stage 3
gamma_stage <- c(gamma_survival_1,gamma_survival_2,gamma_survival_3)

#Set incidence disribution
Incidence_Mortality<-read.csv("Incidence_Mortality_ONS2.csv")

#Set metastatic cancer probabilities by age
metastatic_prob <- data.frame(c(25,35,45,55,65,75,85),
                              c(0.046218154,0.086659039,0.109768116,0.127099924,0.142505975,0.159837783,1.73E-01))

#Create matrix of Nottingham Prognostic Indicator by cancer size
stage_by_size_mat<-data.frame("v1"=c(0.383,0.567,0.611,0.557,0,0),
                            "v2"=c(0.033,0.111,0.180,0.208,0.723,0.563),
                            "v3"=c(0.058,0.057,0.208,0.147,0.206,0.351),
                            "v5"=c(0.525,0.265,0.120,0.088,0.071,0.086))

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
VDG_interval<-c(4.5,7.5,15.5)
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
meta_survival_54 <- exp(-1.787) #age <= 54
meta_survival_74 <- exp(-1.388) #age 55-74
meta_survival_99 <- exp(-1.011) # 75+
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
risk_cutoffs_procas <- c(1.5,3.5,5,8,100) #procas plan
risk_cutoffs_tert <- c(1.946527,2.942792) #tertiles of risk
low_risk_cut<-1.5 #cut off in low risk only strategies

#Cancer size cut-points
ca_size_cut <- c(0.025, 5, 10, 15, 20, 30, 128) #category cut-points from Kolias 1999

#######################Cost Data#########################################

cost_strat<-8.45
cost_screen <- 60.56
cost_follow_up <- 106
cost_biop <- 290
cost_DCIS <- 9840
cost_US <- 52
cost_MRI <-114

# read in and refactor data from Laudicella et al. (2016) https://doi.org/10.1038/bjc.2016.77
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
         DCost.i    = DCost * 1.219312579, # NHSCII inflator for 2010/11-->2020/21
         disc       = 1/1.035^(Yr-0.5),
         DCost.i.d  = DCost.i * disc,
         CDCost.i.d = cumsum(DCost.i.d),
         Yr1        = as.factor(Yr==1),
         Yr2        = as.factor(Yr==2),
         Yr3        = as.factor(Yr==3)) %>%
  filter(Yr > 0) %>%
  arrange(Stage, Age, Yr)

# log-linear model
mod <- lm(data = tbl,
          formula = log(DCost) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age)

# prediction matrix
tblNewDat <- crossing(Yr=1:50, Stage=c("Early", "Late"), Age=c("18.64", "65plus")) %>%
  mutate(Yr1 = as.factor(Yr==1),
         Yr2 = as.factor(Yr==2),
         Yr3 = as.factor(Yr==3))

# generate predictions
tblNewDat %>%
  bind_cols(pred = mod %>% predict(newdata = tblNewDat)) %>%
  mutate(DCost.p = exp(pred)) -> tblPred

# make lookup table
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
recall_rate <- 0.0456 #approx UK recall rate
biopsy_rate <- 0.024 #proporiton of referrals without cancer that have biopsy - Madan

#######################Utility Weights#########################################

#Set age adjusted utility values
utility_ages<-data.frame(c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100),
                         c(0.9383,0.9145,0.9069,0.8824,0.8639,0.8344,0.8222,0.8072,0.8041,0.779,0.7533,0.6985,0.6497,0.6497,0.6497))

#Set time independent utility decrements
#Metastatic cancer 
utility_DCIS <- 1 #assumes no effect

#Set first year utilities: 
utility_stage_cat_y1 <- c("stage1"=0.82/0.822, 
                          "stage2"=0.82/0.822,
                          "stage3"=0.75/0.822,
                          "Metastatic"=0.75/0.822,
                          "DCIS"=utility_DCIS)

#Set following year utilities:
utility_stage_cat_follow <- c("stage1"=0.82/0.822,
                              "stage2"=0.82/0.822,
                              "stage3"=0.75/0.822,
                              "Metastatic"=0.75/0.822,
                              "DCIS"=utility_DCIS)

#########################CREATE SAMPLE OF WOMEN FOR MODEL#######
if(gensample==1){

#Import synthetic dataset of breast density, 
#10 year, and lifetime breast cancer risk derived from 
#PROCAS2 study
risk_mat<-read.csv("synthetic_risk_data.csv")[,2:4]

#If risk based screening is being used then place 
#individual into a risk group

risk_mat[,4]<-numeric(length(risk_mat[,3]))

#Set VDG based on breast density
risk_mat[,5]<-1+findInterval(risk_mat[,1],VDG_interval)

#Breast density cut-offs for supplemental sreening
density_cutoff <- 3 #VDG groups 3 and 4

#Set up data frame of women's lifetimes to simulate
risksample<-risk_mat[sample(nrow(risk_mat),inum,replace=TRUE),]
risksample[,6:14]<-numeric(length=length(risksample[,5]))

###Preload incidence, mortality and clinical detection times for j cases
risksample[,11]<- rweibull(n = length(risksample[,10]),shape = acmmortality_wb_a, scale = acmmortality_wb_b)
for (i in 1:length(risksample[,11])) {
if(risksample[i,11] <= start_age){risksample[i,11]<-qweibull(p = dqrunif(n = 1,min = pweibull(q = start_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b), max = 1),shape = acmmortality_wb_a, scale = acmmortality_wb_b)}}

#Determine if a cancer will develop
risksample[,12]<-ifelse(dqrunif(length(risksample[,11]),0,1)<(risksample[,3]/100),1,0)

#Set number of tumour doublings when there is a cancer
risksample[,13]<-risksample[,12]*(dqrnorm(n = length(risksample[,12]),mean = clin_detection_m,sd = clin_detection_sd))

#Set growth rate for tumours
risksample[,14]<-risksample[,12]*qlnorm(dqrunif(length(risksample[,12]),0,1),meanlog=log_norm_mean,sdlog=sqrt(log_norm_sd))
names(risksample)[4:14]<-paste(c("Risk Group","VDG","MRI Screening","US Screening","Risk Predicted","Feedback","Interval Change","Life Expectancy","Cancer","Clinical Detection Size","Growth Rate"))
risksample[,15]<-(rep(1:10,times=round(length(risksample)/10)))
risksample<-risksample %>% filter(risksample[,11]>=50)

risksplit<-split(risksample,risksample[,15])

#Save risk sample in chunks
for(i in 1:10){
splitsample<-as.data.frame(risksplit[i])
save(splitsample,file = paste("Risksample/risksample_",i,".Rdata",sep=""))
}
rm(risksample)
rm(risksplit)
rm(splitsample)
}
################Outer Individual sampling loop##############################

#Set loop to divide i loop into 10 sub-loops in case of simulation break
for (ii in 1:10) {

load(paste("Risksample/risksample_",ii,".Rdata",sep = ""))
risksample<-splitsample

if(screen_strategy==1 | screen_strategy==9) {
  risksample[,4]<-1+findInterval(risksample[,2],risk_cutoffs_procas)
} else
  if(screen_strategy==2) {
    risksample[,4]<-1+findInterval(risksample[,2],risk_cutoffs_tert)
  } else
    if(screen_strategy==7 | screen_strategy==8) {
      risksample[,4]<-ifelse(risksample[,2]<low_risk_cut,1,2)
    }  

if(supplemental_screening==1){
  for (i in 1:length(risksample[,6])) {
    if(risksample[i,5]>=density_cutoff & risksample[i,2]>=8){risksample[i,6]<1}else
      if(risksample[i,5]>=density_cutoff & risksample[i,2]<8){risksample[i,7]<-1}}}

if(screen_strategy==1 | screen_strategy==2 | (screen_strategy>6 & screen_strategy<10)){
  risksample[,8]<-ifelse(dqrunif(length(risksample[,1]),0,1)<c(rep(risk_uptake,length(risksample[,1]))),1,0)
  risksample[,9]<-ifelse(risksample[,8]==1 & dqrunif(length(risksample[,1]),0,1)<c(rep(risk_feedback)),1,0)
  risksample[,10]<-ifelse(risksample[,9]==1 & dqrunif(length(risksample[,1]),0,1)<c(rep(screen_change)),1,0)
}

itx<-iter(risksample,by="row")

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
results <- foreach(i=itx,.combine = 'rbind',.packages = c('MASS','dqrng','tidyverse')) %dopar% {

#Set up record of age, size, mode of detection of each detected cancer
cancer_diagnostic <- rep(0,10)

#Draw a breast density, 10 year, and lifetime risk of cancer for the individual
risk_data<-as.numeric(i)

###############Screen times###############################

screen_times <- c(999)
if (screen_strategy==1 & risk_data[10]==1) {
  if (risk_data[4]<4) {screen_times<-low_risk_screentimes} else
  if (risk_data[4]>3 & risk_data[4]<5) {screen_times<-med_risk_screentimes} else
  if (risk_data[4]>4) {screen_times<-high_risk_screentimes}
} else if(screen_strategy==1 & risk_data[10]==0) {screen_times<-low_risk_screentimes}
  if(screen_strategy==2 & risk_data[10]==1){
  if(risk_data[4]==1){screen_times<-low_risk_screentimes} else
  if(risk_data[4]==2){screen_times<-med_risk_screentimes} else
  if(risk_data[4]==3){screen_times<-high_risk_screentimes}
  } else if(screen_strategy==1 & risk_data[10]==0) {screen_times<-low_risk_screentimes}
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
  if(screen_strategy==7 & risk_data[10]==1){
    if(risk_data[4]==1){screen_times<-seq(screen_startage, screen_startage+(5*4),5)}
    if(risk_data[4]==2){screen_times<-low_risk_screentimes}
  } else if(screen_strategy==7 & risk_data[10]==0) {screen_times<-low_risk_screentimes}
  if(screen_strategy==8 & risk_data[10]==1){
    if(risk_data[4]==1){screen_times<-seq(screen_startage,screen_startage+(6*3),6)}
    if(risk_data[4]==2){screen_times<-low_risk_screentimes}
  } else if (screen_strategy==8 & risk_data[10]==0) {screen_times<-low_risk_screentimes}
  if(screen_strategy==9 & risk_data[10]==1){
    if (risk_data[4]==1) {screen_times<-seq(screen_startage, screen_startage+(5*4),5)} else
    if (risk_data[4]==2 | risk_data[4]==3) {screen_times<-low_risk_screentimes} else
    if (risk_data[4]==4) {screen_times<-med_risk_screentimes} else
    if (risk_data[4]==5) {screen_times<-high_risk_screentimes}
  } else if(screen_strategy==9 & risk_data[10]==0) {screen_times<-low_risk_screentimes}

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

#######J loop for individual experience of breast cancer screening)
for (j in jnum){
  
#Set J level counters
screen_count <- 0
missed_screen<- 0
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
if (risk_data[12]==1){
  ca_case<-1

#Determine cancer growth rate
grow_rate_i<-risk_data[14]

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
Mort_age <- risk_data[11]

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
  current_discount<-(1/((1+discount_cost)^(Next_event_time+age-screen_startage)))
    
  #Open screening event
  if(Event_place == 1){
       if (screen_count==0 & missed_screen==0 & dqrunif(1,0,1)>uptakefirstscreen |
        screen_count==0 & missed_screen>0 & dqrunif(1,0,1)>uptakenoscreen|
        screen_count>0 & dqrunif(1,0,1)>uptakeotherscreen) {missed_screen<-missed_screen+1}else{
    screen_count<-screen_count+1
    costs<-costs+(cost_screen*current_discount)
    if(screen_count==1 & screen_strategy<3 & risk_data[8]==1){costs<-costs+(cost_strat*current_discount)}
    if(screen_count==1 & screen_strategy==7 & risk_data[8]==1){costs<-costs+(cost_strat*current_discount)}
    if(screen_count==1 & screen_strategy==8 & risk_data[8]==1){costs<-costs+(cost_strat*current_discount)}
    if(screen_count==1 & screen_strategy==9 & risk_data[8]==1){costs<-costs+(cost_strat*current_discount)}
    if(screen_count == length(screen_times)){lastscreen_count <- 1}
    if(risk_data[7] == 1){US_count <- US_count + 1
    costs <- costs + (cost_US*current_discount)
    US_costs<-US_costs+(cost_US*current_discount)}
    if(risk_data[6] == 1){MRI_count <- MRI_count + 1
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
      screen_result <- cmp_screening_result(Ca_size,VDG=risk_data[5],MRI_screening = risk_data[6],US_screening=risk_data[7])
     
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
      }} #End screening event

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
    stage_cat <- cmp_stage_by_size(Ca_size)
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
    if(stage_cat <5){costs<-costs+as.numeric(fnLookupBase(iStage,iAge,min(c(round(Mort_age-age),50)))*current_discount)}
    cancer_diagnostic[9] <- c(Mort_age)
    cancer_diagnostic[2] <- c(stage_cat) 
    
  }else{age <- age + Next_event_time #update age if no cancer
  }
    #update times for next event
    if(screen_count+missed_screen < length(screen_times)){Time_to_screen <- screen_times[screen_count+1] - age}else{Time_to_screen <- 101} #when screen times runs out set time to age 101
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
  QALY_length <- ceiling(Mort_age)-(screen_startage-1)
  if(QALY_length<1){QALY_length <-1}
  if(QALY_length>time_horizon-screen_startage){QALY_length <-time_horizon-screen_startage}
  QALY_vect <- rep(0,QALY_length)
  for (y in 1:length(QALY_vect)){
    QALY_vect[y] <- (utility_ages[match((ceiling(((screen_startage-1)+y)/5)*5),utility_ages[,1]),2])*(1/(1+discount_health)^y)
    QALY_vect[QALY_length]<-QALY_vect[QALY_length]*(1-(ceiling(Mort_age)-Mort_age))
  }
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

print(paste(ii*10,"%"))
} #End 1 million simulation loop
#results #see result if parellel version
#save results
#see results
merged_result <- matrix(0,nrow = 10,ncol = 5)
for (i in 1:10){
  #name of saved files needed
  load(paste("",i,".Rdata",sep = ""))
  results<-results %>% filter(results[,13]>50 | results[,13]==0)
  merged_result[i,1] <- mean(results[,2])
  merged_result[i,2] <- mean(results[,3])
  merged_result[i,3] <- mean(results[,4]) 
  merged_result[i,4] <- mean(results[,5])
  merged_result[i,5] <- mean(results[,9])
}
#store main outputs as csv
write.csv(merged_result,file = "results.csv")

toc()


