DO_INSTALL <- FALSE

if (DO_INSTALL){
  #Install required packages
  install.packages("doParallel")
  install.packages("MASS")
  install.packages("dqrng")
  install.packages("compiler")
  install.packages("tidyverse")
  install.packages("iterators")
}

MISCLASS <- TRUE # Set to TRUE to include impact of errors in risk prediction in model
PREVENTATIVE_DRUG <- TRUE # Set to TRUE to simulate preventative drugs

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
SEPARATE_SAMPLES <- TRUE # If true, only patients who develop cancer have their pathways simulated
if (SEPARATE_SAMPLES){
  sample_fname <- "possample_"
}else{
  sample_fname <- "risksample_"
}

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
screen_strategy<-2

#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Generate new sample? 1=YES, any other number NO
gensample<-1

#Deterministic (0) or Probabilistic Analysis (1)
PSA=0

#Standard (0) or wide (1) distributions for PSA
#Wide intervals recommended for generating data to predict GAM model
intervals=0

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Set loop numbers
chunks<-10 #Number of chunks to split inum into for faster running time
expected_prev <- .12
desired_cases <- 20000
inum <- ceiling((desired_cases / expected_prev)) #Individual women to be sampled to give desired number of positive cancer cases
inum <- chunks * ceiling(inum / chunks) # Make sure number of women is divisible by number of chunks
mcruns<-1 #Monte Carlo runs used if PSA switched on
seed<-set.seed(1) #Set seed for random draws

#Register number of cores for foreach loop
numcores<-19
registerDoParallel(cores=numcores)

#Load file containing required functions for the model
source(file="Functions/MANC_RISK_SCREEN_functions.R")
source(file="Functions/risksample function.R")

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

#### Drug data ####

# First bring in log hazard ratios from networked analysis
loghaz_ests <- readRDS("PreventionOutputs.RDS")
efficacy_ests <- loghaz_ests[1]
dropout_ests <- loghaz_ests[4]

# Extract parameters for multivariate normal draws
efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

# Fit a Weibull distribution to data in Incidence_Mortality. Drug acts to change scale parameter.
weibull_fit <- sample(Incidence_Mortality$age,
                      size=100000,
                      prob=Incidence_Mortality$Cond.on.getting.BC..prob.of.getting.cancer.at.age.t,
                      replace=TRUE) %>%
  fitdistr("weibull")
inc_scale <- weibull_fit$estimate["scale"]
inc_shape <- weibull_fit$estimate["shape"]

ana_eff <- exp(efficacy_ests$AnyBC$means$Anastrozole)
tam_eff <- exp(efficacy_ests$AnyBC$means$Tamoxifen)

full_course_len <- 5

# Assume constant drop out rate with 77% of individuals reaching 5yr mark
tam_completion_prob <- .77

# Estimate Anastrozole completion probability from Tamoxifen estimate and hazard ratios:
ana_completion_prob <- exp(dropout_ests$Adherence$means$Anastrozole - dropout_ests$Adherence$means$Tamoxifen) * tam_completion_prob

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


age_prescribed <- 50

median_age_at_men <- 51
prob_premen <- exp(- age_prescribed * (log(2) / median_age_at_men)) # Assuming exponential distribution (not correct!)

cost_in_full_courses <- TRUE # If TRUE then each person to be prescribed drug incurs cost of full course, otherwise cost is proportional to time taking

#######################Cost Data#########################################

cost_strat<-8.69 #Cost of risk prediction
cost_screen_base <- 62.21 #Cost of Mammography
cost_follow_up_base <- 109.4 #Cost of follow-up
cost_biop_base <- 297.98 #Cost of biopsyy
cost_DCIS_base <- 10107.80 #Cost of treating DCIS
cost_US_base <- 53.41 #Cost of ultrasound
cost_MRI_base <-117.10 #Cost of MRI
cost_drug_base <- c(100., 100.) # Cost of full course of drug

#If deterministic analysis then set costs as base costs
if(PSA==0){
  cost_DCIS<-cost_DCIS_base
  cost_screen<-cost_screen_base
  cost_follow_up <- cost_follow_up_base
  cost_biop <- cost_biop_base
  cost_US <- cost_US_base
  cost_MRI <-cost_MRI_base
  cost_drug <- cost_drug_base
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

      #Total life years
      LY_counter <- 0 #Total life years
      #Total QALYs
      QALY_counter <- 0 #Total QALYs
      
      #Get an all-cause mortality age and make sure this is greater than start
      #age and cancer incidence age
      Mort_age <- risk_data$life_expectancy
      
      #Other individual variables
      age <- start_age
      interval_ca <- 0 #Cancer clinically detected
      screen_detected_ca <- 0 #Cancer screen detected
      
      ##############################DES COMPONENT CANCER ###################################
        
        #Lifetime cancer incidence
        #Determines if a cancer occurs and at what age
        if (risk_data$cancer==1){
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
          if(ca_incidence_i[2] ==1){CD_age <- ca_incidence_i[1]} else{
            CD_age <- ca_incidence_i[1] + ((log((Vm/Vc)^0.25-1)-
                                              log((Vm/((4/3)*pi*(ca_incidence_i[4]/2)^3))^0.25-1))/(0.25*grow_rate_i)) - 
            ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(ca_incidence_i[3]/2)^3))^0.25-1))/(0.25*grow_rate_i))}
          cancer_diagnostic[8] <- c(CD_age)
          
          #Calculate tumour genesis age
          t_gen <- ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(CD_size/2)^3))^0.25-1))/(0.25*grow_rate_i)) #Calculate time to get to clinical detection size
          gen_age <- CD_age - t_gen
          
          #If cancer occurs after age of death, re-draw age of death
          if(Mort_age <= ca_incidence_age){Mort_age <-qweibull(
            p = dqrunif(n = 1,min = pweibull(
              q = CD_age,shape = acmmortality_wb_a,scale = acmmortality_wb_b),max = 1),
            shape = acmmortality_wb_a, scale = acmmortality_wb_b)}
          if(Mort_age >= time_horizon){Mort_age <- 99.99}
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
        current_discount<-(1/((1+discount_cost)^(Next_event_time+age-screen_startage)))
        
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
          
          #Add the cost of DCIS
          if(stage_cat == 5){
          costs = costs + (cost_DCIS*current_discount)}
          
          #Generate a cancer specific survival time, accounting for competing risks
          Ca_mort_age <- cmp_ca_survival_time(stage_cat,Mort_age,age,ca_incidence_age)
          
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
      
      #Update Life-year counter
      LY_counter <- Mort_age-start_age
      
      #Record total QALYs for J loop
      QALY_counter <- sum(cmp_QALY_counter(Mort_age,
                                           incidence_age_record,
                                           stage_cat),na.rm = TRUE)
      }else{
          ca_case <- 0
          
        #For non-cancer individuals
        screen_cost_vec<-rep(cost_screen,length(screen_times))
        follow_up_vec<-rbinom(length(screen_times),1,recall_rate)
        biop_vec<-follow_up_vec*(rbinom(length(follow_up_vec),1,biopsy_rate))
        screen_cost_vec<-screen_cost_vec+(follow_up_vec*cost_follow_up)+(biop_vec*cost_biop)
                                           
        if (screen_strategy==1 | screen_strategy==2 | screen_strategy==7 |
            screen_strategy==8 | screen_strategy==9){screen_cost_vec[1]<-screen_cost_vec[1]+cost_strat}
        discount_vec<-screen_times-rep(screen_startage,length(screen_times))
        for (i in length(screen_cost_vec)){
          screen_cost_vec[i]<-screen_cost_vec[i]*(1/((1+discount_cost)^discount_vec[i]))}
        
        #Screen counter
        screen_count<-length(screen_times)
        
        #Costs  
        costs<-max(0,sum(screen_cost_vec),na.rm=TRUE)
        
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
        
        #Update Life-year counter
        LY_counter <- Mort_age-start_age
        
        #Record total QALYs for J loop
        QALY_counter <- sum(cmp_QALY_counter(Mort_age,
                                             incidence_age_record,
                                             stage_cat),na.rm = TRUE)
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
               LY_counter-(screen_startage-start_age),
               cancer_diagnostic[2:3],
               Mort_age,cancer_diagnostic[10]))}else{
        #If PSA then record outputs + monte carlo draws
        return(as.numeric(c(QALY_counter,
                            costs,
                            screen_count,
                            cancer_diagnostic[8],
                            (screen_detected_ca+interval_ca),
                            screen_detected_ca,
                            screen_strategy,
                            risk_data$growth_rate,
                            LY_counter-(screen_startage-start_age),
                            c(risk_data[15:40]))))
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
    names(results)[10:35]<-c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
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
  write.csv(merged_result,file = paste(det_output_path,
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
    write.csv(merged_result,file = paste(psa_output_path,
                                         "PSAresults_strat_",
                                         screen_strategy,
                                         ".csv",
                                         sep=""))
  }

toc()

