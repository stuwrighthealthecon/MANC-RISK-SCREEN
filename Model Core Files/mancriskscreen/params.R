# Global option for whether to use corrected age-at-detection distribution
CORRECT_BC_AGE <- TRUE

#Screening uptake
uptakefirstscreen <- 0.625 #Uptake for the first screen
uptakeotherscreen <- 0.887 #Uptake if woman has attended >=1 screen
uptakenoscreen <- 0.208 #Uptake if woman has not previously attended screening

#Uptake for risk stratification
risk_uptake <- 1 #Uptake for risk prediction
risk_feedback <- 1 #Uptake for attending risk feedback consultation
screen_change <- 1 #Uptake for changing screening intervals based on risk

#Age of an individual at start of simulation
start_age <- 38 #Default is 38 to give time for tumours to develop pre-screening

#Set time horizon
time_horizon <- 100

#Set health and cost discount rates
discount_health <- 0.035
discount_cost <- 0.035

##############Set clinical parameters for screening########################

#Set parameters of a Weibull survival curve to represent all cause mortality
acmmortality_wb_a <- 9.38
#7.937
acmmortality_wb_b <- 89.35
#86.788

#Set parameters for all cause mortality following breast cancer
gamma_survival_1 <- exp(-5.618) #Exponential distribution scale parameter stage 1
gamma_survival_2 <- exp(-3.808) #Exponential distribution scale parameter stage 2
gamma_survival_3 <- exp(-2.731) #Exponential distribution scale parameter stage 3
gamma_stage <- c(gamma_survival_1, gamma_survival_2, gamma_survival_3)

# Expected population prevalence for doing Bayes correction
expected_prev <- .12

#Read in distribution of cancer incidence by age
Incidence_Mortality <- if (CORRECT_BC_AGE) {
  read.csv("Data/Incidence_Mortality_ONS2.csv") %>%
    mutate(
      BC_age = expected_prev *
        Cond.on.getting.BC..prob.of.getting.cancer.at.age.t /
        (surv / surv[1])
    ) %>%
    mutate(BC_age = replace(BC_age, age == 100, 1 - sum(BC_age[1:100])))
} else {
  read.csv("Data/Incidence_Mortality_ONS2.csv") %>%
    mutate(BC_age = Cond.on.getting.BC..prob.of.getting.cancer.at.age.t)
}

#Set metastatic cancer probabilities by cancer size
metastatic_prob <- data.frame(
  c(25, 35, 45, 55, 65, 75, 85),
  c(
    0.046218154,
    0.086659039,
    0.109768116,
    0.127099924,
    0.142505975,
    0.159837783,
    1.73E-01
  )
)
write.csv(metastatic_prob, "metaprob.csv")

#Create matrix of probability of cancer stage by cancer size
stage_by_size_mat <- data.frame(
  "v1" = c(0.383, 0.567, 0.611, 0.557, 0, 0),
  "v2" = c(0.033, 0.111, 0.180, 0.208, 0.723, 0.563),
  "v3" = c(0.058, 0.057, 0.089, 0.147, 0.206, 0.351),
  "v5" = c(0.525, 0.265, 0.120, 0.088, 0.071, 0.086)
)

#Set mean and sd of tumour doublings at clinical detection
clin_detection_m <- 6.5
clin_detection_sd <- 0.535

#Mammography with sensitivity conditional on tumour diameter parameters
beta1 <- 1.47
beta2 <- 6.51

#Mammography sensitivity by volpara density grade
VDG_interval <- c(4.5, 7.5, 15.5)
Sen_VDG <- c(0.75, 0.735, 0.598, 0.513)
Sen_VDG_av <- 0.666

#Supplemental screening sensitivity parameters
Mammo_cdr <- 4.2 #Cancer detection rate per 1000 high dense screens
MRI_cdr <- 5 #CDR for MRI in Mammogram negative women (incremental)
US_cdr <- 3 #CDR for US in Mammogram negative women (incremental)

density_cutoff <- 3

#Set tumour growth rate parameters
log_norm_mean <- 1.07
log_norm_sd <- 1.31
max_size <- 128 #Maximum size of tumours, mm diameter
start_size <- 0.25 #Starting size of tumours, diameter in mm
Vc = (4 / 3) * pi * (start_size / 2)^3 #Volume at start
Vm = (4 / 3) * pi * (max_size / 2)^3 #Max volume

#Metatstatic survival parameters
meta_survival_54 <- exp(-1.787) #Age <= 54
meta_survival_74 <- exp(-1.388) #Age 55-74
meta_survival_99 <- exp(-1.011) #Age 75+
metastatic_survival <- c(meta_survival_54, meta_survival_74, meta_survival_99)

#Set screening ages
screen_startage <- 50 #Starting age
screen_endage <- 70 #Ending age
high_risk_screentimes <- seq(screen_startage, screen_endage, 1) #Yearly screening
med_risk_screentimes <- seq(screen_startage, screen_endage, 2) #Two yearly screening
low_risk_screentimes <- seq(screen_startage, screen_endage, 3) #Three yearly screening

#Set maximum screening sensitivity
sensitivity_max <- 0.95

#Risk cut-offs for different screening approaches
risk_cutoffs_procas <- c(1.5, 3.5, 5, 8, 100) #PROCAS and PROCAS Full
risk_cutoffs_tert <- c(1.946527, 2.942792) #Tertiles of risk
low_risk_cut <- 1.5 #Cut off in low risk only strategies

#Cancer size cut-points
ca_size_cut <- c(0.025, 5, 10, 15, 20, 30, 128) #Size cut-points for deciding stage

##########False Positive and Overdiagnosis parameters################
recall_rate <- 0.025 #UK false-positive rate
biopsy_rate <- 0.506 #Proportion of referrals without cancer that have biopsy

#### Drug data ####

# First bring in log hazard ratios from networked analysis
loghaz_ests <- readRDS("Data/PreventionOutputs.RDS")
efficacy_ests <- loghaz_ests[1]
dropout_ests <- loghaz_ests[4]

# Extract parameters for multivariate normal draws
efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

# Fit a Weibull distribution to data in Incidence_Mortality. Drug acts to change scale parameter.
weibull_fit <- sample(
  Incidence_Mortality$age,
  size = 100000,
  prob = Incidence_Mortality$BC_age,
  replace = TRUE
) %>%
  fitdistr("weibull", lower = c(0, 0))
inc_scale <- weibull_fit$estimate["scale"]
inc_shape <- weibull_fit$estimate["shape"]

ana_eff <- exp(efficacy_ests$AnyBC$means$Anastrozole)
tam_eff <- exp(efficacy_ests$AnyBC$means$Tamoxifen)

full_course_len <- 5

# Assume constant drop out rate with 77% of individuals reaching 5yr mark
tam_completion_prob <- .77

# Estimate Anastrozole completion probability from Tamoxifen estimate and hazard ratios:
ana_completion_prob <- exp(
  dropout_ests$Adherence$means$Anastrozole -
    dropout_ests$Adherence$means$Tamoxifen
) *
  tam_completion_prob

completion_prob <- c(ana_completion_prob, tam_completion_prob)

# Now estimate per-unit-time dropout rates based on exponential time to dropout:
tam_dropout_rate <- (1. / full_course_len) * log(1 / (tam_completion_prob))
ana_dropout_rate <- (1. / full_course_len) * log(1 / (ana_completion_prob))

# Work out mean time taking each drug
mean_tam_length <- full_course_len *
  (1 - tam_completion_prob) /
  log(1 / tam_completion_prob)
mean_ana_length <- full_course_len *
  (1 - ana_completion_prob) /
  log(1 / ana_completion_prob)

# Quick check: quantities below give efficacy of taking full five-year course
# assuming linear relationship between time taking and hazard ratio. If both of
# these are less than one then the linear model is safe to use because no one
# goes past this point.
tam_full_course_eff <- 1 -
  (1 - tam_eff) * log(1 / tam_completion_prob) / (1 - tam_completion_prob)
ana_full_course_eff <- 1 -
  (1 - ana_eff) * log(1 / ana_completion_prob) / (1 - ana_completion_prob)
if ((tam_full_course_eff > 1) | (ana_full_course_eff > 1)) {
  print(
    "Assumption of linear change to hazard ratio with time taking drug will
        not work for one or both drugs being simulated."
  )
}

course_length <- c(5., 5.)

#Assign women to risk groups based on 10yr risk if using risk-stratified approach
if (screen_strategy == 1 | screen_strategy == 9) {
  risk_red <- matrix(
    c(
      ana_eff,
      tam_eff,
      ana_eff,
      tam_eff,
      ana_eff,
      tam_eff,
      ana_eff,
      tam_eff,
      ana_eff,
      tam_eff
    ),
    nrow = 5,
    ncol = 2
  )

  course_length <- c(5., 5.)

  uptake <- rbind(c(0., 0.), c(0., 0.), c(0., 0.), c(.71, .71), c(.71, .71))

  persistence <- matrix(
    c(
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate
    ),
    nrow = 5,
    ncol = 2
  )
} else if (screen_strategy == 2) {
  risk_red <- matrix(
    c(ana_eff, tam_eff, ana_eff, tam_eff, ana_eff, tam_eff),
    nrow = 3,
    ncol = 2
  )

  uptake <- rbind(c(0., 0.), c(0., 0.), c(.71, .71))

  persistence <- matrix(
    c(
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate,
      ana_dropout_rate,
      tam_dropout_rate
    ),
    nrow = 3,
    ncol = 2
  )
} else if (screen_strategy == 7 | screen_strategy == 8) {
  risk_red <- matrix(c(ana_eff, tam_eff, ana_eff, tam_eff), nrow = 2, ncol = 2)

  uptake <- rbind(c(0., 0.), c(.71, .71))

  persistence <- matrix(
    c(ana_dropout_rate, tam_dropout_rate, ana_dropout_rate, tam_dropout_rate),
    nrow = 2,
    ncol = 2
  )
} else {
  risk_red <- matrix(c(ana_eff, tam_eff), nrow = 1, ncol = 2)

  uptake <- matrix(c(0., 0.), nrow = 1, ncol = 2)

  persistence <- matrix(
    c(ana_dropout_rate, tam_dropout_rate),
    nrow = 1,
    ncol = 2
  )
}

age_prescribed <- 50

median_age_at_men <- 51
prob_premen <- exp(-age_prescribed * (log(2) / median_age_at_men)) # Assuming exponential distribution (not correct!)

cost_in_full_courses <- TRUE # If TRUE then each person to be prescribed drug incurs cost of full course, otherwise cost is proportional to time taking

#######################Cost Data#########################################

cost_strat <- 7.40 #Cost of risk prediction
cost_screen_base <- 39.93 #Cost of Mammography
cost_follow_up_base <- 130.97 #Cost of follow-up
cost_biop_base <- 558 #Cost of biopsyy
cost_DCIS_base <- 12346.83 #Cost of treating DCIS
cost_US_base <- 102.00 #Cost of ultrasound
cost_MRI_base <- 162.00 #Cost of MRI
cost_drug_base <- c(100., 100.) # Cost of full course of drug

#If deterministic analysis then set costs as base costs
if (PSA == 0) {
  cost_DCIS <- cost_DCIS_base
  cost_screen <- cost_screen_base
  cost_follow_up <- cost_follow_up_base
  cost_biop <- cost_biop_base
  cost_US <- cost_US_base
  cost_MRI <- cost_MRI_base
  cost_drug <- cost_drug_base
}

#Set up look-up table for treatment costs
tbl <- tribble(
  ~Yr , ~Early_18.64 , ~Late_18.64 , ~Diff1 , ~Early_65plus , ~Late_65plus , ~Diff2 ,
    0 ,          464 ,         607 ,    143 ,          1086 ,         1324 ,    238 ,
    1 ,        10746 ,       13315 ,   2569 ,          7597 ,         8804 ,   1207 ,
    2 ,         3357 ,        5785 ,   2429 ,          2529 ,         3650 ,   1121 ,
    3 ,         1953 ,        3782 ,   1829 ,          2156 ,         3170 ,   1014 ,
    4 ,         1627 ,        2932 ,   1305 ,          2230 ,         2924 ,    693 ,
    5 ,         1617 ,        2841 ,   1225 ,          2077 ,         2957 ,    880 ,
    6 ,         1547 ,        2645 ,   1099 ,          2174 ,         2783 ,    609 ,
    7 ,         1394 ,        2618 ,   1225 ,          2063 ,         2903 ,    840 ,
    8 ,         1376 ,        2559 ,   1183 ,          2134 ,         2454 ,    320 ,
    9 ,         1279 ,        1848 ,    569 ,          2204 ,         2932 ,    728
) %>%
  dplyr::select(-Diff1, -Diff2) %>%
  pivot_longer(
    cols = contains("6"),
    names_to = c("Stage", "Age"),
    names_sep = "_",
    values_to = "Cost"
  ) %>%
  group_by(Stage, Age) %>%
  mutate(
    DCost = Cost - first(Cost),
    DCost.i = DCost * 1.402093, # NHSCII inflator for 2010/11-->2023/2024
    disc = 1 / 1.035^(Yr - 0.5),
    DCost.i.d = DCost.i * disc,
    CDCost.i.d = cumsum(DCost.i.d),
    Yr1 = as.factor(Yr == 1),
    Yr2 = as.factor(Yr == 2),
    Yr3 = as.factor(Yr == 3)
  ) %>%
  filter(Yr > 0) %>%
  arrange(Stage, Age, Yr)

#Create predictive model of cost by age of diagnosis, stage, life expectancy
mod <- lm(
  data = tbl,
  formula = log(DCost) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age
)

#Prediction matrix
tblNewDat <- crossing(
  Yr = 1:50,
  Stage = c("Early", "Late"),
  Age = c("18.64", "65plus")
) %>%
  mutate(
    Yr1 = as.factor(Yr == 1),
    Yr2 = as.factor(Yr == 2),
    Yr3 = as.factor(Yr == 3)
  )

#Generate cost predictions
tblNewDat %>%
  bind_cols(pred = mod %>% predict(newdata = tblNewDat)) %>%
  mutate(DCost.p = exp(pred)) -> tblPred

#Make lookup table for costs
tblLookup <- tblPred %>%
  filter(Yr == 1) %>%
  mutate(across(c(Yr, pred, DCost.p), ~0)) %>%
  bind_rows(tblPred) %>%
  group_by(Stage, Age) %>%
  mutate(
    DCost.p.i = DCost.p * 1.402093, # NHSCII inflator for 2010/11-->2023/2024
    disc = 1 / 1.035^(Yr - 0.5),
    DCost.p.i.d = DCost.p.i * disc,
    CDCost.p.i.d = cumsum(DCost.p.i.d),
    StageEarly = Stage == "Early",
    AgeYoung = Age == "18.64"
  ) %>%
  arrange(Stage, Age, Yr) %>%
  ungroup()

#######################Utility Weights#########################################

#Set age adjusted utility values
utility_ages <- read.csv("Data/age_related_utilities.csv")

#Set time independent utility decrements
#Utility of DCIS
utility_DCIS <- 1 #assumes no effect

#Set first year cancer utilities:
utility_stage_cat_y1 <- c(
  "stage1" = 0.82 / 0.822,
  "stage2" = 0.82 / 0.822,
  "stage3" = 0.75 / 0.822,
  "Metastatic" = 0.75 / 0.822,
  "DCIS" = utility_DCIS
)

#Set following year cancer utilities:
utility_stage_cat_follow <- c(
  "stage1" = 0.82 / 0.822,
  "stage2" = 0.82 / 0.822,
  "stage3" = 0.75 / 0.822,
  "Metastatic" = 0.75 / 0.822,
  "DCIS" = utility_DCIS
)
