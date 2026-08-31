controls <- list(
  "strategies" = c(3), #A vector of strategies to evaluate
  "gensample" = TRUE, #Whether to generate a new sample to simulate
  "MISCLASS" = TRUE, #whether to include risk misclassification in analysis
  "PREVENTATIVE_DRUG" = FALSE, #whether to include chemoprevention in analysis
  "supplemental_screening" = FALSE, #whether supplemental screening is used for women with dense breasts
  "PSA" = FALSE, #whether to conduct a probabilistic sensitivity analysis
  "intervals" = FALSE, #whether to conduct a PSA with wide intervals for GAM estimations
  "desired_cases" = 30000, #apprximate number of cancer cases required in simulation
  "mcruns" = 1 #number of monte carlo runs in PSA/intervals
) #set number of cores for parallel processing

if (!require("pacman", quietly = TRUE)) install.packages("pacman")

#Install/load required packages
pacman::p_load(
  doParallel,
  MASS,
  dqrng,
  compiler,
  tidyverse,
  iterators,
  here,
  purr,
  tictoc
)

MISCLASS <- controls$MISCLASS # Set to TRUE to include impact of errors in risk prediction in model
PREVENTATIVE_DRUG <- controls$PREVENTATIVE_DRUG # Set to TRUE to simulate preventative drugs

# Add specifiers for output files
det_output_path <- "Deterministic results/"
psa_output_path <- "PSA results/"
if (MISCLASS & PREVENTATIVE_DRUG) {
  dir.create(
    "Deterministic results/misclassification_and_preventative_drug",
    showWarnings = FALSE
  )
  det_output_path <- "Deterministic results/misclassification_and_preventative_drug/"
  dir.create(
    "PSA results/misclassification_and_preventative_drug",
    showWarnings = FALSE
  )
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
} else {
  if (MISCLASS) {
    dir.create("Deterministic results/misclassification", showWarnings = FALSE)
    det_output_path <- "Deterministic results/misclassification/"
    dir.create("PSA results/misclassification", showWarnings = FALSE)
    psa_output_path <- "PSA results/misclassification/"
  }
  if (PREVENTATIVE_DRUG) {
    dir.create("Deterministic results/preventative_drug", showWarnings = FALSE)
    det_output_path <- "Deterministic results/preventative_drug/"
    dir.create(
      "PSA results/misclassification_and_preventative_drug",
      showWarnings = FALSE
    )
    psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
  }
}

sample_fname <- "possample"

#####Choose screening programme and related parameters##########
tic()

#Load file containing required functions for the model
source(file = "MANC_RISK_SCREEN_functions.R")
source(file = "risksample function.R")
source(file = "negsample function.R")
source(file = "des_vectorised_df.R")

#Set the screening strategy: 1=PROCAS, 2=Risk tertiles, 3=3 yearly, 4=2 yearly,
#5=5 yearly, 6=2 rounds at 50 and 60 (10 yearly), 7=Low risk (5 yearly),
#8=Low risk (6 yearly),#9=Fully stratified screening programmes
#Other num=no screening
screen_strategies <- unlist(controls$strategies)
for (r in 1:length(screen_strategies)) {
  screen_strategy <- screen_strategies[r]

  #Turn supplemental Screening (MRI and US) on (1) or off (0)
  supplemental_screening <- ifelse(
    controls$supplemental_screening == TRUE,
    1,
    0
  )

  #Generate new sample? 1=YES, any other number NO
  gensample <- ifelse((controls$gensample == TRUE & r == 1), 1, 0)

  #Deterministic (0) or Probabilistic Analysis (1)
  PSA = ifelse(controls$PSA == TRUE, 1, 0)

  #Standard (0) or wide (1) distributions for PSA
  #Wide intervals recommended for generating data to predict GAM model
  intervals = ifelse(controls$intervals == TRUE, 1, 0)

  #Set loop numbers
  expected_prev <- .12
  desired_cases <- controls$desired_cases
  inum <- ceiling((desired_cases / expected_prev)) #Individual women to be sampled to give desired number of positive cancer cases
  mcruns <- controls$mcruns #Monte Carlo runs used if PSA switched on
  seed <- set.seed(controls$seed) #Set seed for random draws

  #################################Import baseline parameters####################

  source(file = "params.R")

  #########################CREATE SAMPLE OF WOMEN FOR MODEL###################
  if (MISCLASS) {
    if (gensample == 1) {
      dir.create("Risksamplewithmisclass", showWarnings = FALSE)
      create_sample_with_misclass(PSA, intervals, seed, screen_strategy)
    }
  } else {
    if (gensample == 1) {
      dir.create("Risksample", showWarnings = FALSE)
      cmp_create_sample(PSA, intervals, seed, screen_strategy)
    }
  }
  
  ################Outer Individual sampling loop##############################

  start_time <- Sys.time()
  
  if (MISCLASS) {
    load(paste(
      "Risksamplewithmisclass/",
      sample_fname,
      ".Rdata",
      sep = ""
    ))
  } else {
    load(paste("Risksample/", sample_fname, ".Rdata", sep = ""))
  }
  prefix <- paste("^", "X", 1, ".", sep = "")
  names(risksample) <- sub(prefix, "", names(risksample))
  
  #Assign cancer size at diagnosis
  age_mat <- outer(
    risksample$life_expectancy,
    Incidence_Mortality$age,
    FUN = ">"
  ) * Incidence_Mortality$BC_age
  
  # Normalise rows
  age_mat <- age_mat / rowSums(age_mat)
  
  # Vectorised inverse CDF sampling
  cum_mat     <- t(apply(age_mat, 1, cumsum))
  u           <- dqrunif(nrow(risksample), 0, 1)
  col_indices <- rowSums(cum_mat < u) + 1L
  risksample$ca_incidence <- Incidence_Mortality$age[col_indices]
  
  #Add error for month of diagnosis
  risksample$ca_incidence<-risksample$ca_incidence+dqrunif(nrow(risksample), 0, 1)
  
  #Determine size of cancer at diagnosis
  risksample$clin_detect_size_g <- start_size * 2^risksample$clinical_detect_size
  
  #Calculate tumour genesis age
  t_gen <- ((log((Vm / Vc)^0.25 - 1) -
               log((Vm / ((4 / 3) * pi * (risksample$clin_detect_size_g / 2)^3))^0.25 - 1)) /
              (0.25 * risksample$growth_rate)) #Calculate time to get to clinical detection size
  risksample$genage <- risksample$ca_incidence - t_gen
  
  if (MISCLASS) {
    #Assign women to risk groups based on 10yr risk if using risk-stratified approach
    if (screen_strategy == 1 | screen_strategy == 9) {
      risksample$risk_group <- 1 +
        findInterval(risksample$tenyrrisk_est, risk_cutoffs_procas)
    } else if (screen_strategy == 2) {
      risksample$risk_group <- 1 +
        findInterval(risksample$tenyrrisk_est, risk_cutoffs_tert)
    } else if (screen_strategy == 7 | screen_strategy == 8) {
      risksample$risk_group <- ifelse(
        risksample$tenyrrisk_est < low_risk_cut,
        1,
        2
      )
    }
  } else {
    if (screen_strategy == 1 | screen_strategy == 9) {
      risksample$risk_group <- 1 +
        findInterval(risksample$tenyrrisk, risk_cutoffs_procas)
    } else if (screen_strategy == 2) {
      risksample$risk_group <- 1 +
        findInterval(risksample$tenyrrisk, risk_cutoffs_tert)
    } else if (screen_strategy == 7 | screen_strategy == 8) {
      risksample$risk_group <- ifelse(
        risksample$tenyrrisk < low_risk_cut,
        1,
        2
      )
    }
  }
  
  risk_groups<-unique(risksample$risk_group)
  
  if (PREVENTATIVE_DRUG) {
    # # Add extra fields for drug:
    nsample <- nrow(risksample)
    # Use 1, 2 coding for menopause status to match indexing for drug efficacy/uptake
    risksample$starting_menses_status <- ifelse(
      dqrunif(nsample, 0, 1) < prob_premen,
      1,
      2
    )
    risksample$takes_drug <- logical(nsample)
    risksample$time_taking_drug <- numeric(nsample)
    
    if (PSA == 1) {
      # Redraw effects for PSA, assuming Monte Carlo draws of parameters are the same for each individual
      drug_matrix_list <- redraw_drug_pars(risksample[1, ])
      risk_red <- drug_matrix_list[[1]]
      uptake <- drug_matrix_list[[2]]
      persistence <- drug_matrix_list[[3]]
    }
  }
  
  #Assign women to supplemental screening if switched on and criteria met
  if (supplemental_screening == 1) {
    for (i in 1:length(risksample$MRI_screen)) {
      if (
        risksample[i, "VDG"] >= density_cutoff &
        risksample[
          i,
          ifelse(MISCLASS == 1, "tenyrrisk_true", "tenyrrisk_est")
        ] >=
        8
      ) {
        risksample[i, "MRI_screen"] < 1
      } else if (
        risksample[i, "VDG"] >= density_cutoff &
        risksample[
          i,
          ifelse(MISCLASS == 1, "tenyrrisk_true", "tenyrrisk_est")
        ] <
        8
      ) {
        risksample[i, "US_screen"] <- 1
      }
    }
  }

  
  if(MISCLASS){splitmaster<-subset(risksample,select=-c(VBD,cancer,feedback,liferisk_est,liferisk_true,tenyrrisk_est,tenyrrisk_true))
  }else{splitmaster<-subset(risksample,select=-c(VBD,cancer,feedback,liferisk,tenyrrisk))}
  
  #Set loop to divide i loop into a number of sub-loops in case of simulation break
  for (ii in 1:length(risk_groups)) {

    risksample<-splitmaster %>% 
      filter(risk_group==risk_groups[ii])
    
    
    if (PREVENTATIVE_DRUG & risksample$risk_group[1] != 0){
      
      drug_IM <- Incidence_Mortality
      
      idx <- cbind(risksample$risk_group, risksample$starting_menses_status)
      uptake_probs <- uptake[idx]
      risksample$uptake<-dqrunif(nrow(risksample), 0, 1) < uptake_probs
      
      risksample$time_taking_drug <- risksample$uptake*pmin(
        rexp(nrow(risksample), rate = persistence[idx]),
        course_length
      )
      
      risksample$weibullrisk <- (1 -
                                   risk_red[idx]) *
        risksample$time_taking_drug *
        log(1 / completion_prob[risksample$starting_menses_status])
      risksample$weibullrisk <- inc_scale / risksample$weibullrisk
      
      ages <- drug_IM$age[start_age:101]
      n <- nrow(risksample)
      
      # Matrix of weibull densities: rows = ages, cols = individuals
      prob_matrix <- outer(
        ages,
        risksample$weibullrisk,
        FUN = function(age, scale) dweibull(age, shape = inc_shape, scale = scale)
      )
      
      # Transpose to n x ages for consistency with earlier vectorised code
      prob_matrix <- t(prob_matrix)
      
      # Zero out ages >= life_expectancy for each individual
      age_mask <- outer(risksample$life_expectancy, ages, FUN = ">")
      prob_matrix <- prob_matrix * age_mask
      
      # Normalise each row
      row_sums <- rowSums(prob_matrix)
      prob_matrix <- prob_matrix / row_sums
      
      # Inverse CDF sampling
      cum_probs <- t(apply(prob_matrix, 1, cumsum))
      u <- dqrunif(n, 0, 1)
      col_indices <- rowSums(cum_probs < u) + 1L
      col_indices <- pmin(col_indices, length(ages))
      
      risksample$ca_incidence <- ifelse(risksample$time_taking_drug>0, ages[col_indices] + dqrunif(n, 0, 1), risksample$ca_incidence) 
    }
    
    # #If cancer occurs after age of death, take the month jitter off of cancer occurance to get it in the right time
    risksample$ca_incidence<-ifelse(
      risksample$life_expectancy <= risksample$ca_incidence,
           floor(risksample$life_expectancy),
          risksample$ca_incidence)
    risksample$life_expectancy <- ifelse(
      risksample$life_expectancy >= time_horizon,
      99.99,
      risksample$life_expectancy
    )
    
    screen_times <- c(0)
    if (screen_strategy == 1) {
      if (risksample$risk_group[1] < 4) {
        screen_times <- low_risk_screentimes
      } else if (risksample$risk_group[1] > 3 & risksample$risk_group[1] < 5) {
        screen_times <- med_risk_screentimes
      } else if (risksample$risk_group[1] > 4) {
        screen_times <- high_risk_screentimes
      }
    }
    if (screen_strategy == 2) {
      if (risksample$risk_group[1] == 1) {
        screen_times <- low_risk_screentimes
      } else if (risksample$risk_group[1] == 2) {
        screen_times <- med_risk_screentimes
      } else if (risksample$risk_group[1] == 3) {
        screen_times <- high_risk_screentimes
      }
    }
    if (screen_strategy == 3) {
      screen_times <- low_risk_screentimes
    }
    if (screen_strategy == 4) {
      screen_times <- med_risk_screentimes
    }
    if (screen_strategy == 5) {
      screen_times <- seq(screen_startage, screen_startage + (5 * 4), 5)
    }
    if (screen_strategy == 6) {
      screen_times <- seq(screen_startage, screen_startage + 10, 10)
    }
    if (screen_strategy == 7) {
      if (risksample$risk_group[1] == 1) {
        screen_times <- seq(screen_startage, screen_startage + (5 * 4), 5)
      }
      if (risksample$risk_group[1] == 2) {
        screen_times <- low_risk_screentimes
      }
    } else if (screen_strategy == 7) {
      screen_times <- low_risk_screentimes
    }
    if (screen_strategy == 8) {
      if (risksample$risk_group[1] == 1) {
        screen_times <- seq(screen_startage, screen_startage + (6 * 3), 6)
      }
      if (risksample$risk_group[1] == 2) {
        screen_times <- low_risk_screentimes
      }
    } else if (screen_strategy == 8) {
      screen_times <- low_risk_screentimes
    }
    if (screen_strategy == 9) {
      if (risksample$risk_group[1] == 1) {
        screen_times <- seq(screen_startage, screen_startage + (5 * 4), 5)
      } else if (risksample$risk_group[1] == 2 | risksample$risk_group[1] == 3) {
        screen_times <- low_risk_screentimes
      } else if (risksample$risk_group[1] == 4) {
        screen_times <- med_risk_screentimes
      } else if (risksample$risk_group[1] == 5) {
        screen_times <- high_risk_screentimes
      }
    } else if (screen_strategy == 9) {
      screen_times <- low_risk_screentimes
    }
    
    if(length(screen_times)>1){
    #Add blank columns for potential screen times
    risksample[paste0("screen_",seq_len(length(screen_times)))]<-0
 
      #Draw attendance at first screen
      risksample[,"screen_1"] <- rbinom(
        nrow(risksample),
        1,
        uptakefirstscreen
      )
 
      #Loop through remaining screens conditional on previous attendance
      for (i in 1:(length(screen_times) - 1)) {
        risksample[, paste0("screen_", i+1)] <- ifelse(
          rowSums(risksample[, paste0("screen_", 1:i), drop = FALSE]) >= 1,
          rbinom(nrow(risksample), 1, uptakeotherscreen),
          rbinom(nrow(risksample), 1, uptakenoscreen)
        )
      }
    }else{risksample$screen_1<-999}
    
        #If PSA switched on, replace base case parameter values with Monte Carlo draws
        if (PSA == 1) {
          beta1 <- risk_data$PSA_beta_1
          beta2 <- risk_data$PSA_beta_2

          log_norm_mean <- risk_data$PSA_log_norm_mean
          log_norm_sd <- risk_data$PSA_log_norm_sd

          gamma_survival_1 <- exp(risk_data$PSA_gamma_survival_1)
          gamma_survival_2 <- exp(risk_data$PSA_gamma_survival_2)
          gamma_survival_3 <- exp(risk_data$PSA_gamma_survival_3)
          gamma_stage <- c(gamma_survival_1, gamma_survival_2, gamma_survival_3)

          meta_survival_54 <- exp(risk_data$PSA_meta_survival_54)
          meta_survival_74 <- exp(risk_data$PSA_meta_survival_74)
          meta_survival_99 <- exp(risk_data$PSA_meta_survival_99)
          metastatic_survival <- c(
            meta_survival_54,
            meta_survival_74,
            meta_survival_99
          )

          Sen_VDG <- c(
            risk_data$PSA_VDG1_sen,
            risk_data$PSA_VDG2_sen,
            risk_data$PSA_VDG3_sen,
            risk_data$PSA_VDG4_sen
          )
          Sen_VDG_av <- mean(Sen_VDG)

          MRI_cdr <- risk_data$PSA_MRI_cdr
          US_cdr <- risk_data$PSA_US_cdr

          risk_data$growth_rate <- qlnorm(
              dqrunif(1, 0, 1),
              meanlog = log_norm_mean,
              sdlog = sqrt(log_norm_sd)
            )

          utility_stage_cat_y1 <- c(
            "stage1" = risk_data$PSA_util_1to3 / 0.822,
            "stage2" = risk_data$PSA_util_1to3 / 0.822,
            "stage3" = risk_data$PSA_util_1to3 / 0.822,
            "Metastatic" = risk_data$PSA_util_4 / 0.822,
            "DCIS" = utility_DCIS
          )

          utility_stage_cat_follow <- c(
            "stage1" = risk_data$PSA_util_1to3 / 0.822,
            "stage2" = risk_data$PSA_util_1to3 / 0.822,
            "stage3" = risk_data$PSA_util_1to3 / 0.822,
            "Metastatic" = risk_data$PSA_util_4 / 0.822,
            "DCIS" = utility_DCIS
          )

          cost_strat <- risk_data$PSA_cost_strat
          cost_DCIS <- cost_DCIS_base * (1 + risk_data$PSA_costvar)
          cost_screen <- cost_screen_base * (1 + risk_data$PSA_costscreen)
          cost_follow_up <- cost_follow_up_base *
            (1 + risk_data$PSA_cost_follow_up)
          cost_biop <- cost_biop_base * (1 + risk_data$PSA_cost_biop)
          cost_US <- cost_US_base * (1 + risk_data$PSA_cost_US)
          cost_MRI <- cost_MRI_base * (1 + risk_data$PSA_cost_MRI)
          cost_drug <- cost_drug_base * (1 + risk_data$PSA_cost_drug)
        }

    run_des_vectorised(risksample, start_age=49,
                       PREVENTATIVE_DRUG = controls$PREVENTATIVE_DRUG,
                       discount_cost = discount_cost, screen_startage = 50,
                       cost_screen = cost_screen,cost_strat=cost_strat,cost_US=cost_US,
                       cost_MRI=cost_MRI,cost_follow_up=cost_follow_up,cost_DCIS=cost_DCIS,Vm=Vm,Vc=Vc,
                       screen_strategy = screen_strategy,PSA=PSA,recall_rate=recall_rate,
                       biopsy_rate = biopsy_rate,cost_biop=cost_biop)
  } #End i loop

  negsamplefn(screen_strategy, MISCLASS, PSA)

 print(paste("Strategy ",r," Complete")) 
}
toc()

