#########################Lookup function for treatment costs############################

fnLookupBase <- function(iStage, iAge, iLE) {
  as.numeric(tblLookup$CDCost.p.i.d[
    tblLookup$Stage == iStage & tblLookup$Age == iAge & tblLookup$Yr == iLE
  ])
}

########################stage calculator#######################################

stage_by_size <- function(Ca_size) {
  stage_cat <- 0

  #First determine if advanced cancer or not based on metastatic prob by size (categorical)
  if (Ca_size <= 25) {
    m_size <- 25
  } else {
    m_size <- ceiling((Ca_size - 25) / 10) * 10 + 25
  }
  if (m_size > 85) {
    m_size <- 85
  }
  if (
    dqrunif(1, 0, 1) <
      metastatic_prob[match(m_size, metastatic_prob[, 1]), 2] &&
      stage_cat == 0
  ) {
    stage_cat <- 4
  }

  #Sample from stage 1,2 & 3 with probability of each based on the correct row of stagebysize matrix
  #Ca_size is continuous, need to match to closest larger value in stage_by_size column 1
  #NB stage 5=DCIS
  if (stage_cat == 0) {
    size_cat <- findInterval(Ca_size, ca_size_cut)
    stage_cat <- sample(
      x = c(1, 2, 3, 5),
      size = 1,
      prob = c(stage_by_size_mat[size_cat, ])
    ) #1 best 3 worst prognosis
  }
  #Return the stage category
  result <- stage_cat

  return(result)
}
cmp_stage_by_size <- cmpfun(stage_by_size)

#######################Screening test results simulation##########################

#Inputs are tumour diameter, VDG, MRI_screening(0/1), US_screening(0/1)
screening_result <- function(Ca_size, VDG, MRI_screening, US_screening) {
  #Calculate size specific sensitivity
  Sensitivity <- if (
    exp((Ca_size - beta2) / beta1) / (1 + exp((Ca_size - beta2) / beta1)) >
      sensitivity_max
  ) {
    sensitivity_max
  } else {
    exp((Ca_size - beta2) / beta1) / (1 + exp((Ca_size - beta2) / beta1))
  } #use to set max sensitivity 0.95

  #Adjust sensitivity for breast density
  dense_OR <- (Sen_VDG[VDG] / (1 - Sen_VDG[VDG])) /
    (Sen_VDG_av / (1 - Sen_VDG_av))
  Sensitivity <- ((Sensitivity / (1 - Sensitivity)) * dense_OR) /
    (1 + ((Sensitivity / (1 - Sensitivity)) * dense_OR))

  #Draw random number used to compare to Sensitivity with and without supplemental screening
  rnd_1 <- dqrunif(1, 0, 1)

  #Is cancer detected by mammogram?
  if (rnd_1 < Sensitivity) {
    Screen_detected_ca <- 1
    Mammo_detected_ca <- 1 #Keep track of which stage it is detected
  } else {
    Screen_detected_ca <- 0
    Mammo_detected_ca <- 0
  }

  #Is cancer detected by supplemental tests?
  if (Screen_detected_ca == 0) {
    if (MRI_screening == 1) {
      MRI_supp_odds <- (Sensitivity / (1 - Sensitivity)) *
        ((MRI_cdr + Mammo_cdr) / Mammo_cdr)
      MRI_supp_sens <- MRI_supp_odds / (MRI_supp_odds + 1)
      if (rnd_1 < MRI_supp_sens) {
        Screen_detected_ca <- 1
        MRI_detected_ca <- 1
      } else {
        MRI_detected_ca <- 0
      }
    } else {
      MRI_detected_ca <- 0
    }

    if (US_screening == 1) {
      US_supp_odds <- (Sensitivity / (1 - Sensitivity)) *
        ((US_cdr + Mammo_cdr) / Mammo_cdr)
      US_supp_sens <- US_supp_odds / (US_supp_odds + 1)
      if (rnd_1 < US_supp_sens) {
        Screen_detected_ca <- 1
        US_detected_ca <- 1
      } else {
        US_detected_ca <- 0
      }
    } else {
      US_detected_ca <- 0
    }
  } else {
    US_detected_ca <- 0
    MRI_detected_ca <- 0
  }

  #Uses estimate based on increased cancer detection rate

  result <- c(
    Screen_detected_ca,
    Mammo_detected_ca,
    MRI_detected_ca,
    US_detected_ca
  )
  return(result)
}
cmp_screening_result <- cmpfun(screening_result)

############################Simulate survival by stage##########################

Ca_survival_time <- function(stage_cat, Mort_age, age, ca_incidence_age) {
  #Assign survival for non-metastatic cancer
  if (stage_cat < 4) {
    survival_time <- -(log(x = dqrunif(1, 0, 1)) / gamma_stage[stage_cat]) #inverse of cdf when rate is gamma_stage[x]

    #Adjust for additional mortality at ages above 65
    if (ca_incidence_age > 65) {
      survival_time <- -(log(x = dqrunif(1, 0, 1)) /
        ((Incidence_Mortality$X10year.mort.prob[min(
          (floor(ca_incidence_age) + 1),
          100
        )] /
          Incidence_Mortality$X10year.mort.prob[66]) *
          gamma_stage[stage_cat]))
    }

    #Data are for 10-year survival, after 10 years assume that pop mortality rates apply
    if (survival_time > 10) {
      Mort_age <- qweibull(
        p = dqrunif(
          n = 1,
          min = pweibull(
            q = ca_incidence_age + 10,
            shape = acmmortality_wb_a,
            scale = acmmortality_wb_b
          ),
          max = 1
        ),
        shape = acmmortality_wb_a,
        scale = acmmortality_wb_b
      )
      if (Mort_age > time_horizon) {
        Mort_age <- time_horizon
      }
      survival_time <- Mort_age - age
    }
  }

  #Assign survival for metastatic cancer
  if (stage_cat == 4) {
    if (age < 55) {
      age_cat_M <- 1
    }
    if (age >= 55 && age < 75) {
      age_cat_M <- 2
    }
    if (age >= 75) {
      age_cat_M <- 3
    }
    survival_time <- -(log(dqrunif(1, 0, 1)) / metastatic_survival[age_cat_M])
  }

  #Assign survival for DCIS i.e. no effect
  if (stage_cat == 5) {
    survival_time <- (Mort_age - ca_incidence_age)
  }

  if (ca_incidence_age + survival_time > time_horizon) {
    survival_time <- time_horizon - ca_incidence_age
  }
  result <- ca_incidence_age + survival_time
  #Reduce age of death if cancer causes woman to die earlier
  if (result < Mort_age) {
    Mort_age <- result
  }
  return(result)
}
cmp_ca_survival_time <- cmpfun(Ca_survival_time)

###################################QALY Counter##########################
QALY_counter <- function(Mort_age, incidence_age_record, stage_cat) {
  #QALY counter
  #Set up a QALY vector of length equal to life years
  QALY_length <- ceiling(Mort_age) - (screen_startage - 1)

  #If less than 1 life year lived, set length to 1
  if (QALY_length < 1) {
    QALY_length <- 1
  }

  #Fill QALY vector with 0's
  QALY_vect <- rep(0, QALY_length)

  #Fill QALY vector with discounted age related utility values
  for (y in 1:length(QALY_vect)) {
    QALY_vect[y] <- (utility_ages[
      match((ceiling((screen_startage - 1) + y)), utility_ages[, 1]),
      2
    ]) *
      (1 / (1 + discount_health)^(y))
    QALY_vect[QALY_length] <- QALY_vect[QALY_length] *
      (1 - (ceiling(Mort_age) - Mort_age))
  }

  #If cancer occurs then fill QALY vector with discounted cancer utilities from incidence age
  #NB this code accounts for partial years spent in different health states
  if (incidence_age_record > 0) {
    QALY_vect[
      floor(incidence_age_record) - screen_startage
    ] <- utility_stage_cat_y1[stage_cat] *
      QALY_vect[floor(incidence_age_record) - screen_startage] *
      (1 - (incidence_age_record - floor(incidence_age_record)))
  }
  if (incidence_age_record > 0 & Mort_age - incidence_age_record > 1) {
    QALY_vect[
      (floor(incidence_age_record) - screen_startage) + 1
    ] <- (utility_stage_cat_y1[
      stage_cat
    ] *
      QALY_vect[(floor(incidence_age_record) - screen_startage) + 1] *
      (incidence_age_record - floor(incidence_age_record))) +
      (utility_stage_cat_follow[stage_cat] *
        QALY_vect[(floor(incidence_age_record) - screen_startage) + 1] *
        (1 - (incidence_age_record - floor(incidence_age_record))))
  }
  if (
    incidence_age_record > 0 &&
      ceiling(
        if (Mort_age < 100) {
          Mort_age
        } else {
          100
        }
      ) >
        incidence_age_record + 2
  ) {
    for (y in (incidence_age_record + 2):min(
      (incidence_age_record + 8),
      ceiling(
        if (Mort_age < 100) {
          Mort_age
        } else {
          100
        }
      )
    )) {
      QALY_vect[y - screen_startage] <- QALY_vect[y - screen_startage] *
        utility_stage_cat_follow[stage_cat]
    }
  }
  return(QALY_vect)
}
cmp_QALY_counter <- cmpfun(QALY_counter)
