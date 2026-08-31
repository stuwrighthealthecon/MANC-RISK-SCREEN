
# VECTORISED DES - Breast Cancer Screening Model (base R data.frame)

# Run_des_vectorised
run_des_vectorised <- function(risk_df,
                               screen_strategy,
                               start_age,
                               screen_startage,
                               discount_cost,
                               recall_rate,
                               Vm, Vc,
                               cost_screen, cost_strat, cost_US,
                               cost_MRI, cost_follow_up, cost_biop,
                               cost_DCIS, cost_drug,
                               biopsy_rate,
                               course_length,
                               PREVENTATIVE_DRUG,
                               PSA,
                               cost_in_full_courses,
                               uptake,
                               persistence) {
  
  # 0a. Initialise working data.frame and extract per-woman screen schedules
  df <- risk_df
  n  <- nrow(df)
  
  df$Mort_age <- df$life_expectancy
  df$CD_age   <- df$ca_incidence
  df$CD_size  <- df$clin_detect_size_g
  
  screen_cols <- grep("^screen_", names(df), value = TRUE)
  screen_cols <- screen_cols[order(as.integer(sub("screen_", "", screen_cols)))]
  attend_mat  <- as.matrix(df[, screen_cols])
  max_screens <- ncol(attend_mat)
  
  screen_mat                  <- t(t(attend_mat) * screen_times)
  screen_mat[screen_mat == 0] <- NA_real_
  
  # 0b. Initialise state and counter columns
  df$age                  <- start_age
  df$active               <- TRUE #Is person active in simulation (alive and cancer free)
  df$interval_ca          <- 0L #Interval cancer detected
  df$screen_detected_ca   <- 0L #Cancer found by screening
  df$screen_count         <- 0L #Number of screens (mammo) attended
  df$US_count             <- 0L #Number of ultrasounds attended
  df$MRI_count            <- 0L #Number of MRIs attended
  df$recall_count         <- 0L #Number of false-positive results
  df$lastscreen_count     <- 0L #Last screen attended
  df$sdfirst_cancer       <- 0L #Cancer found at first screen
  df$sdlast_cancer        <- 0L #Cancer found at last screen
  
  df$costs                <- 0.0 #Total cost counter
  df$drug_costs           <- 0.0 #Preventative drug cost counter
  df$US_costs             <- 0.0 #Ultrasound cost counter
  df$MRI_costs            <- 0.0 #MRI cost counter
  df$costs_follow_up      <- 0.0 #Follow-up cost counter
  
  df$cd_age_record        <- NA_real_ #Cancer diagnosis age
  df$cd_stage             <- NA_real_ #Cancer stage
  df$cd_size              <- NA_real_ #Cancer size
  df$cd_screen_flag       <- NA_integer_ 
  df$cd_screen_sens       <- NA_real_ #Screen sensitivity (based on size)
  df$cd_screen_spec       <- NA_real_ #Screen specificity
  df$cd_mort_age_raw      <- df$Mort_age #Base age of death
  df$cd_CD_age            <- df$CD_age #Cancer diagnosis age
  df$cd_death_age         <- NA_real_ #Cancer death age
  df$cd_screen_number     <- NA_integer_ #Screen at which cancer found
  
  df$Ca_mort_age          <- df$Mort_age #Cancer death age
  df$incidence_age_record <- NA_real_ #Cancer diagnosis age
  df$stage_cat            <- NA_real_ #Cancer stage
  df$Ca_size_screen       <- NA_real_ #Cancer size at screen
  df$LY_counter           <- NA_real_ #Life year counter
  df$QALY_counter         <- NA_real_ #QALY counter
  
  # Pre-calculate discount lookup vector
  max_age     <- ceiling(max(df$Mort_age))
  disc_lookup <- 1 / (1 + discount_cost)^(seq(0, max_age - screen_startage))
  
  get_disc <- function(age_vec) {
    disc_lookup[pmax(floor(age_vec - screen_startage) + 1L, 1L)]
  }
  
  # 1. Pre-screening preventative drug costs
  if (PREVENTATIVE_DRUG) {
    drug_idx <- df$risk_group != 0
    if (any(drug_idx)) {
      ms   <- df$starting_menses_status[drug_idx]
      ttd  <- df$time_taking_drug[drug_idx]
      cl   <- course_length[ms]
      prop <- if (cost_in_full_courses) rep(1.0, sum(drug_idx)) else ttd / cl
      cd   <- cost_drug[ms]
      
      df$costs[drug_idx]      <- df$costs[drug_idx]      + prop * cd
      df$drug_costs[drug_idx] <- df$drug_costs[drug_idx] + prop * cd
    }
  }
  
  # 2. Main event loop
  last_screen_slot                                    <- max.col(!is.na(screen_mat),
                                                                 ties.method = "last")
  last_screen_slot[rowSums(!is.na(screen_mat)) == 0L] <- 0L
  
  #Loop over maximum number of screens
  for (s in seq_len(max_screens)) {
    
    #Stop if no-one is left alive without cancer
    if (!any(df$active)) break
    
    screen_age_vec <- screen_mat[, s]
    
    # 2a. Did anyone die before this screen event?
    death_idx <- df$active &
      !is.na(screen_age_vec) &
      df$Mort_age <= screen_age_vec &
      df$Mort_age <= df$CD_age
    
    if (any(death_idx)) {
      df$active[death_idx]        <- FALSE
      df$cd_age_record[death_idx] <- df$Mort_age[death_idx]
    }
    
    # 2b. Did anyone have an interval cancer detected before this screen?
    cd_idx <- df$active &
      !is.na(screen_age_vec) &
      df$CD_age < df$Mort_age &
      df$CD_age < screen_age_vec
    
    if (any(cd_idx)) {
      df$interval_ca[cd_idx]          <- 1L
      df$active[cd_idx]               <- FALSE
      df$incidence_age_record[cd_idx] <- df$CD_age[cd_idx]
      df$cd_age_record[cd_idx]        <- df$CD_age[cd_idx]
      df$cd_size[cd_idx]              <- df$CD_size[cd_idx]
      df$age[cd_idx]                  <- df$CD_age[cd_idx]
      
      # Lookup discount rate
      disc <- get_disc(df$CD_age[cd_idx])
      
      #Add cost of follow-up to total costs
      df$costs[cd_idx]           <- df$costs[cd_idx]           + cost_follow_up * disc
      #Add cost of follow-up to follow-up costs
      df$costs_follow_up[cd_idx] <- df$costs_follow_up[cd_idx] + cost_follow_up * disc
      
      #Work out cancer stage
      df$stage_cat[cd_idx] <- vec_stage_by_size(df$CD_size[cd_idx])
      
      #Re-calculate age of death due to cancer
      df$Ca_mort_age[cd_idx] <- vec_ca_survival_time(
        df$stage_cat[cd_idx],
        df$Mort_age[cd_idx],
        df$CD_age[cd_idx],
        df$CD_age[cd_idx]
      )
      
      # Add cost if DCIS
      dcis_cd <- which(cd_idx & !is.na(df$stage_cat) & df$stage_cat == 5L)
      if (length(dcis_cd) > 0) {
        disc_dcis         <- get_disc(df$CD_age[dcis_cd])
        df$costs[dcis_cd] <- df$costs[dcis_cd] + cost_DCIS * disc_dcis
      }
      
      # Calculate cancer costs if invasive cancer
      stage_flag <- as.integer(df$stage_cat[cd_idx] >= 3L)
      age_flag   <- as.integer(df$CD_age[cd_idx] >= 65)
      surv_yrs   <- pmin(round(df$Ca_mort_age[cd_idx] - df$CD_age[cd_idx]), 9L)
      
      tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
      tx_costs[!is.na(df$stage_cat[cd_idx]) & df$stage_cat[cd_idx] >= 5L] <- 0.0
      if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[cd_idx])
      
      # Add cancer costs to total costs
      df$costs[cd_idx]        <- df$costs[cd_idx] + tx_costs * disc
      #Record stage
      df$cd_stage[cd_idx]     <- df$stage_cat[cd_idx]
      #Record age of death (min of all-cause or cancer death)
      df$cd_death_age[cd_idx] <- pmin(df$Ca_mort_age[cd_idx], df$Mort_age[cd_idx])
      # Record life years
      df$LY_counter[cd_idx]   <- df$Ca_mort_age[cd_idx] - start_age
    }
    
    # 2c. Screen event
    act_idx <- df$active &
      !is.na(screen_age_vec) &
      df$Mort_age > screen_age_vec &
      df$CD_age   > screen_age_vec
    
    if (!any(act_idx)) next
    
    # Lookup discount rates
    disc_screen          <- rep(NA_real_, n)
    disc_screen[act_idx] <- get_disc(screen_age_vec[act_idx])
    
    #Add a screen to count
    df$screen_count[act_idx] <- df$screen_count[act_idx] + 1L
    #Add the cost of a screen to total costs
    df$costs[act_idx]        <- df$costs[act_idx] + cost_screen * disc_screen[act_idx]
    
    #If stratified screening used add cost of stratification on first screen
    strat_idx <- act_idx &
      df$screen_count == 1L &
      screen_strategy %in% c(1L, 2L, 7L, 8L, 9L) &
      df$risk_predicted == 1L
    df$costs[strat_idx] <- df$costs[strat_idx] + cost_strat * disc_screen[strat_idx]
    
    #If last screen add counter
    is_last_screen <- act_idx & (last_screen_slot == s)
    df$lastscreen_count[is_last_screen] <- 1L
    
    #Add costs for ultrasound or MRI if woman has these
    us_idx  <- act_idx & df$US_screen  == 1L
    mri_idx <- act_idx & df$MRI_screen == 1L
    
    df$US_count[us_idx]    <- df$US_count[us_idx]    + 1L
    df$costs[us_idx]       <- df$costs[us_idx]       + cost_US  * disc_screen[us_idx]
    df$US_costs[us_idx]    <- df$US_costs[us_idx]    + cost_US  * disc_screen[us_idx]
    
    df$MRI_count[mri_idx]  <- df$MRI_count[mri_idx]  + 1L
    df$costs[mri_idx]      <- df$costs[mri_idx]       + cost_MRI * disc_screen[mri_idx]
    df$MRI_costs[mri_idx]  <- df$MRI_costs[mri_idx]   + cost_MRI * disc_screen[mri_idx]
    
    # Tumour detection
    t_tumour   <- screen_age_vec - df$genage
    tumour_idx <- act_idx & !is.na(t_tumour) & t_tumour > 0
    
    #Calculate current cancer sizes
    if (any(tumour_idx)) {
      t_i  <- t_tumour[tumour_idx]
      gr_i <- df$growth_rate[tumour_idx]
      vol  <- Vm / (1 + ((Vm / Vc)^0.25 - 1) * exp(-0.25 * gr_i * t_i))^4
      df$Ca_size_screen[tumour_idx] <- 2 * (vol / (4/3 * pi))^(1/3)
      
      #Determine if cancer is detected by screen
      screen_results <- vec_screening_result(
        df$Ca_size_screen[tumour_idx],
        df$VDG[tumour_idx],
        df$MRI_screen[tumour_idx],
        df$US_screen[tumour_idx]
      )
      #Update cancer detection counters
      sr_detected <- screen_results[, 1]
      sr_mammo    <- screen_results[, 2]
      sr_MRI      <- screen_results[, 3]
      sr_US       <- screen_results[, 4]
      
      tumour_rows   <- which(tumour_idx)
      detected_rows <- tumour_rows[sr_detected == 1L]
      
      if (length(detected_rows) > 0) {
        det_screen_age <- screen_age_vec[detected_rows]
        
        #Lookup discount rate
        det_disc <- get_disc(det_screen_age)
        
        #Update counters
        df$screen_detected_ca[detected_rows]   <- 1L
        df$active[detected_rows]               <- FALSE
        df$cd_age_record[detected_rows]        <- det_screen_age
        df$cd_size[detected_rows]              <- df$Ca_size_screen[detected_rows]
        df$cd_screen_flag[detected_rows]       <- 1L
        df$cd_screen_sens[detected_rows]       <- sr_MRI[sr_detected == 1L]
        df$cd_screen_spec[detected_rows]       <- sr_US[sr_detected == 1L]
        df$cd_screen_number[detected_rows]     <- df$screen_count[detected_rows]
        df$incidence_age_record[detected_rows] <- det_screen_age
        df$age[detected_rows]                  <- det_screen_age
        
        #Add follow-up costs
        df$costs[detected_rows]           <- df$costs[detected_rows] +
          cost_follow_up * det_disc
        df$costs_follow_up[detected_rows] <- df$costs_follow_up[detected_rows] +
          cost_follow_up * det_disc
        
        first_det <- detected_rows[df$screen_count[detected_rows] == 1L]
        if (length(first_det) > 0) df$sdfirst_cancer[first_det] <- 1L
        if (s == max_screens)      df$sdlast_cancer[detected_rows] <- 1L
        
        #Determine cancer stage
        df$stage_cat[detected_rows] <- vec_stage_by_size(
          df$Ca_size_screen[detected_rows]
        )
        
        #Determine age of death
        df$Ca_mort_age[detected_rows] <- vec_ca_survival_time(
          df$stage_cat[detected_rows],
          df$Mort_age[detected_rows],
          det_screen_age,
          df$CD_age[detected_rows]
        )
        
        # Add DCIS costs
        dcis_det <- detected_rows[
          !is.na(df$stage_cat[detected_rows]) & df$stage_cat[detected_rows] == 5L
        ]
        if (length(dcis_det) > 0) {
          det_disc_dcis      <- get_disc(screen_age_vec[dcis_det])
          df$costs[dcis_det] <- df$costs[dcis_det] + cost_DCIS * det_disc_dcis
        }
        
        # Calculate invasive cancer costs
        stage_flag <- as.integer(df$stage_cat[detected_rows] >= 3L)
        age_flag   <- as.integer(df$age[detected_rows] >= 65)
        surv_yrs   <- pmin(round(df$Ca_mort_age[detected_rows] -
                                   df$age[detected_rows]), 9L)
        
        tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
        tx_costs[!is.na(df$stage_cat[detected_rows]) &
                   df$stage_cat[detected_rows] >= 5L] <- 0.0
        if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[detected_rows])
        
        #Add invasive cancer treatment costs
        df$costs[detected_rows]        <- df$costs[detected_rows] + tx_costs * det_disc
        df$cd_stage[detected_rows]     <- df$stage_cat[detected_rows]
        df$cd_death_age[detected_rows] <- pmin(df$Ca_mort_age[detected_rows],
                                               df$Mort_age[detected_rows])
        df$LY_counter[detected_rows]   <- df$Ca_mort_age[detected_rows] - start_age
      }
    }
    
    #False-positive recalls
    fp_idx     <- act_idx & df$screen_detected_ca == 0L
    fp_draw    <- dqrunif(n, 0, 1)
    recall_idx <- fp_idx & fp_draw < recall_rate
    
    #Add false-positive costs and add to counters
    if (any(recall_idx)) {
      fp_cost <- (cost_follow_up + biopsy_rate * cost_biop) * disc_screen[recall_idx]
      df$recall_count[recall_idx]    <- df$recall_count[recall_idx]    + 1L
      df$costs[recall_idx]           <- df$costs[recall_idx]           + fp_cost
      df$costs_follow_up[recall_idx] <- df$costs_follow_up[recall_idx] + fp_cost
    }
    
    #Update column of active women in simulation
    still_active <- act_idx & df$active
    df$age[still_active] <- screen_age_vec[still_active]
    
  } # end screen loop
  
  # 3. Resolve women still active after all screens
  
  # 3a. Clinical detection after last screen
  post_cd_idx <- df$active & df$CD_age < df$Mort_age
  
  if (any(post_cd_idx)) {
    df$interval_ca[post_cd_idx]          <- 1L
    df$active[post_cd_idx]               <- FALSE
    df$incidence_age_record[post_cd_idx] <- df$CD_age[post_cd_idx]
    df$cd_age_record[post_cd_idx]        <- df$CD_age[post_cd_idx]
    df$cd_size[post_cd_idx]              <- df$CD_size[post_cd_idx]
    df$age[post_cd_idx]                  <- df$CD_age[post_cd_idx]
    
    # Lookup discount rate
    disc <- get_disc(df$CD_age[post_cd_idx])
    
    #Add follow-up costs
    df$costs[post_cd_idx]           <- df$costs[post_cd_idx]           + cost_follow_up * disc
    df$costs_follow_up[post_cd_idx] <- df$costs_follow_up[post_cd_idx] + cost_follow_up * disc
    
    #Calculate cancer stage
    df$stage_cat[post_cd_idx] <- vec_stage_by_size(df$CD_size[post_cd_idx])
    
    # Calculate age of death
    df$Ca_mort_age[post_cd_idx] <- vec_ca_survival_time(
      df$stage_cat[post_cd_idx],
      df$Mort_age[post_cd_idx],
      df$CD_age[post_cd_idx],
      df$CD_age[post_cd_idx]
    )
    
    # Add DCIS costs
    dcis_post <- which(post_cd_idx & !is.na(df$stage_cat) & df$stage_cat == 5L)
    if (length(dcis_post) > 0) {
      disc_dcis           <- get_disc(df$CD_age[dcis_post])
      df$costs[dcis_post] <- df$costs[dcis_post] + cost_DCIS * disc_dcis
    }
    
    # Calulcate invasive cancer treatment costs
    stage_flag <- as.integer(df$stage_cat[post_cd_idx] >= 3L)
    age_flag   <- as.integer(df$CD_age[post_cd_idx] >= 65)
    surv_yrs   <- pmin(round(df$Ca_mort_age[post_cd_idx] -
                               df$CD_age[post_cd_idx]), 9L)
    
    tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
    tx_costs[!is.na(df$stage_cat[post_cd_idx]) &
               df$stage_cat[post_cd_idx] >= 5L] <- 0.0
    if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[post_cd_idx])
    
    #Add invasive cancer treatment costs
    df$costs[post_cd_idx]        <- df$costs[post_cd_idx] + tx_costs * disc
    df$cd_stage[post_cd_idx]     <- df$stage_cat[post_cd_idx]
    df$cd_death_age[post_cd_idx] <- pmin(df$Ca_mort_age[post_cd_idx],
                                         df$Mort_age[post_cd_idx])
    df$LY_counter[post_cd_idx]   <- df$Ca_mort_age[post_cd_idx] - start_age
  }
  
  # 3b. Cancer-free deaths
  df$cd_age_record[df$active] <- df$Mort_age[df$active]
  
  # 4. QALY and life-year calculations
  df$LY_counter <- df$Ca_mort_age - start_age
  
  #Create warning if women never had cancer in the simulation
  no_cancer <- is.na(df$incidence_age_record)
  if (any(no_cancer)) {
    warning(paste(sum(no_cancer), "women have NA incidence_age_record after",
                  "event loop — check CD_age and genage values for these rows."))
    df$incidence_age_record[no_cancer] <- 0
    df$stage_cat[no_cancer]            <- 1L
  }
  
  #Run QALY counter
  df$QALY_counter <- vec_QALY_counter(
    df$Ca_mort_age,
    df$incidence_age_record,
    df$stage_cat
  )
  
  # 5. Post-screening preventative drug block
  if (PREVENTATIVE_DRUG) {
    drug_idx    <- df$risk_group != 0
    uptake_draw <- dqrunif(n, 0, 1)
    
    if (any(drug_idx)) {
      df$time_taking_drug2 <- 0.0
      
      # OPT 2: bulk rexp draw replaces mapply loop
      rg <- df$risk_group[drug_idx]
      ms <- df$starting_menses_status[drug_idx]
      ud <- uptake_draw[drug_idx]
      
      cl <- course_length[ms]
      up <- uptake[cbind(rg, ms)]
      pe <- persistence[cbind(rg, ms)]
      
      took_drug <- ud < up
      exp_draw  <- rexp(sum(drug_idx), rate = pe)
      df$time_taking_drug2[drug_idx] <- ifelse(took_drug, pmin(exp_draw, cl), 0.0)
      
      prop <- if (cost_in_full_courses) rep(1.0, sum(drug_idx)) else
        df$time_taking_drug2[drug_idx] / cl
      cd   <- cost_drug[ms]
      
      #Add costs of taking preventative drugs
      df$costs[drug_idx]      <- df$costs[drug_idx]      + prop * cd
      df$drug_costs[drug_idx] <- df$drug_costs[drug_idx] + prop * cd
    }
  }
  
  # 6. Assemble output
  results <- data.frame(
    QALY                   = df$QALY_counter,
    Cost                   = df$costs,
    Screens                = df$screen_count,
    `Cancer Diagnosed Age` = df$cd_age_record,
    Cancer                 = df$screen_detected_ca + df$interval_ca,
    `screen detected`      = df$screen_detected_ca,
    alternative            = screen_strategy,
    `Growth rate`          = df$growth_rate,
    `Life Years`           = df$LY_counter,
    Stage                  = df$cd_stage,
    `Cancer Size`          = df$cd_size,
    `Death Age`            = df$cd_death_age,
    `Cancer Screen Number` = df$cd_screen_number,
    check.names            = FALSE
  )
  
  if (PSA == 1L) {
    psa_col_names <- c(
      "PSA_gamma_survival_1", "PSA_gamma_survival_2", "PSA_gamma_survival_3",
      "PSA_meta_survival_54", "PSA_meta_survival_74", "PSA_meta_survival_99",
      "PSA_beta_1", "PSA_beta_2",
      "PSA_VDG1_sen", "PSA_VDG2_sen", "PSA_VDG3_sen", "PSA_VDG4_sen",
      "PSA_MRI_cdr", "PSA_US_cdr",
      "PSA_log_norm_mean", "PSA_log_norm_sd",
      "PSA_eff_ana", "PSA_eff_tam", "PSA_dropout_ana", "PSA_dropout_tam",
      "PSA_uptake_1", "PSA_uptake_2",
      "PSA_cost_strat", "PSA_costvar",
      "PSA_util_1to3", "PSA_util_4",
      "PSA_costscreen", "PSA_cost_follow_up", "PSA_cost_biop",
      "PSA_cost_US", "PSA_cost_MRI", "PSA_cost_drug",
      "mcid"
    )
    results <- cbind(results, risk_df[, psa_col_names])
  }
  
  if (PSA == 0) {
    save(results,
         file = paste(det_output_path, "Determ_", screen_strategy, "_", ii,
                      ".Rdata", sep = ""))
  } else {
    save(results,
         file = paste(psa_output_path, "PSA_", screen_strategy, "_", ii,
                      ".Rdata", sep = ""))
  }
  
  return(results)
}