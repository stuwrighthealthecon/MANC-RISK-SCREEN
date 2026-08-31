#################Compute strategy-dependent drug matrices###################
##################Used for prevantative drugs############################
get_drug_matrices <- function(screen_strategy) {
  if (screen_strategy %in% c(1, 9)) {
    risk_red <- matrix(
      rep(c(ana_eff, tam_eff), 5), nrow = 5, ncol = 2
    )
    course_length <- c(5., 5.)
    uptake <- rbind(c(0.,0.), c(0.,0.), c(0.,0.), c(.71,.71), c(.71,.71))
    persistence <- matrix(
      rep(c(ana_dropout_rate, tam_dropout_rate), 5), nrow = 5, ncol = 2
    )
  } else if (screen_strategy == 2) {
    risk_red <- matrix(
      rep(c(ana_eff, tam_eff), 3), nrow = 3, ncol = 2
    )
    uptake <- rbind(c(0.,0.), c(0.,0.), c(.71,.71))
    persistence <- matrix(
      rep(c(ana_dropout_rate, tam_dropout_rate), 3), nrow = 3, ncol = 2
    )
    course_length <- c(5., 5.)
  } else if (screen_strategy %in% c(7, 8)) {
    risk_red    <- matrix(c(ana_eff, tam_eff, ana_eff, tam_eff), nrow = 2, ncol = 2)
    uptake      <- rbind(c(0.,0.), c(.71,.71))
    persistence <- matrix(
      c(ana_dropout_rate, tam_dropout_rate, ana_dropout_rate, tam_dropout_rate),
      nrow = 2, ncol = 2
    )
    course_length <- c(5., 5.)
  } else {
    risk_red      <- matrix(c(ana_eff, tam_eff), nrow = 1, ncol = 2)
    uptake        <- matrix(c(0., 0.),           nrow = 1, ncol = 2)
    persistence   <- matrix(
      c(ana_dropout_rate, tam_dropout_rate), nrow = 1, ncol = 2
    )
    course_length <- c(5., 5.)
  }
  list(risk_red = risk_red, uptake = uptake,
       persistence = persistence, course_length = course_length)
}

###############Derive screen_times for a given strategy + risk group#############
get_screen_times <- function(screen_strategy, risk_group) {
  
  # Check: NA or missing risk_group falls back to low risk
  if (is.na(risk_group)) {
    warning(sprintf("get_screen_times: NA risk_group for strategy %d, using low_risk_screentimes",
                    screen_strategy))
    return(low_risk_screentimes)
  }
  if (screen_strategy == 0) return(c(0))     # no screening
  # PROCAS: risk_cutoffs_procas has 5 cutpoints -> groups 1-6
  if (screen_strategy == 1) {
    if      (risk_group <= 3)                      return(low_risk_screentimes)
    else if (risk_group == 4)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  # Tertiles: groups 1-3
  if (screen_strategy == 2) {
    if      (risk_group == 1)                      return(low_risk_screentimes)
    else if (risk_group == 2)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  #3 yearly
  if (screen_strategy == 3)                        return(low_risk_screentimes)
  #2 yearly
  if (screen_strategy == 4)                        return(med_risk_screentimes)
  #5 yearly
  if (screen_strategy == 5)                        return(seq(screen_startage, screen_startage + 5*4, 5))
  #10 yearly
  if (screen_strategy == 6)                        return(seq(screen_startage, screen_startage + 10, 10))
  #Low risk 5 yearly
  if (screen_strategy == 7) {
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 5*4, 5))
    else                                           return(low_risk_screentimes)
  }
  #Low risk 6 yearly
  if (screen_strategy == 8) {
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 6*3, 6))
    else                                           return(low_risk_screentimes)
  }
  # Fully stratified: risk_cutoffs_procas -> groups 1-6
  if (screen_strategy == 9) {
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 5*4, 5))
    else if (risk_group %in% c(2, 3))              return(low_risk_screentimes)
    else if (risk_group == 4)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  return(low_risk_screentimes)  # fallback
}

#######################Assign risk groups for a given strategy##################
assign_risk_groups <- function(df, screen_strategy, MISCLASS) {
  risk_col <- if (MISCLASS) "tenyrrisk_est" else "tenyrrisk"
  
  #Return error if no risk groups column
  if (!risk_col %in% names(df))
    stop(sprintf("assign_risk_groups: column '%s' not found in df", risk_col))
  
  rv <- df[[risk_col]]
  rg <- rep(0L, nrow(df))
  
  if (screen_strategy %in% c(1, 9)) {
    rg <- 1L + findInterval(rv, risk_cutoffs_procas)
  } else if (screen_strategy == 2) {
    rg <- 1L + findInterval(rv, risk_cutoffs_tert)
  } else if (screen_strategy %in% c(7, 8)) {
    rg <- ifelse(rv < low_risk_cut, 1L, 2L)
  }
  
  #Return error if NA values in risk group allocation
  if (any(is.na(rg)))
    warning(sprintf("assign_risk_groups: %d NA values in risk_group for strategy %d",
                    sum(is.na(rg)), screen_strategy))
  rg
}

#########################Lookup function for treatment costs############################

vec_fnLookupBase <- function(iStage_vec, iAge_vec, iLE_vec) {
  # Build a single composite key for both the lookup table and the query
  lookup_key <- paste(tblLookup$Stage, tblLookup$Age, tblLookup$Yr, sep = "|")
  #Create key for specific cancer to look up
  query_key  <- paste(iStage_vec,      iAge_vec,      iLE_vec,      sep = "|")
  #Look up total discounted cost for cancer
  as.numeric(tblLookup$CDCost.p.i.d[match(query_key, lookup_key)])
}

########################Stage calculator#######################################

vec_stage_by_size <- function(Ca_size_vec) {
  n <- length(Ca_size_vec)
  
  # Metastatic probability lookup
  # Bin each tumour size to the nearest metastatic_prob breakpoint
  m_size <- ifelse(Ca_size_vec <= 25,
                   25,
                   pmin(ceiling((Ca_size_vec - 25) / 10) * 10 + 25, 85))
  
  met_prob <- metastatic_prob[match(m_size, metastatic_prob[, 1]), 2]
  
  # Bulk random draws
  draw_met   <- dqrunif(n, 0, 1)   # one draw per woman for metastatic check
  draw_stage <- dqrunif(n, 0, 1)   # one draw per woman for stage sampling
  
  # Metastatic assignment
  stage_cat <- integer(n)
  is_met    <- draw_met < met_prob
  stage_cat[is_met] <- 4L
  
  # Stage sampling for non-metastatic women
  non_met  <- !is_met
  size_cat <- findInterval(Ca_size_vec[non_met], ca_size_cut)
  
  # Build cumulative probability matrix for non-metastatic women
  # stage_by_size_mat rows correspond to size_cat; cols are probs for c(1,2,3,5)
  prob_mat <- stage_by_size_mat[size_cat, , drop = FALSE]
  cum_mat <- matrixStats::rowCumsums(prob_mat)
  
  # Use the pre-drawn uniform to select stage via row-wise interval lookup
  u        <- draw_stage[non_met]
  stage_choices <- c(1L, 2L, 3L, 5L)
  
  # For each woman, find which cumulative probability bin her draw falls in
  stage_idx <- rowSums(cum_mat < u) + 1L   # gives index 1-4 into stage_choices
  stage_cat[non_met] <- stage_choices[stage_idx]
  
  return(stage_cat)
}

#######################Screening test results simulation##########################

vec_screening_result <- function(Ca_size_vec, VDG_vec, MRI_vec, US_vec) {
  n <- length(Ca_size_vec)
  
  # Base mammographic sensitivity based on sizes of cancers
  logit_val   <- exp((Ca_size_vec - beta2) / beta1)
  base_sens   <- logit_val / (1 + logit_val)
  base_sens   <- pmin(base_sens, sensitivity_max)
  
  # Density adjustment via odds ratio
  dense_OR  <- (Sen_VDG[VDG_vec] / (1 - Sen_VDG[VDG_vec])) /
    (Sen_VDG_av      / (1 - Sen_VDG_av))
  sens_odds <- (base_sens / (1 - base_sens)) * dense_OR
  Sensitivity <- sens_odds / (1 + sens_odds)
  
  # Single bulk random draw (one per woman)
  rnd_1 <- dqrunif(n, 0, 1)
  
  # Is the cancer detected by mammography?
  Mammo_detected <- rnd_1 < Sensitivity          
  Screen_detected <- Mammo_detected
  
  # MRI supplemental (only for non-mammo-detected, MRI-eligible women)
  mri_eligible <- !Screen_detected & (MRI_vec == 1L)
  MRI_detected <- logical(n)
  
  if (any(mri_eligible)) {
    mri_odds <- (Sensitivity[mri_eligible] / (1 - Sensitivity[mri_eligible])) *
      ((MRI_cdr + Mammo_cdr) / Mammo_cdr)
    mri_sens <- mri_odds / (1 + mri_odds)
    MRI_detected[mri_eligible] <- rnd_1[mri_eligible] < mri_sens
    Screen_detected[mri_eligible] <- Screen_detected[mri_eligible] |
      MRI_detected[mri_eligible]
  }
  
  # US supplemental (only for non-mammo-detected, US-eligible women)
  us_eligible <- !Mammo_detected & (US_vec == 1L)
  US_detected <- logical(n)
  
  if (any(us_eligible)) {
    us_odds <- (Sensitivity[us_eligible] / (1 - Sensitivity[us_eligible])) *
      ((US_cdr + Mammo_cdr) / Mammo_cdr)
    us_sens <- us_odds / (1 + us_odds)
    US_detected[us_eligible] <- rnd_1[us_eligible] < us_sens
    Screen_detected[us_eligible] <- Screen_detected[us_eligible] |
      US_detected[us_eligible]
  }
  
  # --- 7. Return matrix (rows = women, cols mirror original result vector) ---
  # col 1: Screen_detected, col 2: Mammo_detected,
  # col 3: MRI_detected,    col 4: US_detected
  cbind(
    as.integer(Screen_detected),
    as.integer(Mammo_detected),
    as.integer(MRI_detected),
    as.integer(US_detected)
  )
}

############################Simulate survival by stage##########################

vec_ca_survival_time <- function(stage_cat_vec, Mort_age_vec, age_vec, ca_incidence_age_vec) {
  n <- length(stage_cat_vec)
  
  # Pre-allocate output
  result_age <- numeric(n)
  
  # --- Stage masks ---
  is_early <- stage_cat_vec < 4L
  is_met   <- stage_cat_vec == 4L
  is_dcis  <- stage_cat_vec == 5L
  
  # BRANCH 1: Non-metastatic (stages 1-3)
  if (any(is_early)) {
    idx       <- which(is_early)
    n_early   <- length(idx)
    stage_i   <- stage_cat_vec[idx]
    age_i     <- age_vec[idx]
    inc_age_i <- ca_incidence_age_vec[idx]
    mort_i    <- Mort_age_vec[idx]
    
    # Base exponential survival draw 
    rate_base <- gamma_stage[stage_i]
    u1        <- dqrunif(n_early, 0, 1)
    surv_time <- -log(u1) / rate_base
    
    #Age >65 mortality adjustment
    old_idx <- idx[inc_age_i > 65]
    if (length(old_idx) > 0) {
      inc_old   <- ca_incidence_age_vec[old_idx]
      stage_old <- stage_cat_vec[old_idx]
      
      mort_row  <- pmin(floor(inc_old) + 1L, 100L)
      mort_mult <- Incidence_Mortality$X10year.mort.prob[mort_row] /
        Incidence_Mortality$X10year.mort.prob[66]
      rate_old  <- mort_mult * gamma_stage[stage_old]
      
      u1_old              <- dqrunif(length(old_idx), 0, 1)
      surv_time[inc_age_i > 65] <- -log(u1_old) / rate_old
    }
    
    # Survival > 10 years: switch to population Weibull mortality
    long_mask <- surv_time > 10
    if (any(long_mask)) {
      long_idx  <- idx[long_mask]
      inc_long  <- ca_incidence_age_vec[long_idx]
      
      p_lower <- pweibull(inc_long + 10,
                          shape = acmmortality_wb_a,
                          scale = acmmortality_wb_b)
      u2      <- dqrunif(length(long_idx), min = 0, max = 1)
      p_draw  <- p_lower + u2 * (1 - p_lower) 
      
      new_mort <- qweibull(p_draw,
                           shape = acmmortality_wb_a,
                           scale = acmmortality_wb_b)
      new_mort <- pmin(new_mort, time_horizon)
      
      # surv_time for these women is implicitly replaced via result_age below
      result_age[long_idx] <- pmin(new_mort, time_horizon)
      
      # Mark as resolved so we don't overwrite in the final step
      is_early[long_idx] <- FALSE
    }
    
    # Women not resolved by Weibull path
    still_early <- which(is_early & stage_cat_vec < 4L)
    result_age[still_early] <- ca_incidence_age_vec[still_early] +
      surv_time[match(still_early, idx)]
  }
  
  # BRANCH 2: Metastatic (stage 4)
  if (any(is_met)) {
    idx     <- which(is_met)
    age_i   <- age_vec[idx]
    
    age_cat_M        <- integer(length(idx))
    age_cat_M[age_i < 55]              <- 1L
    age_cat_M[age_i >= 55 & age_i < 75] <- 2L
    age_cat_M[age_i >= 75]             <- 3L
    
    rate_met       <- metastatic_survival[age_cat_M]
    u_met          <- dqrunif(length(idx), 0, 1)
    surv_time_met  <- -log(u_met) / rate_met
    
    result_age[idx] <- ca_incidence_age_vec[idx] + surv_time_met
  }
  
  # BRANCH 3: DCIS (stage 5) — no cancer mortality effect
  if (any(is_dcis)) {
    idx             <- which(is_dcis)
    result_age[idx] <- Mort_age_vec[idx]
  }
  
  # Final cap at time horizon and competing-risk floor
  result_age <- pmin(result_age, time_horizon)
  result_age <- pmin(result_age, 
                     pmax(result_age, ca_incidence_age_vec))  # never before incidence
  
  return(result_age)
}

#############################################QALY Counter##################
vec_QALY_counter_core <- function(Mort_age_vec, incidence_age_record_vec, stage_cat_vec) {
  n          <- length(Mort_age_vec)
  max_years  <- ceiling(max(Mort_age_vec)) - (screen_startage - 1)
  
  # Row = woman, col = year index y (1-indexed from screen_startage)
  # Initialised to 0
  QALY_mat <- matrix(0.0, nrow = n, ncol = max_years)
  
  # Fill with discounted age-related utility values
  y_idx      <- seq_len(max_years)
  
  # Age at each column: (screen_startage - 1) + y
  age_at_y   <- pmin(ceiling((screen_startage - 1) + y_idx),
                     max(utility_ages[, 1]))
  util_at_y  <- utility_ages[match(age_at_y, utility_ages[, 1]), 2]
  disc_at_y  <- 1 / (1 + discount_health)^y_idx
  
  # Base QALY weight for each column (same row pattern for all women)
  base_weight <- util_at_y * disc_at_y
  
  # Fill matrix: each row gets base_weight, then we mask/trim per-woman
  QALY_mat <- matrix(base_weight, nrow = n, ncol = max_years, byrow = TRUE)
  
  # --- Per-woman QALY_length and final-year partial adjustment
  QALY_length <- pmax(ceiling(Mort_age_vec) - (screen_startage - 1L), 1L)
  frac_last   <- 1 - (ceiling(Mort_age_vec) - Mort_age_vec)
  
  # Zero out columns beyond each woman's life and apply partial-year fraction
  # to her final column
  col_idx <- matrix(seq_len(max_years), nrow = n, ncol = max_years, byrow = TRUE)
  QALY_mat[col_idx > QALY_length] <- 0.0
  
  last_cell <- cbind(seq_len(n), QALY_length)
  QALY_mat[last_cell] <- QALY_mat[last_cell] * frac_last
  
  # Cancer utility adjustments
  has_cancer <- incidence_age_record_vec > 0
  
  if (any(has_cancer)) {
    ci  <- which(has_cancer)
    iar <- incidence_age_record_vec[ci]
    sc  <- stage_cat_vec[ci]
    ma  <- Mort_age_vec[ci]
    
    frac_into_year <- iar - floor(iar)           #to deal with months into year when cancer occurs
    col_y1         <- floor(iar) - screen_startage      # year-1 column index
    col_y2         <- col_y1 + 1L                       # year-2 column index
    
    u_y1     <- utility_stage_cat_y1[sc]
    u_follow <- utility_stage_cat_follow[sc]
    
    # Partial year at incidence (year 1 of cancer)
    valid1 <- col_y1 >= 1L & col_y1 <= max_years
    if (any(valid1)) {
      cells1 <- cbind(ci[valid1], col_y1[valid1])
      QALY_mat[cells1] <- u_y1[valid1] * QALY_mat[cells1] * (1 - frac_into_year[valid1])
    }
    
    # Transition year (straddles y1 and follow-up utility)
    long_enough <- (ma - iar) > 1
    if (any(long_enough)) {
      ci2 <- ci[long_enough]; c2v <- col_y2[long_enough]
      fv  <- frac_into_year[long_enough]
      u1v <- u_y1[long_enough]; uFv <- u_follow[long_enough]
      
      valid2 <- c2v >= 1L & c2v <= max_years
      if (any(valid2)) {
        cells2 <- cbind(ci2[valid2], c2v[valid2])
        QALY_mat[cells2] <- (u1v[valid2] * QALY_mat[cells2] * fv[valid2]) +
          (uFv[valid2] * QALY_mat[cells2] * (1 - fv[valid2]))
      }
    }
    
    #Follow-up years (y+2 to min(y+8, Mort_age))
    mort_cap <- pmin(ma, 100)
    has_followup <- ceiling(mort_cap) > (iar + 2)
    
    if (any(has_followup)) {
      ci3  <- ci[has_followup]
      iar3 <- iar[has_followup]; mc3 <- mort_cap[has_followup]; uF3 <- u_follow[has_followup]
      
      col_start3 <- floor(iar3) + 2L - screen_startage
      col_end3   <- pmin(floor(iar3) + 8L, ceiling(mc3)) - screen_startage
      
      col_idx_sub <- matrix(seq_len(max_years), nrow = length(ci3), ncol = max_years, byrow = TRUE)
      mask   <- col_idx_sub >= col_start3 & col_idx_sub <= col_end3
      uF_full <- matrix(uF3, nrow = length(ci3), ncol = max_years)  # recycles uF3 down each column, i.e. per row
      
      QALY_mat[ci3, ] <- QALY_mat[ci3, ] * ifelse(mask, uF_full, 1)
    }
  }
  #Return row sums (one total QALY per woman)
  rowSums(QALY_mat)
}

########################Split QALY calculation when n is very large#########
vec_QALY_counter <- function(Mort_age_vec, incidence_age_record_vec,
                             stage_cat_vec, chunk_size = 10000) {
  n <- length(Mort_age_vec)
  
  if (n <= chunk_size) {
    return(vec_QALY_counter_core(Mort_age_vec, incidence_age_record_vec,
                                 stage_cat_vec))
  }
  
  # Process in chunks to keep matrix sizes manageable
  result <- numeric(n)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  
  for (idx in chunks) {
    result[idx] <- vec_QALY_counter_core(
      Mort_age_vec[idx],
      incidence_age_record_vec[idx],
      stage_cat_vec[idx]
    )
  }
  result
}
