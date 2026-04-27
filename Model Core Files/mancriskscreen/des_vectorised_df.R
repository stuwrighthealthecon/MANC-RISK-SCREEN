# =============================================================================
# VECTORISED DES - Breast Cancer Screening Model (base R data.frame)
# =============================================================================
# Optimisations applied vs previous version:
#   1. last_screen_slot   - apply() replaced with max.col()
#   2. Drug persistence   - mapply() replaced with bulk rexp() draw
#   3. screen_mat         - k-loop replaced with matrix multiply broadcast
#   4. QALY trim loop     - for-loop replaced with vectorised matrix masking
#   5. disc recomputation - precomputed lookup vector used throughout
#   6. iStage/iAge        - integer flag keys replace character ifelse()
#
# Vectorised helper functions:
#   vec_stage_by_size(), vec_screening_result(), vec_ca_survival_time(),
#   vec_fnLookupBase(), vec_QALY_counter()
#
# Pre-requisite: run once at model initialisation before calling the DES:
#   tblLookup$stage_flag <- ifelse(tblLookup$Stage == "Early", 0L, 1L)
#   tblLookup$age_flag   <- ifelse(tblLookup$Age   == "18.64", 0L, 1L)
#   tblLookup$lookup_key <- paste(tblLookup$stage_flag,
#                                  tblLookup$age_flag,
#                                  tblLookup$Yr, sep = "|")
# =============================================================================


# =============================================================================
# 1. vec_stage_by_size
# =============================================================================
vec_stage_by_size <- function(Ca_size_vec) {
  n <- length(Ca_size_vec)
  
  # Metastatic probability lookup
  m_size   <- ifelse(Ca_size_vec <= 25,
                     25,
                     pmin(ceiling((Ca_size_vec - 25) / 10) * 10 + 25, 85))
  met_prob <- metastatic_prob[match(m_size, metastatic_prob[, 1]), 2]
  
  # Bulk random draws
  draw_met   <- dqrunif(n, 0, 1)
  draw_stage <- dqrunif(n, 0, 1)
  
  # Metastatic assignment
  stage_cat         <- integer(n)
  is_met            <- draw_met < met_prob
  stage_cat[is_met] <- 4L
  
  # Stage sampling for non-metastatic women
  non_met       <- !is_met
  size_cat      <- findInterval(Ca_size_vec[non_met], ca_size_cut)
  prob_mat      <- stage_by_size_mat[size_cat, , drop = FALSE]
  cum_mat       <- t(apply(prob_mat, 1, cumsum))
  u             <- draw_stage[non_met]
  stage_choices <- c(1L, 2L, 3L, 5L)
  stage_idx     <- rowSums(cum_mat < u) + 1L
  stage_cat[non_met] <- stage_choices[stage_idx]
  
  return(stage_cat)
}


# =============================================================================
# 2. vec_screening_result
# =============================================================================
vec_screening_result <- function(Ca_size_vec, VDG_vec, MRI_vec, US_vec) {
  n <- length(Ca_size_vec)
  
  # Base mammographic sensitivity (logistic, capped)
  logit_val   <- exp((Ca_size_vec - beta2) / beta1)
  base_sens   <- logit_val / (1 + logit_val)
  base_sens   <- pmin(base_sens, sensitivity_max)
  
  # Density adjustment via odds ratio
  dense_OR    <- (Sen_VDG[VDG_vec] / (1 - Sen_VDG[VDG_vec])) /
    (Sen_VDG_av      / (1 - Sen_VDG_av))
  sens_odds   <- (base_sens / (1 - base_sens)) * dense_OR
  Sensitivity <- sens_odds / (1 + sens_odds)
  
  # Single bulk random draw (one per woman, reused across modalities)
  rnd_1 <- dqrunif(n, 0, 1)
  
  # Mammography detection
  Mammo_detected  <- rnd_1 < Sensitivity
  Screen_detected <- Mammo_detected
  
  # MRI supplemental (gated on no mammography detection)
  mri_eligible <- !Screen_detected & (MRI_vec == 1L)
  MRI_detected <- logical(n)
  if (any(mri_eligible)) {
    mri_odds <- (Sensitivity[mri_eligible] / (1 - Sensitivity[mri_eligible])) *
      ((MRI_cdr + Mammo_cdr) / Mammo_cdr)
    mri_sens <- mri_odds / (1 + mri_odds)
    MRI_detected[mri_eligible]    <- rnd_1[mri_eligible] < mri_sens
    Screen_detected[mri_eligible] <- Screen_detected[mri_eligible] |
      MRI_detected[mri_eligible]
  }
  
  # US supplemental (gated on mammography miss, independent of MRI result)
  us_eligible <- !Mammo_detected & (US_vec == 1L)
  US_detected <- logical(n)
  if (any(us_eligible)) {
    us_odds <- (Sensitivity[us_eligible] / (1 - Sensitivity[us_eligible])) *
      ((US_cdr + Mammo_cdr) / Mammo_cdr)
    us_sens <- us_odds / (1 + us_odds)
    US_detected[us_eligible]     <- rnd_1[us_eligible] < us_sens
    Screen_detected[us_eligible] <- Screen_detected[us_eligible] |
      US_detected[us_eligible]
  }
  
  # Return n x 4 integer matrix
  # col 1: Screen_detected, col 2: Mammo_detected,
  # col 3: MRI_detected,    col 4: US_detected
  cbind(
    as.integer(Screen_detected),
    as.integer(Mammo_detected),
    as.integer(MRI_detected),
    as.integer(US_detected)
  )
}


# =============================================================================
# 3. vec_ca_survival_time
# =============================================================================
vec_ca_survival_time <- function(stage_cat_vec, Mort_age_vec, age_vec,
                                 ca_incidence_age_vec) {
  n <- length(stage_cat_vec)
  
  result_age <- numeric(n)
  
  is_early <- stage_cat_vec < 4L
  is_met   <- stage_cat_vec == 4L
  is_dcis  <- stage_cat_vec == 5L
  
  # ---------------------------------------------------------------------------
  # BRANCH 1: Non-metastatic (stages 1-3)
  # ---------------------------------------------------------------------------
  if (any(is_early)) {
    idx       <- which(is_early)
    n_early   <- length(idx)
    stage_i   <- stage_cat_vec[idx]
    inc_age_i <- ca_incidence_age_vec[idx]
    
    # Base exponential survival draw
    rate_base <- gamma_stage[stage_i]
    u1        <- dqrunif(n_early, 0, 1)
    surv_time <- -log(u1) / rate_base
    
    # Age >65 mortality adjustment
    old_mask <- inc_age_i > 65
    if (any(old_mask)) {
      old_idx   <- idx[old_mask]
      inc_old   <- ca_incidence_age_vec[old_idx]
      stage_old <- stage_cat_vec[old_idx]
      
      mort_row  <- pmin(floor(inc_old) + 1L, 100L)
      mort_mult <- Incidence_Mortality$X10year.mort.prob[mort_row] /
        Incidence_Mortality$X10year.mort.prob[66]
      rate_old  <- mort_mult * gamma_stage[stage_old]
      
      u1_old              <- dqrunif(length(old_idx), 0, 1)
      surv_time[old_mask] <- -log(u1_old) / rate_old
    }
    
    # Survival > 10 years: switch to population Weibull mortality
    long_mask <- surv_time > 10
    if (any(long_mask)) {
      long_idx <- idx[long_mask]
      inc_long <- ca_incidence_age_vec[long_idx]
      
      p_lower  <- pweibull(inc_long + 10,
                           shape = acmmortality_wb_a,
                           scale = acmmortality_wb_b)
      u2       <- dqrunif(length(long_idx), 0, 1)
      p_draw   <- p_lower + u2 * (1 - p_lower)
      
      new_mort <- qweibull(p_draw,
                           shape = acmmortality_wb_a,
                           scale = acmmortality_wb_b)
      new_mort <- pmin(new_mort, time_horizon)
      
      result_age[long_idx] <- new_mort
      is_early[long_idx]   <- FALSE
    }
    
    # Women not resolved by Weibull path
    still_early <- which(is_early & stage_cat_vec < 4L)
    result_age[still_early] <- ca_incidence_age_vec[still_early] +
      surv_time[match(still_early, idx)]
  }
  
  # ---------------------------------------------------------------------------
  # BRANCH 2: Metastatic (stage 4)
  # ---------------------------------------------------------------------------
  if (any(is_met)) {
    idx   <- which(is_met)
    age_i <- age_vec[idx]
    
    age_cat_M                            <- integer(length(idx))
    age_cat_M[age_i < 55]               <- 1L
    age_cat_M[age_i >= 55 & age_i < 75] <- 2L
    age_cat_M[age_i >= 75]              <- 3L
    
    rate_met      <- metastatic_survival[age_cat_M]
    u_met         <- dqrunif(length(idx), 0, 1)
    surv_time_met <- -log(u_met) / rate_met
    
    result_age[idx] <- ca_incidence_age_vec[idx] + surv_time_met
  }
  
  # ---------------------------------------------------------------------------
  # BRANCH 3: DCIS (stage 5) — no cancer mortality effect
  # ---------------------------------------------------------------------------
  if (any(is_dcis)) {
    idx             <- which(is_dcis)
    result_age[idx] <- Mort_age_vec[idx]
  }
  
  # Cap at time horizon
  result_age <- pmin(result_age, time_horizon)
  
  return(result_age)
}


# =============================================================================
# 4. vec_fnLookupBase
# =============================================================================
# Expects tblLookup$lookup_key pre-computed at initialisation using integer
# stage_flag (0=Early, 1=Late) and age_flag (0=18.64, 1=65plus) keys.
vec_fnLookupBase <- function(stage_flag_vec, age_flag_vec, iLE_vec) {
  query_key <- paste(stage_flag_vec, age_flag_vec, iLE_vec, sep = "|")
  as.numeric(tblLookup$CDCost.p.i.d[match(query_key, tblLookup$lookup_key)])
}


# =============================================================================
# 5. vec_QALY_counter
# =============================================================================
vec_QALY_counter <- function(Mort_age_vec, incidence_age_record_vec, stage_cat_vec) {
  n <- length(Mort_age_vec)
  
  # --- 1. Precompute cumulative discounted utility vector (no n×T matrix) ---
  mort_yr    <- ceiling(Mort_age_vec)
  frac_last  <- 1 - (mort_yr - Mort_age_vec)
  years_lived <- pmax(mort_yr - (screen_startage - 1L), 1L)
  
  max_years  <- max(years_lived)
  y_idx      <- seq_len(max_years)
  age_at_y   <- pmin(ceiling((screen_startage - 1) + y_idx), max(utility_ages[, 1]))
  util_at_y  <- utility_ages[match(age_at_y, utility_ages[, 1]), 2]
  disc_at_y  <- 1 / (1 + discount_health)^y_idx
  base_weight <- util_at_y * disc_at_y
  cum_base    <- cumsum(base_weight)   # length max_years — O(T) not O(n*T)
  
  # Each woman's base QALY = sum of base_weight[1..years_lived],
  # with the final year prorated by frac_last
  full_yrs   <- years_lived - 1L  # complete years before the partial final year
  qaly_base  <- ifelse(full_yrs >= 1L, cum_base[full_yrs], 0) +
    frac_last * base_weight[years_lived]
  
  # --- 2. Cancer utility adjustments (operates only on cancer subset) ---
  qaly_total <- qaly_base
  has_cancer <- incidence_age_record_vec > 0
  
  if (any(has_cancer)) {
    ci  <- which(has_cancer)
    iar <- incidence_age_record_vec[ci]
    sc  <- stage_cat_vec[ci]
    ma  <- Mort_age_vec[ci]
    
    frac_into_year <- iar - floor(iar)
    col_y1         <- floor(iar) - (screen_startage - 1L)   # 1-based year index at incidence
    col_y2         <- col_y1 + 1L
    
    u_y1     <- utility_stage_cat_y1[sc]
    u_follow <- utility_stage_cat_follow[sc]
    
    # 2a. Partial year at incidence
    # Woman loses (1 - u_y1) of the base_weight for the fraction of year AFTER incidence.
    # Weight for col_y1: normal for frac_into_year portion, u_y1-scaled for remainder.
    valid_y1 <- col_y1 >= 1L & col_y1 <= max_years
    if (any(valid_y1)) {
      v   <- which(valid_y1)
      w   <- base_weight[col_y1[v]]                # base weight for that year
      # Original full-year contribution already in qaly_base; apply delta
      delta_y1 <- w * (1 - frac_into_year[v]) * (u_y1[v] - 1)
      qaly_total[ci[v]] <- qaly_total[ci[v]] + delta_y1
    }
    
    # 2b. Transition year (col_y2): mix of y1 utility and follow-up utility
    long_enough <- (ma - iar) > 1
    if (any(long_enough)) {
      v   <- which(long_enough)
      c2v <- col_y2[v]
      valid_y2 <- c2v >= 1L & c2v <= max_years
      v2  <- v[valid_y2]
      if (length(v2) > 0) {
        c2    <- col_y2[v2]
        w     <- base_weight[c2]
        frac  <- frac_into_year[v2]
        # Transition year: frac portion at u_y1, remainder at u_follow
        mixed_util  <- frac * u_y1[v2] + (1 - frac) * u_follow[v2]
        delta_y2    <- w * (mixed_util - 1)
        qaly_total[ci[v2]] <- qaly_total[ci[v2]] + delta_y2
      }
    }
    
    # 2c. Follow-up years (incidence+2 to min(incidence+8, Mort_age))
    mort_cap     <- pmin(ma, 100)
    has_followup <- ceiling(mort_cap) > (iar + 2)
    
    if (any(has_followup)) {
      v <- which(has_followup)
      for (j in seq_along(v)) {
        jj      <- v[j]
        iar_j   <- iar[jj]
        mc_j    <- mort_cap[jj]
        uF_j    <- u_follow[jj]
        
        y_start <- floor(iar_j) + 2L
        y_end   <- min(floor(iar_j) + 8L, ceiling(mc_j))
        cols    <- (y_start:y_end) - (screen_startage - 1L)
        valid   <- cols >= 1L & cols <= max_years
        
        if (any(valid)) {
          cv <- cols[valid]
          delta <- base_weight[cv] * (uF_j - 1)
          qaly_total[ci[jj]] <- qaly_total[ci[jj]] + sum(delta)
        }
      }
    }
  }
  
  qaly_total
}


# =============================================================================
# 6. run_des_vectorised
# =============================================================================
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
  
  # ---------------------------------------------------------------------------
  # 0. Initialise working data.frame and extract per-woman screen schedules
  # ---------------------------------------------------------------------------
  df <- risk_df
  n  <- nrow(df)
  
  df$Mort_age <- df$life_expectancy
  df$CD_age   <- df$ca_incidence
  df$CD_size  <- df$clin_detect_size_g
  
  screen_cols <- grep("^screen_", names(df), value = TRUE)
  screen_cols <- screen_cols[order(as.integer(sub("screen_", "", screen_cols)))]
  attend_mat  <- as.matrix(df[, screen_cols])
  max_screens <- ncol(attend_mat)
  
  # OPT 3: matrix broadcast replaces k-loop for screen_mat construction
  screen_mat                  <- t(t(attend_mat) * screen_times)
  screen_mat[screen_mat == 0] <- NA_real_
  
  # ---------------------------------------------------------------------------
  # 0b. Initialise state and accumulator columns
  # ---------------------------------------------------------------------------
  df$age                  <- start_age
  df$active               <- TRUE
  df$interval_ca          <- 0L
  df$screen_detected_ca   <- 0L
  df$screen_count         <- 0L
  df$US_count             <- 0L
  df$MRI_count            <- 0L
  df$recall_count         <- 0L
  df$lastscreen_count     <- 0L
  df$sdfirst_cancer       <- 0L
  df$sdlast_cancer        <- 0L
  
  df$costs                <- 0.0
  df$drug_costs           <- 0.0
  df$US_costs             <- 0.0
  df$MRI_costs            <- 0.0
  df$costs_follow_up      <- 0.0
  
  df$cd_age_record        <- NA_real_
  df$cd_stage             <- NA_real_
  df$cd_size              <- NA_real_
  df$cd_screen_flag       <- NA_integer_
  df$cd_screen_sens       <- NA_real_
  df$cd_screen_spec       <- NA_real_
  df$cd_mort_age_raw      <- df$Mort_age
  df$cd_CD_age            <- df$CD_age
  df$cd_death_age         <- NA_real_
  df$cd_screen_number     <- NA_integer_
  
  df$Ca_mort_age          <- df$Mort_age
  df$incidence_age_record <- NA_real_
  df$stage_cat            <- NA_real_
  df$Ca_size_screen       <- NA_real_
  df$LY_counter           <- NA_real_
  df$QALY_counter         <- NA_real_
  
  # ---------------------------------------------------------------------------
  # OPT 5: Precompute discount lookup vector
  # Covers all possible ages from screen_startage to max competing-risk death
  # ---------------------------------------------------------------------------
  max_age     <- ceiling(max(df$Mort_age))
  disc_lookup <- 1 / (1 + discount_cost)^(seq(0, max_age - screen_startage))
  
  get_disc <- function(age_vec) {
    disc_lookup[pmax(floor(age_vec - screen_startage) + 1L, 1L)]
  }
  
  # ---------------------------------------------------------------------------
  # 1. Pre-screening preventative drug costs
  # ---------------------------------------------------------------------------
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
  
  # ---------------------------------------------------------------------------
  # 2. Main event loop
  # ---------------------------------------------------------------------------
  # OPT 1: max.col replaces apply() over rows
  last_screen_slot                                    <- max.col(!is.na(screen_mat),
                                                                 ties.method = "last")
  last_screen_slot[rowSums(!is.na(screen_mat)) == 0L] <- 0L
  
  for (s in seq_len(max_screens)) {
    
    if (!any(df$active)) break
    
    screen_age_vec <- screen_mat[, s]
    
    # ------------------------------------------------------------------
    # 2a. Death events before this screen slot
    # ------------------------------------------------------------------
    death_idx <- df$active &
      !is.na(screen_age_vec) &
      df$Mort_age <= screen_age_vec &
      df$Mort_age <= df$CD_age
    
    if (any(death_idx)) {
      df$active[death_idx]        <- FALSE
      df$cd_age_record[death_idx] <- df$Mort_age[death_idx]
    }
    
    # ------------------------------------------------------------------
    # 2b. Clinical detection before this screen slot
    # ------------------------------------------------------------------
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
      
      # OPT 5: discount lookup
      disc <- get_disc(df$CD_age[cd_idx])
      
      df$costs[cd_idx]           <- df$costs[cd_idx]           + cost_follow_up * disc
      df$costs_follow_up[cd_idx] <- df$costs_follow_up[cd_idx] + cost_follow_up * disc
      
      df$stage_cat[cd_idx] <- vec_stage_by_size(df$CD_size[cd_idx])
      
      df$Ca_mort_age[cd_idx] <- vec_ca_survival_time(
        df$stage_cat[cd_idx],
        df$Mort_age[cd_idx],
        df$CD_age[cd_idx],
        df$CD_age[cd_idx]
      )
      
      # DCIS costs
      dcis_cd <- which(cd_idx & !is.na(df$stage_cat) & df$stage_cat == 5L)
      if (length(dcis_cd) > 0) {
        disc_dcis         <- get_disc(df$CD_age[dcis_cd])
        df$costs[dcis_cd] <- df$costs[dcis_cd] + cost_DCIS * disc_dcis
      }
      
      # OPT 6: integer flag keys replace character ifelse vectors
      stage_flag <- as.integer(df$stage_cat[cd_idx] >= 3L)
      age_flag   <- as.integer(df$CD_age[cd_idx] >= 65)
      surv_yrs   <- pmin(round(df$Ca_mort_age[cd_idx] - df$CD_age[cd_idx]), 9L)
      
      tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
      tx_costs[!is.na(df$stage_cat[cd_idx]) & df$stage_cat[cd_idx] >= 5L] <- 0.0
      if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[cd_idx])
      
      df$costs[cd_idx]        <- df$costs[cd_idx] + tx_costs * disc
      df$cd_stage[cd_idx]     <- df$stage_cat[cd_idx]
      df$cd_death_age[cd_idx] <- pmin(df$Ca_mort_age[cd_idx], df$Mort_age[cd_idx])
      df$LY_counter[cd_idx]   <- df$Ca_mort_age[cd_idx] - start_age
    }
    
    # ------------------------------------------------------------------
    # 2c. Screen event
    # ------------------------------------------------------------------
    act_idx <- df$active &
      !is.na(screen_age_vec) &
      df$Mort_age > screen_age_vec &
      df$CD_age   > screen_age_vec
    
    if (!any(act_idx)) next
    
    # OPT 5: discount lookup
    disc_screen          <- rep(NA_real_, n)
    disc_screen[act_idx] <- get_disc(screen_age_vec[act_idx])
    
    df$screen_count[act_idx] <- df$screen_count[act_idx] + 1L
    df$costs[act_idx]        <- df$costs[act_idx] + cost_screen * disc_screen[act_idx]
    
    strat_idx <- act_idx &
      df$screen_count == 1L &
      screen_strategy %in% c(1L, 2L, 7L, 8L, 9L) &
      df$risk_predicted == 1L
    df$costs[strat_idx] <- df$costs[strat_idx] + cost_strat * disc_screen[strat_idx]
    
    is_last_screen <- act_idx & (last_screen_slot == s)
    df$lastscreen_count[is_last_screen] <- 1L
    
    us_idx  <- act_idx & df$US_screen  == 1L
    mri_idx <- act_idx & df$MRI_screen == 1L
    
    df$US_count[us_idx]    <- df$US_count[us_idx]    + 1L
    df$costs[us_idx]       <- df$costs[us_idx]       + cost_US  * disc_screen[us_idx]
    df$US_costs[us_idx]    <- df$US_costs[us_idx]    + cost_US  * disc_screen[us_idx]
    
    df$MRI_count[mri_idx]  <- df$MRI_count[mri_idx]  + 1L
    df$costs[mri_idx]      <- df$costs[mri_idx]       + cost_MRI * disc_screen[mri_idx]
    df$MRI_costs[mri_idx]  <- df$MRI_costs[mri_idx]   + cost_MRI * disc_screen[mri_idx]
    
    # ---- Tumour detection --------------------------------------------------
    t_tumour   <- screen_age_vec - df$genage
    tumour_idx <- act_idx & !is.na(t_tumour) & t_tumour > 0
    
    if (any(tumour_idx)) {
      t_i  <- t_tumour[tumour_idx]
      gr_i <- df$growth_rate[tumour_idx]
      vol  <- Vm / (1 + ((Vm / Vc)^0.25 - 1) * exp(-0.25 * gr_i * t_i))^4
      df$Ca_size_screen[tumour_idx] <- 2 * (vol / (4/3 * pi))^(1/3)
      
      screen_results <- vec_screening_result(
        df$Ca_size_screen[tumour_idx],
        df$VDG[tumour_idx],
        df$MRI_screen[tumour_idx],
        df$US_screen[tumour_idx]
      )
      sr_detected <- screen_results[, 1]
      sr_mammo    <- screen_results[, 2]
      sr_MRI      <- screen_results[, 3]
      sr_US       <- screen_results[, 4]
      
      tumour_rows   <- which(tumour_idx)
      detected_rows <- tumour_rows[sr_detected == 1L]
      
      if (length(detected_rows) > 0) {
        det_screen_age <- screen_age_vec[detected_rows]
        
        # OPT 5: discount lookup
        det_disc <- get_disc(det_screen_age)
        
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
        
        df$costs[detected_rows]           <- df$costs[detected_rows] +
          cost_follow_up * det_disc
        df$costs_follow_up[detected_rows] <- df$costs_follow_up[detected_rows] +
          cost_follow_up * det_disc
        
        first_det <- detected_rows[df$screen_count[detected_rows] == 1L]
        if (length(first_det) > 0) df$sdfirst_cancer[first_det] <- 1L
        if (s == max_screens)      df$sdlast_cancer[detected_rows] <- 1L
        
        df$stage_cat[detected_rows] <- vec_stage_by_size(
          df$Ca_size_screen[detected_rows]
        )
        
        df$Ca_mort_age[detected_rows] <- vec_ca_survival_time(
          df$stage_cat[detected_rows],
          df$Mort_age[detected_rows],
          det_screen_age,
          df$CD_age[detected_rows]
        )
        
        # DCIS costs
        dcis_det <- detected_rows[
          !is.na(df$stage_cat[detected_rows]) & df$stage_cat[detected_rows] == 5L
        ]
        if (length(dcis_det) > 0) {
          det_disc_dcis      <- get_disc(screen_age_vec[dcis_det])
          df$costs[dcis_det] <- df$costs[dcis_det] + cost_DCIS * det_disc_dcis
        }
        
        # OPT 6: integer flag keys
        stage_flag <- as.integer(df$stage_cat[detected_rows] >= 3L)
        age_flag   <- as.integer(df$age[detected_rows] >= 65)
        surv_yrs   <- pmin(round(df$Ca_mort_age[detected_rows] -
                                   df$age[detected_rows]), 9L)
        
        tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
        tx_costs[!is.na(df$stage_cat[detected_rows]) &
                   df$stage_cat[detected_rows] >= 5L] <- 0.0
        if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[detected_rows])
        
        df$costs[detected_rows]        <- df$costs[detected_rows] + tx_costs * det_disc
        df$cd_stage[detected_rows]     <- df$stage_cat[detected_rows]
        df$cd_death_age[detected_rows] <- pmin(df$Ca_mort_age[detected_rows],
                                               df$Mort_age[detected_rows])
        df$LY_counter[detected_rows]   <- df$Ca_mort_age[detected_rows] - start_age
      }
    }
    
    # ---- False-positive recalls --------------------------------------------
    fp_idx     <- act_idx & df$screen_detected_ca == 0L
    fp_draw    <- dqrunif(n, 0, 1)
    recall_idx <- fp_idx & fp_draw < recall_rate
    
    if (any(recall_idx)) {
      fp_cost <- (cost_follow_up + biopsy_rate * cost_biop) * disc_screen[recall_idx]
      df$recall_count[recall_idx]    <- df$recall_count[recall_idx]    + 1L
      df$costs[recall_idx]           <- df$costs[recall_idx]           + fp_cost
      df$costs_follow_up[recall_idx] <- df$costs_follow_up[recall_idx] + fp_cost
    }
    
    still_active <- act_idx & df$active
    df$age[still_active] <- screen_age_vec[still_active]
    
  } # end screen loop
  
  # ---------------------------------------------------------------------------
  # 3. Resolve women still active after all screens
  # ---------------------------------------------------------------------------
  
  # 3a. Clinical detection after last screen
  post_cd_idx <- df$active & df$CD_age < df$Mort_age
  
  if (any(post_cd_idx)) {
    df$interval_ca[post_cd_idx]          <- 1L
    df$active[post_cd_idx]               <- FALSE
    df$incidence_age_record[post_cd_idx] <- df$CD_age[post_cd_idx]
    df$cd_age_record[post_cd_idx]        <- df$CD_age[post_cd_idx]
    df$cd_size[post_cd_idx]              <- df$CD_size[post_cd_idx]
    df$age[post_cd_idx]                  <- df$CD_age[post_cd_idx]
    
    # OPT 5: discount lookup
    disc <- get_disc(df$CD_age[post_cd_idx])
    
    df$costs[post_cd_idx]           <- df$costs[post_cd_idx]           + cost_follow_up * disc
    df$costs_follow_up[post_cd_idx] <- df$costs_follow_up[post_cd_idx] + cost_follow_up * disc
    
    df$stage_cat[post_cd_idx] <- vec_stage_by_size(df$CD_size[post_cd_idx])
    
    df$Ca_mort_age[post_cd_idx] <- vec_ca_survival_time(
      df$stage_cat[post_cd_idx],
      df$Mort_age[post_cd_idx],
      df$CD_age[post_cd_idx],
      df$CD_age[post_cd_idx]
    )
    
    # DCIS costs
    dcis_post <- which(post_cd_idx & !is.na(df$stage_cat) & df$stage_cat == 5L)
    if (length(dcis_post) > 0) {
      disc_dcis           <- get_disc(df$CD_age[dcis_post])
      df$costs[dcis_post] <- df$costs[dcis_post] + cost_DCIS * disc_dcis
    }
    
    # OPT 6: integer flag keys
    stage_flag <- as.integer(df$stage_cat[post_cd_idx] >= 3L)
    age_flag   <- as.integer(df$CD_age[post_cd_idx] >= 65)
    surv_yrs   <- pmin(round(df$Ca_mort_age[post_cd_idx] -
                               df$CD_age[post_cd_idx]), 9L)
    
    tx_costs <- vec_fnLookupBase(stage_flag, age_flag, surv_yrs)
    tx_costs[!is.na(df$stage_cat[post_cd_idx]) &
               df$stage_cat[post_cd_idx] >= 5L] <- 0.0
    if (PSA == 1L) tx_costs <- tx_costs * (1 + df$PSA_costvar[post_cd_idx])
    
    df$costs[post_cd_idx]        <- df$costs[post_cd_idx] + tx_costs * disc
    df$cd_stage[post_cd_idx]     <- df$stage_cat[post_cd_idx]
    df$cd_death_age[post_cd_idx] <- pmin(df$Ca_mort_age[post_cd_idx],
                                         df$Mort_age[post_cd_idx])
    df$LY_counter[post_cd_idx]   <- df$Ca_mort_age[post_cd_idx] - start_age
  }
  
  # 3b. Cancer-free deaths
  df$cd_age_record[df$active] <- df$Mort_age[df$active]
  
  # ---------------------------------------------------------------------------
  # 4. QALY and life-year calculations
  # ---------------------------------------------------------------------------
  df$LY_counter <- df$Ca_mort_age - start_age
  
  no_cancer <- is.na(df$incidence_age_record)
  if (any(no_cancer)) {
    warning(paste(sum(no_cancer), "women have NA incidence_age_record after",
                  "event loop — check CD_age and genage values for these rows."))
    df$incidence_age_record[no_cancer] <- 0
    df$stage_cat[no_cancer]            <- 1L
  }
  
  df$QALY_counter <- vec_QALY_counter(
    df$Ca_mort_age,
    df$incidence_age_record,
    df$stage_cat
  )
  
  # ---------------------------------------------------------------------------
  # 5. Post-screening preventative drug block
  # ---------------------------------------------------------------------------
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
      
      df$costs[drug_idx]      <- df$costs[drug_idx]      + prop * cd
      df$drug_costs[drug_idx] <- df$drug_costs[drug_idx] + prop * cd
    }
  }
  
  # ---------------------------------------------------------------------------
  # 6. Assemble output
  # ---------------------------------------------------------------------------
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