# =============================================================================
# MANC-RISK-SCREEN Main Model Script
# =============================================================================
# Optimisations vs original:
#   1.  source() calls moved outside strategy loop
#   2.  load() moved outside strategy loop; master copy kept
#   3.  names() regex fix moved outside strategy loop
#   4.  dqset.seed() replaces set.seed() for dqrng compatibility
#   5.  ca_incidence vectorised with outer()/rowSums() replacing map_dbl()
#   6.  supplemental screening loop replaced with vectorised indexing + bug fix
#   7.  screen attendance uses direct column names not fragile match()
#   8.  risk group split() pre-computed once before ii loop
#   9.  outer() transpose removed in drug block
#   10. t(apply(..., cumsum)) replaced with matrixStats::rowCumsums()
#   11. pacman replaced with explicit library() calls
#   12. Strategy loop parallelised via foreach/doParallel
#   13. params.R split: fixed params sourced once, strategy-dependent params
#       computed per strategy inside loop
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
library(doParallel)
library(MASS)
library(dqrng)
library(compiler)
library(tidyverse)
library(iterators)
library(here)
library(purrr)        # fixed spelling from original 'purr'
library(tictoc)
library(matrixStats)  # for rowCumsums()

# -----------------------------------------------------------------------------
# 1. Controls
# -----------------------------------------------------------------------------
controls <- list(
  strategies       = c(0,1,2,3,4,9),
  gensample        = TRUE,
  MISCLASS         = TRUE,
  PREVENTATIVE_DRUG = FALSE,
  supplemental_screening = FALSE,
  PSA              = FALSE,
  intervals        = FALSE,
  desired_cases    = 3000,
  mcruns           = 1,
  seed             = 42,
  n_cores          = max(1L, parallel::detectCores() - 1L)
)

MISCLASS          <- controls$MISCLASS
PREVENTATIVE_DRUG <- controls$PREVENTATIVE_DRUG
PSA               <- as.integer(controls$PSA)
intervals         <- as.integer(controls$intervals)
#Place control on number of cores for memory management
controls$n_cores <- if (controls$desired_cases >= 200000) 3L else
  if (controls$desired_cases >= 100000) 6L else
    max(1L, parallel::detectCores() - 1L)

# -----------------------------------------------------------------------------
# 2. Output directories
# -----------------------------------------------------------------------------
det_output_path <- "Deterministic results/"
psa_output_path <- "PSA results/"

if (MISCLASS & PREVENTATIVE_DRUG) {
  dir.create("Deterministic results/misclassification_and_preventative_drug",
             showWarnings = FALSE)
  dir.create("PSA results/misclassification_and_preventative_drug",
             showWarnings = FALSE)
  det_output_path <- "Deterministic results/misclassification_and_preventative_drug/"
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
} else if (MISCLASS) {
  dir.create("Deterministic results/misclassification", showWarnings = FALSE)
  dir.create("PSA results/misclassification",           showWarnings = FALSE)
  det_output_path <- "Deterministic results/misclassification/"
  psa_output_path <- "PSA results/misclassification/"
} else if (PREVENTATIVE_DRUG) {
  dir.create("Deterministic results/preventative_drug",                    showWarnings = FALSE)
  dir.create("PSA results/misclassification_and_preventative_drug",        showWarnings = FALSE)
  det_output_path <- "Deterministic results/preventative_drug/"
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
}

sample_fname <- "possample"

# -----------------------------------------------------------------------------
# 3. Fixed parameters and functions — sourced ONCE outside all loops
# -----------------------------------------------------------------------------
expected_prev <- 0.12
desired_cases <- controls$desired_cases
mcruns        <- controls$mcruns

inum <- ceiling(desired_cases / expected_prev)

# OPT 4: dqset.seed controls dqrng draws; base set.seed does not
dqset.seed(controls$seed)

# Source function files once — they never change between strategies
source("MANC_RISK_SCREEN_functions.R")
source("risksample function.R")
source("negsample function.R")
source("des_vectorised_df.R")

# Source fixed params (everything in params.R that does NOT depend on
# params.R contains a strategy-dependent block that references screen_strategy.
# params.R now guards against this with exists("screen_strategy"), so it can
# be sourced safely here before the strategy loop.
source("params.R")

# -----------------------------------------------------------------------------
# 4. Generate sample (once, before strategy loop)
# -----------------------------------------------------------------------------
screen_strategies <- unlist(controls$strategies)

if (controls$gensample) {
  # Use first strategy for sample generation — sample is shared across strategies
  screen_strategy <- screen_strategies[1]

  if (MISCLASS) {
    dir.create("Risksamplewithmisclass", showWarnings = FALSE)
    create_sample_with_misclass(PSA, intervals, controls$seed, screen_strategy)
  } else {
    dir.create("Risksample", showWarnings = FALSE)
    cmp_create_sample(PSA, intervals, controls$seed, screen_strategy)
  }
}

# OPT 2: load sample ONCE outside the strategy loop
if (MISCLASS) {
  load(paste0("Risksamplewithmisclass/", sample_fname, ".Rdata"))
} else {
  load(paste0("Risksample/", sample_fname, ".Rdata"))
}

# OPT 3: fix column names ONCE after loading
prefix <- paste0("^", "X", 1, ".")
names(risksample) <- sub(prefix, "", names(risksample))

# OPT 5: Vectorised cancer incidence age assignment
# Replaces map_dbl loop over every woman
age_cols <- Incidence_Mortality$age  # length 101

# Matrix: rows = women, cols = ages. Zero out ages >= life_expectancy per woman.
age_mat  <- outer(rep(1, nrow(risksample)), Incidence_Mortality$BC_age) *
            (outer(risksample$life_expectancy, age_cols, FUN = ">"))
age_mat  <- age_mat / rowSums(age_mat)

cum_age_mat <- matrixStats::rowCumsums(age_mat)   # OPT 10
u_age       <- dqrunif(nrow(risksample), 0, 1)
col_idx     <- rowSums(cum_age_mat < u_age) + 1L
col_idx     <- pmin(col_idx, length(age_cols))
risksample$ca_incidence <- age_cols[col_idx]

# Add fractional month jitter
risksample$ca_incidence <- risksample$ca_incidence + dqrunif(nrow(risksample), 0, 1)

# Tumour size and genesis age
risksample$clin_detect_size_g <- start_size * 2^risksample$clinical_detect_size

t_gen <- ((log((Vm / Vc)^0.25 - 1) -
             log((Vm / ((4/3) * pi * (risksample$clin_detect_size_g / 2)^3))^0.25 - 1)) /
            (0.25 * risksample$growth_rate))
risksample$genage <- risksample$ca_incidence - t_gen

# Risk group assignment
if (MISCLASS) {
  risk_col <- "tenyrrisk_est"
} else {
  risk_col <- "tenyrrisk"
}

# Initialise risk_group — will be overwritten per-strategy inside loop
# for strategies that use risk stratification; set to 0 (average) as default
risksample$risk_group <- 0L

# OPT 6: vectorised supplemental screening assignment (fixes original bug)
if (controls$supplemental_screening) {
  supp_col        <- ifelse(MISCLASS, "tenyrrisk_true", "tenyrrisk_est")
  high_risk_dense <- risksample$VDG >= density_cutoff &
                     risksample[[supp_col]] >= 8
  low_risk_dense  <- risksample$VDG >= density_cutoff &
                     risksample[[supp_col]] < 8
  risksample$MRI_screen[high_risk_dense] <- 1L   # was bug: used < not <-
  risksample$US_screen[low_risk_dense]   <- 1L
}

# Drug setup for PREVENTATIVE_DRUG runs
if (PREVENTATIVE_DRUG) {
  nsample <- nrow(risksample)
  risksample$starting_menses_status <- ifelse(
    dqrunif(nsample, 0, 1) < prob_premen, 1L, 2L
  )
  risksample$takes_drug        <- logical(nsample)
  risksample$time_taking_drug  <- numeric(nsample)
}

# Keep a clean master copy for reuse across strategies
risksample_master <- risksample

# Keep risk columns in splitmaster_base so assign_risk_groups() can use them
# inside the strategy loop. Only drop columns that are truly never needed.
if (MISCLASS) {
  drop_cols_now  <- c("VBD", "cancer", "feedback", "liferisk_est", "liferisk_true")
  drop_cols_late <- c("tenyrrisk_est", "tenyrrisk_true")
} else {
  drop_cols_now  <- c("VBD", "cancer", "feedback", "liferisk")
  drop_cols_late <- c("tenyrrisk")
}
if (PSA == 0L) {
  psa_cols   <- grep("^PSA_", names(splitmaster), value = TRUE)
  splitmaster <- splitmaster[, !names(splitmaster) %in% psa_cols]
}
splitmaster_base <- risksample_master[, !names(risksample_master) %in% drop_cols_now]

# -----------------------------------------------------------------------------
# 5a. Helper: compute strategy-dependent drug matrices
# Replicates the screen_strategy-dependent block from params.R so workers
# can call this with their own screen_strategy value.
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

# 5b. Helper: derive screen_times for a given strategy + risk group
# -----------------------------------------------------------------------------
get_screen_times <- function(screen_strategy, risk_group) {
  # Guard: NA or missing risk_group falls back to low risk
  if (is.na(risk_group)) {
    warning(sprintf("get_screen_times: NA risk_group for strategy %d, using low_risk_screentimes",
                    screen_strategy))
    return(low_risk_screentimes)
  }
  if (screen_strategy == 0) return(c(0))     # no screening
  if (screen_strategy == 1) {
    # PROCAS: risk_cutoffs_procas has 5 cutpoints -> groups 1-6
    if      (risk_group <= 3)                      return(low_risk_screentimes)
    else if (risk_group == 4)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  if (screen_strategy == 2) {
    # Tertiles: groups 1-3
    if      (risk_group == 1)                      return(low_risk_screentimes)
    else if (risk_group == 2)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  if (screen_strategy == 3)                        return(low_risk_screentimes)
  if (screen_strategy == 4)                        return(med_risk_screentimes)
  if (screen_strategy == 5)                        return(seq(screen_startage, screen_startage + 5*4, 5))
  if (screen_strategy == 6)                        return(seq(screen_startage, screen_startage + 10, 10))
  if (screen_strategy == 7) {
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 5*4, 5))
    else                                           return(low_risk_screentimes)
  }
  if (screen_strategy == 8) {
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 6*3, 6))
    else                                           return(low_risk_screentimes)
  }
  if (screen_strategy == 9) {
    # Fully stratified: risk_cutoffs_procas -> groups 1-6
    if      (risk_group == 1)                      return(seq(screen_startage, screen_startage + 5*4, 5))
    else if (risk_group %in% c(2, 3))              return(low_risk_screentimes)
    else if (risk_group == 4)                      return(med_risk_screentimes)
    else                                           return(high_risk_screentimes)
  }
  return(low_risk_screentimes)  # fallback
}

# -----------------------------------------------------------------------------
# 6. Helper: assign risk groups for a given strategy
# -----------------------------------------------------------------------------
assign_risk_groups <- function(df, screen_strategy, MISCLASS) {
  risk_col <- if (MISCLASS) "tenyrrisk_est" else "tenyrrisk"

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

  if (any(is.na(rg)))
    warning(sprintf("assign_risk_groups: %d NA values in risk_group for strategy %d",
                    sum(is.na(rg)), screen_strategy))
  rg
}

# -----------------------------------------------------------------------------
# 7. Main strategy loop — parallelised via foreach
# -----------------------------------------------------------------------------
tic()

cl <- makeCluster(controls$n_cores)
registerDoParallel(cl)

# Export globals needed by workers
worker_globals <- c(
  # data frames / matrices
  "splitmaster_base", "Incidence_Mortality", "tblLookup",
  "metastatic_prob", "stage_by_size_mat", "utility_ages",
  # scalar / vector params
  "screen_startage", "screen_endage", "start_age", "time_horizon",
  "discount_cost", "discount_health",
  "acmmortality_wb_a", "acmmortality_wb_b",
  "gamma_stage", "metastatic_survival",
  "beta1", "beta2", "sensitivity_max",
  "Sen_VDG", "Sen_VDG_av", "Mammo_cdr", "MRI_cdr", "US_cdr",
  "density_cutoff", "ca_size_cut",
  "log_norm_mean", "log_norm_sd",
  "Vc", "Vm", "start_size", "max_size",
  "recall_rate", "biopsy_rate",
  "risk_cutoffs_procas", "risk_cutoffs_tert", "low_risk_cut",
  "cost_strat", "cost_screen", "cost_follow_up", "cost_biop",
  "cost_DCIS", "cost_US", "cost_MRI", "cost_drug",
  "cost_screen_base", "cost_follow_up_base", "cost_biop_base",
  "cost_DCIS_base", "cost_US_base", "cost_MRI_base", "cost_drug_base",
  "low_risk_screentimes", "med_risk_screentimes", "high_risk_screentimes",
  "uptakefirstscreen", "uptakeotherscreen", "uptakenoscreen",
  "utility_stage_cat_y1", "utility_stage_cat_follow", "utility_DCIS",
  "MISCLASS", "PREVENTATIVE_DRUG", "PSA", "intervals",
  "det_output_path", "psa_output_path", "mcruns",
  "course_length", "uptake", "persistence", "risk_red",
  "inc_scale", "inc_shape", "completion_prob", "prob_premen",
  "cost_in_full_courses", "drop_cols_late",
  "screen_strategies", "controls",
  # functions
  "run_des_vectorised", "vec_stage_by_size", "vec_screening_result",
  "vec_ca_survival_time", "vec_fnLookupBase",
  "vec_QALY_counter", "vec_QALY_counter_core",
  "get_screen_times", "assign_risk_groups",
  "negsamplefn", "redraw_drug_pars"
)

foreach(
  r         = seq_along(screen_strategies),
  .packages = c("dqrng", "matrixStats", "tidyverse"),
  .export   = worker_globals
) %dopar% {

  screen_strategy        <- screen_strategies[r]
  supplemental_screening <- as.integer(controls$supplemental_screening)

  # ---- Recompute strategy-dependent drug matrices --------------------------
  dm          <- get_drug_matrices(screen_strategy)
  risk_red    <- dm$risk_red
  uptake      <- dm$uptake
  persistence <- dm$persistence
  course_length <- dm$course_length

  # ---- Assign risk groups for this strategy --------------------------------
  splitmaster <- splitmaster_base
  splitmaster$risk_group <- assign_risk_groups(splitmaster, screen_strategy,
                                               MISCLASS)
  risk_groups <- sort(unique(splitmaster$risk_group))

  # Drop the risk score columns now that risk groups are assigned
  splitmaster <- splitmaster[, !names(splitmaster) %in% drop_cols_late]

  # OPT 8: pre-split data.frame once before ii loop
  risk_group_list <- split(splitmaster, splitmaster$risk_group)

  # ---- ii loop: one iteration per risk group --------------------------------
  for (ii in seq_along(risk_groups)) {

    risksample <- risk_group_list[[as.character(risk_groups[ii])]]

    # Defensive check — catch NA risk_group before it causes cryptic errors
    if (is.na(risksample$risk_group[1]))
      stop(sprintf("strategy %d ii %d: risk_group[1] is NA",
                   screen_strategy, ii))

    # ---- Drug incidence adjustment (PREVENTATIVE_DRUG only) ----------------
    if (PREVENTATIVE_DRUG && risksample$risk_group[1] != 0) {

      if (PSA == 1L) {
        drug_matrix_list <- redraw_drug_pars(risksample[1, ])
        risk_red    <- drug_matrix_list[[1]]
        uptake      <- drug_matrix_list[[2]]
        persistence <- drug_matrix_list[[3]]
      }

      idx         <- cbind(risksample$risk_group, risksample$starting_menses_status)
      uptake_probs <- uptake[idx]
      risksample$uptake <- dqrunif(nrow(risksample), 0, 1) < uptake_probs

      risksample$time_taking_drug <- risksample$uptake * pmin(
        rexp(nrow(risksample), rate = persistence[idx]),
        course_length
      )

      risksample$weibullrisk <- (1 - risk_red[idx]) *
        risksample$time_taking_drug *
        log(1 / completion_prob[risksample$starting_menses_status])
      risksample$weibullrisk <- inc_scale / risksample$weibullrisk

      ages <- Incidence_Mortality$age[start_age:101]
      n    <- nrow(risksample)

      # OPT 9: outer() argument order avoids transpose
      prob_matrix <- outer(
        risksample$weibullrisk, ages,
        FUN = function(scale, age) dweibull(age, shape = inc_shape, scale = scale)
      )

      # Zero out ages >= life_expectancy
      age_mask    <- outer(risksample$life_expectancy, ages, FUN = ">")
      prob_matrix <- prob_matrix * age_mask

      # Normalise rows
      row_sums    <- rowSums(prob_matrix)
      prob_matrix <- prob_matrix / row_sums

      # OPT 10: rowCumsums() replaces t(apply(..., cumsum))
      cum_probs   <- matrixStats::rowCumsums(prob_matrix)
      u           <- dqrunif(n, 0, 1)
      col_indices <- rowSums(cum_probs < u) + 1L
      col_indices <- pmin(col_indices, length(ages))

      risksample$ca_incidence <- ifelse(
        risksample$time_taking_drug > 0,
        ages[col_indices] + dqrunif(n, 0, 1),
        risksample$ca_incidence
      )
    }

    # ---- Cap incidence and life expectancy ----------------------------------
    risksample$ca_incidence <- ifelse(
      risksample$life_expectancy <= risksample$ca_incidence,
      floor(risksample$life_expectancy),
      risksample$ca_incidence
    )
    risksample$life_expectancy <- ifelse(
      risksample$life_expectancy >= time_horizon,
      99.99,
      risksample$life_expectancy
    )

    # ---- Screen times for this strategy + risk group -----------------------
    screen_times <- get_screen_times(screen_strategy, risksample$risk_group[1])

    # ---- Screen attendance columns -----------------------------------------
    if (length(screen_times) > 1) {
      n_screens <- length(screen_times)

      # Initialise columns one at a time
      for (k in seq_len(n_screens)) {
        col_k <- paste0("screen_", k)
        risksample[[col_k]] <- rep(0L, nrow(risksample))
      }

      # First screen attendance
      risksample[["screen_1"]] <- rbinom(nrow(risksample), 1L, uptakefirstscreen)

      ever_attended <- risksample[["screen_1"]]

      for (i in seq_len(n_screens - 1L)) {
        next_col <- paste0("screen_", i + 1L)
        new_vals <- ifelse(
          ever_attended >= 1L,
          rbinom(nrow(risksample), 1L, uptakeotherscreen),
          rbinom(nrow(risksample), 1L, uptakenoscreen)
        )
        risksample[[next_col]] <- new_vals
        ever_attended <- pmax(ever_attended, risksample[[next_col]])
      }

    } else {
      risksample[["screen_1"]] <- 999L
    }

    # ---- PSA parameter overrides -------------------------------------------
    if (PSA == 1L) {
      risk_data <- risksample[1, ]

      beta1 <- risk_data$PSA_beta_1
      beta2 <- risk_data$PSA_beta_2

      log_norm_mean <- risk_data$PSA_log_norm_mean
      log_norm_sd   <- risk_data$PSA_log_norm_sd

      gamma_stage <- c(
        exp(risk_data$PSA_gamma_survival_1),
        exp(risk_data$PSA_gamma_survival_2),
        exp(risk_data$PSA_gamma_survival_3)
      )
      metastatic_survival <- c(
        exp(risk_data$PSA_meta_survival_54),
        exp(risk_data$PSA_meta_survival_74),
        exp(risk_data$PSA_meta_survival_99)
      )

      Sen_VDG <- c(risk_data$PSA_VDG1_sen, risk_data$PSA_VDG2_sen,
                   risk_data$PSA_VDG3_sen, risk_data$PSA_VDG4_sen)
      Sen_VDG_av <- mean(Sen_VDG)

      MRI_cdr <- risk_data$PSA_MRI_cdr
      US_cdr  <- risk_data$PSA_US_cdr

      utility_stage_cat_y1 <- c(
        "stage1"     = risk_data$PSA_util_1to3 / 0.822,
        "stage2"     = risk_data$PSA_util_1to3 / 0.822,
        "stage3"     = risk_data$PSA_util_1to3 / 0.822,
        "Metastatic" = risk_data$PSA_util_4    / 0.822,
        "DCIS"       = utility_DCIS
      )
      utility_stage_cat_follow <- utility_stage_cat_y1

      cost_strat    <- risk_data$PSA_cost_strat
      cost_DCIS     <- cost_DCIS_base     * (1 + risk_data$PSA_costvar)
      cost_screen   <- cost_screen_base   * (1 + risk_data$PSA_costscreen)
      cost_follow_up <- cost_follow_up_base * (1 + risk_data$PSA_cost_follow_up)
      cost_biop     <- cost_biop_base     * (1 + risk_data$PSA_cost_biop)
      cost_US       <- cost_US_base       * (1 + risk_data$PSA_cost_US)
      cost_MRI      <- cost_MRI_base      * (1 + risk_data$PSA_cost_MRI)
      cost_drug     <- cost_drug_base     * (1 + risk_data$PSA_cost_drug)
    }

    if (PSA == 1L) {
      # ... existing scalar extraction ...
      psa_cols <- grep("^PSA_", names(risksample), value = TRUE)
      risksample <- risksample[, !names(risksample) %in% psa_cols]
    }
    
    # ---- Run DES ------------------------------------------------------------
    run_des_vectorised(
      risk_df            = risksample,
      screen_strategy    = screen_strategy,
      start_age          = start_age,
      screen_startage    = screen_startage,
      discount_cost      = discount_cost,
      recall_rate        = recall_rate,
      Vm                 = Vm,
      Vc                 = Vc,
      cost_screen        = cost_screen,
      cost_strat         = cost_strat,
      cost_US            = cost_US,
      cost_MRI           = cost_MRI,
      cost_follow_up     = cost_follow_up,
      cost_biop          = cost_biop,
      cost_DCIS          = cost_DCIS,
      cost_drug          = cost_drug,
      biopsy_rate        = biopsy_rate,
      course_length      = course_length,
      PREVENTATIVE_DRUG  = PREVENTATIVE_DRUG,
      PSA                = PSA,
      cost_in_full_courses = cost_in_full_courses,
      uptake             = uptake,
      persistence        = persistence
    )

  } # end ii loop

  negsamplefn(screen_strategy, MISCLASS, PSA)
  message(paste("Strategy", r, "complete"))

} # end foreach strategy loop

stopCluster(cl)
toc()
