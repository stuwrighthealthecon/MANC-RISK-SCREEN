#Load require packages
#-----------------------------------------------------------------------------
library(doParallel)
library(MASS)
library(dqrng)
library(compiler)
library(tidyverse)
library(iterators)
library(here)
library(purrr)        
library(tictoc)
library(matrixStats)

# ----------------------------------------------------------------------------
#Set model controls
#-----------------------------------------------------------------------------
controls <- list(
  strategies       = c(0,1,2,3,4,9), #vector of strategies to simulate
  gensample        = TRUE, #create new sample to simulate?
  MISCLASS         = TRUE, #should error in risk prediction be included?
  PREVENTATIVE_DRUG = FALSE, #should risk reducing medicines be used for high risk?
  supplemental_screening = FALSE, #should ultrasound and MRI be used as supplemental screening?
  PSA              = TRUE, #run PSA?
  intervals        = FALSE, #run PSA with wide distributions for GAM estimation?
  desired_cases    = 100, #number of cancer cases required
  mcruns           = 100, #Number of Monte Carlo runs
  seed             = 42, #Set seed for random number generation
  n_cores          = max(1L, parallel::detectCores() - 1L) #Select number of computer cores to use
)

MISCLASS          <- controls$MISCLASS
PREVENTATIVE_DRUG <- controls$PREVENTATIVE_DRUG
PSA               <- as.integer(controls$PSA)
intervals         <- as.integer(controls$intervals)

#Place control on number of cores for memory management at large sample sizes
controls$n_cores <- if (controls$desired_cases >= 200000) 3L else
  if (controls$desired_cases >= 100000) 6L else
    max(1L, parallel::detectCores() - 1L)

# -----------------------------------------------------------------------------
# Set output directories
# -----------------------------------------------------------------------------
det_output_path <- "Deterministic results/"
psa_output_path <- "PSA results/"

if (MISCLASS & PREVENTATIVE_DRUG) {
  dir.create("Deterministic results/misclassification_and_preventative_drug",
             showWarnings = FALSE, recursive = TRUE)
  dir.create("PSA results/misclassification_and_preventative_drug",
             showWarnings = FALSE, recursive = TRUE)
  det_output_path <- "Deterministic results/misclassification_and_preventative_drug/"
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
} else if (MISCLASS) {
  dir.create("Deterministic results/misclassification", showWarnings = FALSE, recursive = TRUE)
  dir.create("PSA results/misclassification",           showWarnings = FALSE, recursive = TRUE)
  det_output_path <- "Deterministic results/misclassification/"
  psa_output_path <- "PSA results/misclassification/"
} else if (PREVENTATIVE_DRUG) {
  dir.create("Deterministic results/preventative_drug",             showWarnings = FALSE, recursive = TRUE)
  dir.create("PSA results/misclassification_and_preventative_drug", showWarnings = FALSE, recursive = TRUE)
  det_output_path <- "Deterministic results/preventative_drug/"
  psa_output_path <- "PSA results/misclassification_and_preventative_drug/"
}

sample_fname <- "possample"

# -----------------------------------------------------------------------------
# Set parameters and functions
# -----------------------------------------------------------------------------
expected_prev <- 0.12
desired_cases <- controls$desired_cases
mcruns        <- controls$mcruns

#Set total number of women to simulate to get desired case numbers
inum <- ceiling(desired_cases / expected_prev)

#Set seed for dq based random number generators
dqset.seed(controls$seed)

# Source function files
source("MANC_RISK_SCREEN_functions.R")
source("risksample function.R")
source("negsample function.R")
source("des_vectorised_df.R")

# Source fixed parameters
source("params.R")

# -----------------------------------------------------------------------------
# Generate sample
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

# Load sample
if (MISCLASS) {
  load(paste0("Risksamplewithmisclass/", sample_fname, ".Rdata"))
} else {
  load(paste0("Risksample/", sample_fname, ".Rdata"))
}

#Fix column names for risksample
prefix <- paste0("^", "X", 1, ".")
names(risksample) <- sub(prefix, "", names(risksample))

# Risk group assignment
if (MISCLASS) {
  risk_col <- "tenyrrisk_est"
} else {
  risk_col <- "tenyrrisk"
}

# Initialise risk_group — will be overwritten per-strategy inside loop
# for strategies that use risk stratification; set to 0 (average) as default
risksample$risk_group <- 0L

#Assign supplemental screening if used
if (controls$supplemental_screening) {
  supp_col        <- ifelse(MISCLASS, "tenyrrisk_true", "tenyrrisk_est")
  high_risk_dense <- risksample$VDG >= density_cutoff &
                     risksample[[supp_col]] >= 8
  low_risk_dense  <- risksample$VDG >= density_cutoff &
                     risksample[[supp_col]] < 8
  risksample$MRI_screen[high_risk_dense] <- 1L
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

#Drop un-needed columns for memory management, some now some later
if (MISCLASS) {
  drop_cols_now  <- c("VBD", "cancer", "feedback", "liferisk_est", "liferisk_true")
  drop_cols_late <- c("tenyrrisk_est", "tenyrrisk_true")
} else {
  drop_cols_now  <- c("VBD", "cancer", "feedback", "liferisk")
  drop_cols_late <- c("tenyrrisk")
}
splitmaster_base <- risksample_master[, !names(risksample_master) %in% drop_cols_now]

if (PSA == 0L) {
  psa_cols     <- grep("^PSA_", names(splitmaster_base), value = TRUE)
  splitmaster_base <- splitmaster_base[, !names(splitmaster_base) %in% psa_cols]
}
rm(risksample_master); gc()

# -----------------------------------------------------------------------------
# Main strategy loop — parallelised via foreach
# -----------------------------------------------------------------------------
tic() #Start timer

cl <- makeCluster(controls$n_cores, outfile = "")
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

  # Pre-split data.frame before ii loop (for memory management)
  risk_group_list <- split(splitmaster, splitmaster$risk_group)
  rm(splitmaster); gc()

  # ---- ii loop: one iteration per risk group --------------------------------
  for (ii in seq_along(risk_groups)) {
    tryCatch({
      risksample <- risk_group_list[[as.character(risk_groups[ii])]]

    # Catch NA risk_group before it causes errors
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

      #Create column of time taking drugs for those who take them
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

      prob_matrix <- outer(
        risksample$weibullrisk, ages,
        FUN = function(scale, age) dweibull(age, shape = inc_shape, scale = scale)
      )

      age_mask    <- outer(risksample$life_expectancy, ages, FUN = ">")
      prob_matrix <- prob_matrix * age_mask

      row_sums    <- rowSums(prob_matrix)
      prob_matrix <- prob_matrix / row_sums

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

    }, error = function(e) {
      message(sprintf(
        "FAILED: strategy=%s ii=%d risk_group=%s nrow=%d PSA=%d mcruns=%d — %s",
        screen_strategy, ii, risk_groups[ii],
        nrow(risk_group_list[[as.character(risk_groups[ii])]]),
        PSA, mcruns, conditionMessage(e)
      ))
      stop(e)
    })
  } # end ii loop

  #Run the model for people without cancer
  negsamplefn(screen_strategy, MISCLASS, PSA)
  
  message(paste("Strategy", r, "complete"))

} # end foreach strategy loop

stopCluster(cl)
toc()
