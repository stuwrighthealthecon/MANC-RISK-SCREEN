create_sample <- function(PSA = 0, intervals = 0, seed = 1, screen_strategy) {
  #Remove existing samples to avoid errors
  do.call(
    file.remove,
    list(list.files("Risksamplewithmisclass/", full.names = TRUE))
  )
  do.call(file.remove, list(list.files("Risksample/", full.names = TRUE)))

  #Import synthetic dataset derived from PROCAS2 study
  risk_mat <- read.csv("Data/synthetic_risk_data.csv")[, 2:4]
  colnames(risk_mat) <- c("VBD", "tenyrrisk", "liferisk")

  #Creat column for risk group to be entered
  risk_mat$risk_group <- numeric(length(risk_mat$VBD))

  #Set VDG based on breast density
  risk_mat$VDG <- 1 + findInterval(risk_mat[, 1], VDG_interval)

  #Breast density cut-offs for supplemental screening
  density_cutoff <- 3

  #Set up data frame of women's lifetimes to simulate
  risksample <- risk_mat[
    sample(nrow(risk_mat), inum * ifelse(PSA == 1, mcruns, 1), replace = TRUE),
  ]
  risksample[, c(
    "MRI_screen",
    "US_screen",
    "risk_predicted",
    "feedback",
    "interval_change",
    "life_expectancy",
    "cancer",
    "clinical_detect_size",
    "growth_rate"
  )] <- numeric(length = length(risksample$VBD))

  #If risk-stratified screening used then determine if each woman chooses to have
  #risk predicted, attends risk consultation, and changes interval

  risksample$risk_predicted <- rbinom(length(nrow(risksample)), 1, risk_uptake)
  risksample$feedback <- ifelse(
    risksample$risk_predicted == 1 &
      rbinom(length(nrow(risksample)), 1, (risk_feedback)) == 1,
    1,
    0
  )
  risksample$interval_change <- ifelse(
    risksample$feedback == 1 &
      rbinom(length(nrow(risksample)), 1, risk_feedback) == 1,
    1,
    0
  )

  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy <- rweibull(
    n = length(risksample$life_expectancy),
    shape = acmmortality_wb_a,
    scale = acmmortality_wb_b
  )
  risksample$life_expectancy <- ifelse(
    risksample$life_expectancy <= start_age,
    qweibull(
      p = dqrunif(
        n = 1,
        min = pweibull(
          q = start_age,
          shape = acmmortality_wb_a,
          scale = acmmortality_wb_b
        ),
        max = 1
      ),
      shape = acmmortality_wb_a,
      scale = acmmortality_wb_b
    ),
    risksample$life_expectancy
  )
  risksample$life_expectancy <- ifelse(
    risksample$life_expectancy >=
      rep(100, length = length(risksample$life_expectancy)),
    99.99,
    risksample$life_expectancy
  )

  #Determine if a cancer will develop
  risksample$cancer <- ifelse(
    dqrunif(length(risksample$cancer), 0, 1) < (risksample$liferisk / 100),
    1,
    0
  )

  #Set clinical detection size for cancer
  risksample$clinical_detect_size <- risksample$cancer *
    (dqrnorm(
      n = length(risksample$cancer),
      mean = clin_detection_m,
      sd = clin_detection_sd
    ))

  #Set growth rate
  risksample$growth_rate <- risksample$cancer *
    qlnorm(
      dqrunif(length(risksample$cancer), 0, 1),
      meanlog = log_norm_mean,
      sdlog = sqrt(log_norm_sd)
    )

  if (PSA == 1) {
    if (intervals == 0) {
      #########################Add Monte Carlo Draws into Sample##############

      source(file = "psa_params.R")

      #Draw stage I to III survival parameters
      PSA_gamma_survival <- mvrnorm(mcruns, survmeans, survcovmat) %>%
        data.frame()

      # Draw Metatstatic survival parameters
      PSA_meta_survival <- mvrnorm(mcruns, metmeans, metmat) %>%
        data.frame()

      #Draw Mammography with sensitivity conditional on tumour diameter parameters W-F
      PSA_beta1 <- rnorm(mcruns, PSA_beta1[1], PSA_beta1[2])
      PSA_beta2 <- rnorm(mcruns, PSA_beta2[1], PSA_beta2[2])

      #Draw Mammography sensitivity by volpara density grade from PREVENTICON
      PSA_Sen_VDG <- data.frame(
        rbeta(mcruns, PSA_Sen_VDG1[1], PSA_Sen_VDG1[2]),
        rbeta(mcruns, PSA_Sen_VDG2[1], PSA_Sen_VDG2[2]),
        rbeta(mcruns, PSA_Sen_VDG3[1], PSA_Sen_VDG3[2]),
        rbeta(mcruns, PSA_Sen_VDG4[1], PSA_Sen_VDG4[2])
      )

      #Draw supplemental Screening CDRs
      PSA_MRI_cdr <- rbeta(mcruns, PSA_MRI_cdr[1], PSA_MRI_cdr[2]) #CDR for MRI in Mammo negative women (incremental)
      PSA_US_cdr <- rbeta(mcruns, PSA_US_cdr[1], PSA_US_cdr[2]) #CDR for US in Mammo negative women (incremental)

      #Draw tumour growth rate parameters
      PSA_log_norm_mean <- rnorm(
        mcruns,
        PSA_log_norm_mean[1],
        PSA_log_norm_mean[2]
      )
      PSA_log_norm_sd <- rnorm(mcruns, PSA_log_norm_sd[1], PSA_log_norm_sd[2])

      #Chemoprevention drug parameters

      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame()

      #Drug uptake parameters

      PSA_uptake_1 <- rnorm(mcruns, PSA_uptake_1[1], PSA_uptake_1[2])
      PSA_uptake_2 <- rnorm(mcruns, PSA_uptake_2[1], PSA_uptake_2[2])

      #Draw costs
      PSA_cost_strat <- (rlnorm(mcruns, PSA_cost_strat[1], PSA_cost_strat[2]))
      PSA_costvar <- rnorm(mcruns, PSA_costvar[1], PSA_costvar[2])
      PSA_costscreen <- rnorm(mcruns, PSA_costscreen[1], PSA_costscreen[2])
      PSA_cost_follow_up <- rnorm(
        mcruns,
        PSA_cost_follow_up[1],
        PSA_cost_follow_up[2]
      )
      PSA_cost_biop <- rnorm(mcruns, PSA_cost_biop[1], PSA_cost_biop[2])
      PSA_cost_US <- rnorm(mcruns, PSA_cost_US[1], PSA_cost_US[2])
      PSA_cost_MRI <- rnorm(mcruns, PSA_cost_MRI[1], PSA_cost_MRI[2])
      PSA_cost_drug <- rnorm(mcruns, PSA_cost_drug[1], PSA_cost_drug[2])

      #Generate utility values
      PSA_util <- (1 - exp(mvrnorm(mcruns, utilmeans, covutil))) %>%
        data.frame()
    } else {
      source(file = "psa_params_intervals.R")

      PSA_gamma_survival <- data.frame(
        runif(mcruns, PSA_gamma_survival1[1], PSA_gamma_survival1[2]),
        runif(mcruns, PSA_gamma_survival2[1], PSA_gamma_survival2[2]),
        runif(mcruns, PSA_gamma_survival3[1], PSA_gamma_survival3[2])
      )
      PSA_meta_survival <- data.frame(
        runif(mcruns, PSA_meta_survival1[1], PSA_meta_survival1[2]),
        runif(mcruns, PSA_meta_survival2[1], PSA_meta_survival2[2]),
        runif(mcruns, PSA_meta_survival3[1], PSA_meta_survival3[2])
      )

      #Mammography with sensitivity conditional on tumour diameter parameters W-F
      PSA_beta1 <- rnorm(mcruns, PSA_beta1[1], PSA_beta1[2])
      PSA_beta2 <- rnorm(mcruns, PSA_beta2[1], PSA_beta2[2])

      #Mammography sensitivity by volpara density grade from PREVENTICON
      PSA_Sen_VDG <- data.frame(
        rbeta(mcruns, PSA_Sen_VDG1[1], PSA_Sen_VDG1[2]),
        rbeta(mcruns, PSA_Sen_VDG2[1], PSA_Sen_VDG2[2]),
        rbeta(mcruns, PSA_Sen_VDG3[1], PSA_Sen_VDG3[2]),
        rbeta(mcruns, PSA_Sen_VDG4[1], PSA_Sen_VDG4[2])
      )

      #Supplemental Screening
      PSA_MRI_cdr <- rbeta(mcruns, PSA_MRI_cdr[1], PSA_MRI_cdr[2]) #CDR for MRI in Mammo negative women (incremental)
      PSA_US_cdr <- rbeta(mcruns, PSA_US_cdr[1], PSA_US_cdr[2]) #CDR for US in Mammo negative women (incremental)

      #Draw tumour growth rate parameters
      PSA_log_norm_mean <- rnorm(
        mcruns,
        PSA_log_norm_mean[1],
        PSA_log_norm_mean[2]
      )
      PSA_log_norm_sd <- rnorm(mcruns, PSA_log_norm_sd[1], PSA_log_norm_sd[2])

      #Drug parameters

      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame()

      PSA_uptake_1 <- rnorm(mcruns, PSA_uptake_1[1], PSA_uptake_1[2])
      PSA_uptake_2 <- rnorm(mcruns, PSA_uptake_2[1], PSA_uptake_2[2])

      #Costs
      #Draw costs
      PSA_cost_strat <- (rlnorm(mcruns, PSA_cost_strat[1], PSA_cost_strat[2]))
      PSA_costvar <- rnorm(mcruns, PSA_costvar[1], PSA_costvar[2])
      PSA_costscreen <- rnorm(mcruns, PSA_costscreen[1], PSA_costscreen[2])
      PSA_cost_follow_up <- rnorm(
        mcruns,
        PSA_cost_follow_up[1],
        PSA_cost_follow_up[2]
      )
      PSA_cost_biop <- rnorm(mcruns, PSA_cost_biop[1], PSA_cost_biop[2])
      PSA_cost_US <- rnorm(mcruns, PSA_cost_US[1], PSA_cost_US[2])
      PSA_cost_MRI <- rnorm(mcruns, PSA_cost_MRI[1], PSA_cost_MRI[2])
      PSA_cost_drug <- rnorm(mcruns, PSA_cost_drug[1], PSA_cost_drug[2])

      PSA_util <- data.frame(
        runif(mcruns, PSA_util_1[1], PSA_util_1[2]),
        runif(mcruns, PSA_util_4[1], PSA_util_4[2])
      )
    }

    #Generate id for monte carlo set
    mcid <- c(1:mcruns)

    #Bind monte carlo draws
    PSA_all_p <- cbind(
      PSA_gamma_survival,
      PSA_meta_survival,
      PSA_beta1,
      PSA_beta2,
      PSA_Sen_VDG,
      PSA_MRI_cdr,
      PSA_US_cdr,
      PSA_log_norm_mean,
      PSA_log_norm_sd,
      PSA_eff,
      PSA_dropout,
      PSA_uptake_1,
      PSA_uptake_2,
      PSA_cost_strat,
      PSA_costvar,
      PSA_util,
      PSA_costscreen,
      PSA_cost_follow_up,
      PSA_cost_biop,
      PSA_cost_US,
      PSA_cost_MRI,
      PSA_cost_drug,
      mcid
    )
    PSA_all_p <- as.data.frame(PSA_all_p)
    colnames(PSA_all_p) <- c(
      "PSA_gamma_survival_1",
      "PSA_gamma_survival_2",
      "PSA_gamma_survival_3",
      "PSA_meta_survival_54",
      "PSA_meta_survival_74",
      "PSA_meta_survival_99",
      "PSA_beta_1",
      "PSA_beta_2",
      'PSA_VDG1_sen',
      'PSA_VDG2_sen',
      'PSA_VDG3_sen',
      'PSA_VDG4_sen',
      "PSA_MRI_cdr",
      "PSA_US_cdr",
      "PSA_log_norm_mean",
      "PSA_log_norm_sd",
      "PSA_eff_ana",
      "PSA_eff_tam",
      "PSA_dropout_ana",
      "PSA_dropout_tam",
      "PSA_uptake_1",
      "PSA_uptake_2",
      "PSA_cost_strat",
      "PSA_costvar",
      "PSA_util_1to3",
      "PSA_util_4",
      "PSA_costscreen",
      "PSA_cost_follow_up",
      "PSA_cost_biop",
      "PSA_cost_US",
      "PSA_cost_MRI",
      "PSA_cost_drug",
      "mcid"
    )

    #Bind individual level parameters and monte carlo draws
    masterframe <- data.frame(matrix(
      nrow = inum * mcruns,
      ncol = length(risksample[1, ]) + length(PSA_all_p[1, ])
    ))
    masterframe[, 1:14] <- risksample
    masterframe[, 15:47] <- PSA_all_p
    colnames(masterframe)[1:14] <- colnames(risksample)
    colnames(masterframe)[15:47] <- colnames(PSA_all_p)

    #Split the dataframe into chunks for easier computation
    masterframe$split <- (rep(
      1:chunks,
      times = round(length(masterframe$VBD) / chunks)
    ))
    masterframe <- masterframe %>% filter(masterframe$life_expectancy >= 50)

    negsample <- masterframe %>% filter(masterframe$cancer == 0)
    save(negsample, file = paste("Risksample/negsample.Rdata", sep = ""))
    masterframe <- masterframe %>% filter(masterframe$cancer == 1)

    risksplit <- split(masterframe, masterframe$split)

    #Clean up redundant inputs
    rm(masterframe, risksample, PSA_all_p, risk_mat)
    gc()

    #Save risk sample in chunks
    for (i in 1:chunks) {
      cancer_col <- paste("X", i, ".cancer", sep = "") %>% as.name()
      splitsample <- risksplit[i] %>%
        as.data.frame() %>%
        filter(!!cancer_col == 1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
      save(
        splitsample,
        file = paste("Risksample/possample_", i, ".Rdata", sep = "")
      )
    }
  } else {
    risksample$split <- (rep(
      1:chunks,
      times = round(length(risksample$VBD) / chunks)
    ))
    risksample <- risksample %>% filter(risksample$life_expectancy >= 50)

    negsample <- risksample %>% filter(risksample$cancer == 0)
    save(negsample, file = paste("Risksample/negsample.Rdata", sep = ""))
    risksample <- risksample %>% filter(risksample$cancer == 1)
    risksplit <- split(risksample, risksample$split)
  }

  #Save risk sample in chunks
  for (i in 1:chunks) {
    cancer_col <- paste("X", i, ".cancer", sep = "") %>% as.name()
    splitsample <- risksplit[i] %>%
      as.data.frame() %>%
      filter(!!cancer_col == 1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
    save(
      splitsample,
      file = paste("Risksample/possample_", i, ".Rdata", sep = "")
    )
  }
}
cmp_create_sample <- cmpfun(create_sample)

# Version of create_sample that stores separate "true" and predicted ten year
# risk:
create_sample_with_misclass <- function(
  PSA = 0,
  intervals = 0,
  seed = 1,
  screen_strategy
) {
  #Import synthetic dataset derived from PROCAS2 study
  risk_mat <- read.csv("Data/synthetic_risk_data_with_misclassification.csv")[,
    2:7
  ] %>%
    dplyr::select(!syn.Age)
  colnames(risk_mat) <- c(
    "VBD",
    "tenyrrisk_est",
    "liferisk_est",
    "tenyrrisk_true",
    "liferisk_true"
  )

  #Creat column for risk group to be entered
  risk_mat$risk_group <- numeric(length(risk_mat$VBD))

  #Set VDG based on breast density
  risk_mat$VDG <- 1 + findInterval(risk_mat[, 1], VDG_interval)

  #Breast density cut-offs for supplemental sreening
  density_cutoff <- 3

  #Set up data frame of women's lifetimes to simulate
  risksample <- risk_mat[
    sample(nrow(risk_mat), inum * ifelse(PSA == 1, mcruns, 1), replace = TRUE),
  ]
  risksample[, c(
    "MRI_screen",
    "US_screen",
    "risk_predicted",
    "feedback",
    "interval_change",
    "life_expectancy",
    "cancer",
    "clinical_detect_size",
    "growth_rate"
  )] <- numeric(length = length(risksample$VBD))

  #If risk-stratified screening used then determine if each woman chooses to have
  #risk predicted, attends risk consultation, and changes interval

  risksample$risk_predicted <- rbinom(length(risksample$VBD), 1, risk_uptake)
  risksample$feedback <- ifelse(
    risksample$risk_predicted == 1 &
      rbinom(length(risksample$VBD), 1, (risk_feedback)) == 1,
    1,
    0
  )
  risksample$interval_change <- ifelse(
    risksample$feedback == 1 &
      rbinom(length(risksample$VBD), 1, risk_feedback) == 1,
    1,
    0
  )

  ###Preload incidence, mortality and clinical detection times
  risksample$life_expectancy <- rweibull(
    n = length(risksample$life_expectancy),
    shape = acmmortality_wb_a,
    scale = acmmortality_wb_b
  )
  risksample$life_expectancy <- ifelse(
    risksample$life_expectancy <= start_age,
    qweibull(
      p = dqrunif(
        n = 1,
        min = pweibull(
          q = start_age,
          shape = acmmortality_wb_a,
          scale = acmmortality_wb_b
        ),
        max = 1
      ),
      shape = acmmortality_wb_a,
      scale = acmmortality_wb_b
    ),
    risksample$life_expectancy
  )
  risksample$life_expectancy <- ifelse(
    risksample$life_expectancy >=
      rep(100, length = length(risksample$life_expectancy)),
    100,
    risksample$life_expectancy
  )
  #Determine if a cancer will develop
  risksample$cancer <- ifelse(
    dqrunif(length(risksample$cancer), 0, 1) < (risksample$liferisk_true / 100),
    1,
    0
  )

  #Set clinical detection size for cancer
  risksample$clinical_detect_size <- risksample$cancer *
    (dqrnorm(
      n = length(risksample$cancer),
      mean = clin_detection_m,
      sd = clin_detection_sd
    ))

  #Set growth rate
  risksample$growth_rate <- risksample$cancer *
    qlnorm(
      dqrunif(length(risksample$cancer), 0, 1),
      meanlog = log_norm_mean,
      sdlog = sqrt(log_norm_sd)
    )

  if (PSA == 1) {
    if (intervals == 0) {
      #########################Add Monte Carlo Draws into Sample##############

      source(file = "psa_params.R")

      #Draw stage I to III survival parameters
      PSA_gamma_survival <- mvrnorm(mcruns, survmeans, survcovmat) %>%
        data.frame()

      # Draw Metatstatic survival parameters
      PSA_meta_survival <- mvrnorm(mcruns, metmeans, metmat) %>%
        data.frame()

      #Draw Mammography with sensitivity conditional on tumour diameter parameters W-F
      PSA_beta1 <- rnorm(mcruns, PSA_beta1[1], PSA_beta1[2])
      PSA_beta2 <- rnorm(mcruns, PSA_beta2[1], PSA_beta2[2])

      #Draw Mammography sensitivity by volpara density grade from PREVENTICON
      PSA_Sen_VDG <- data.frame(
        rbeta(mcruns, PSA_Sen_VDG1[1], PSA_Sen_VDG1[2]),
        rbeta(mcruns, PSA_Sen_VDG2[1], PSA_Sen_VDG2[2]),
        rbeta(mcruns, PSA_Sen_VDG3[1], PSA_Sen_VDG3[2]),
        rbeta(mcruns, PSA_Sen_VDG4[1], PSA_Sen_VDG4[2])
      )

      #Draw supplemental Screening CDRs
      PSA_MRI_cdr <- rbeta(mcruns, PSA_MRI_cdr[1], PSA_MRI_cdr[2]) #CDR for MRI in Mammo negative women (incremental)
      PSA_US_cdr <- rbeta(mcruns, PSA_US_cdr[1], PSA_US_cdr[2]) #CDR for US in Mammo negative women (incremental)

      #Draw tumour growth rate parameters
      PSA_log_norm_mean <- rnorm(
        mcruns,
        PSA_log_norm_mean[1],
        PSA_log_norm_mean[2]
      )
      PSA_log_norm_sd <- rnorm(mcruns, PSA_log_norm_sd[1], PSA_log_norm_sd[2])

      #Chemoprevention drug parameters

      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame()

      #Drug uptake parameters

      PSA_uptake_1 <- rnorm(mcruns, PSA_uptake_1[1], PSA_uptake_1[2])
      PSA_uptake_2 <- rnorm(mcruns, PSA_uptake_2[1], PSA_uptake_2[2])

      #Draw costs
      PSA_cost_strat <- (rlnorm(mcruns, PSA_cost_strat[1], PSA_cost_strat[2]))
      PSA_costvar <- rnorm(mcruns, PSA_costvar[1], PSA_costvar[2])
      PSA_costscreen <- rnorm(mcruns, PSA_costscreen[1], PSA_costscreen[2])
      PSA_cost_follow_up <- rnorm(
        mcruns,
        PSA_cost_follow_up[1],
        PSA_cost_follow_up[2]
      )
      PSA_cost_biop <- rnorm(mcruns, PSA_cost_biop[1], PSA_cost_biop[2])
      PSA_cost_US <- rnorm(mcruns, PSA_cost_US[1], PSA_cost_US[2])
      PSA_cost_MRI <- rnorm(mcruns, PSA_cost_MRI[1], PSA_cost_MRI[2])
      PSA_cost_drug <- rnorm(mcruns, PSA_cost_drug[1], PSA_cost_drug[2])

      #Generate utility values
      PSA_util <- (1 - exp(mvrnorm(mcruns, utilmeans, covutil))) %>%
        data.frame()
    } else {
      source(file = "psa_params_intervals.R")

      PSA_gamma_survival <- data.frame(
        runif(mcruns, PSA_gamma_survival1[1], PSA_gamma_survival1[2]),
        runif(mcruns, PSA_gamma_survival2[1], PSA_gamma_survival2[2]),
        runif(mcruns, PSA_gamma_survival3[1], PSA_gamma_survival3[2])
      )
      PSA_meta_survival <- data.frame(
        runif(mcruns, PSA_meta_survival1[1], PSA_meta_survival1[2]),
        runif(mcruns, PSA_meta_survival2[1], PSA_meta_survival2[2]),
        runif(mcruns, PSA_meta_survival3[1], PSA_meta_survival3[2])
      )

      #Mammography with sensitivity conditional on tumour diameter parameters W-F
      PSA_beta1 <- rnorm(mcruns, PSA_beta1[1], PSA_beta1[2])
      PSA_beta2 <- rnorm(mcruns, PSA_beta2[1], PSA_beta2[2])

      #Mammography sensitivity by volpara density grade from PREVENTICON
      PSA_Sen_VDG <- data.frame(
        rbeta(mcruns, PSA_Sen_VDG1[1], PSA_Sen_VDG1[2]),
        rbeta(mcruns, PSA_Sen_VDG2[1], PSA_Sen_VDG2[2]),
        rbeta(mcruns, PSA_Sen_VDG3[1], PSA_Sen_VDG3[2]),
        rbeta(mcruns, PSA_Sen_VDG4[1], PSA_Sen_VDG4[2])
      )

      #Supplemental Screening
      PSA_MRI_cdr <- rbeta(mcruns, PSA_MRI_cdr[1], PSA_MRI_cdr[2]) #CDR for MRI in Mammo negative women (incremental)
      PSA_US_cdr <- rbeta(mcruns, PSA_US_cdr[1], PSA_US_cdr[2]) #CDR for US in Mammo negative women (incremental)

      #Draw tumour growth rate parameters
      PSA_log_norm_mean <- rnorm(
        mcruns,
        PSA_log_norm_mean[1],
        PSA_log_norm_mean[2]
      )
      PSA_log_norm_sd <- rnorm(mcruns, PSA_log_norm_sd[1], PSA_log_norm_sd[2])

      #Drug parameters

      PSA_eff <- mvrnorm(mcruns, efficacy_mu, efficacy_sigma) %>%
        data.frame()
      PSA_dropout <- mvrnorm(mcruns, dropout_mu, dropout_sigma) %>%
        data.frame()

      PSA_uptake_1 <- rnorm(mcruns, PSA_uptake_1[1], PSA_uptake_1[2])
      PSA_uptake_2 <- rnorm(mcruns, PSA_uptake_2[1], PSA_uptake_2[2])

      #Costs
      #Draw costs
      PSA_cost_strat <- (rlnorm(mcruns, PSA_cost_strat[1], PSA_cost_strat[2]))
      PSA_costvar <- rnorm(mcruns, PSA_costvar[1], PSA_costvar[2])
      PSA_costscreen <- rnorm(mcruns, PSA_costscreen[1], PSA_costscreen[2])
      PSA_cost_follow_up <- rnorm(
        mcruns,
        PSA_cost_follow_up[1],
        PSA_cost_follow_up[2]
      )
      PSA_cost_biop <- rnorm(mcruns, PSA_cost_biop[1], PSA_cost_biop[2])
      PSA_cost_US <- rnorm(mcruns, PSA_cost_US[1], PSA_cost_US[2])
      PSA_cost_MRI <- rnorm(mcruns, PSA_cost_MRI[1], PSA_cost_MRI[2])
      PSA_cost_drug <- rnorm(mcruns, PSA_cost_drug[1], PSA_cost_drug[2])

      PSA_util <- data.frame(
        runif(mcruns, PSA_util_1[1], PSA_util_1[2]),
        runif(mcruns, PSA_util_4[1], PSA_util_4[2])
      )
    }

    #Generate id for monte carlo set
    mcid <- c(1:mcruns)

    #Bind monte carlo draws
    PSA_all_p <- cbind(
      PSA_gamma_survival,
      PSA_meta_survival,
      PSA_beta1,
      PSA_beta2,
      PSA_Sen_VDG,
      PSA_MRI_cdr,
      PSA_US_cdr,
      PSA_log_norm_mean,
      PSA_log_norm_sd,
      PSA_eff,
      PSA_dropout,
      PSA_uptake_1,
      PSA_uptake_2,
      PSA_cost_strat,
      PSA_costvar,
      PSA_util,
      PSA_costscreen,
      PSA_cost_follow_up,
      PSA_cost_biop,
      PSA_cost_US,
      PSA_cost_MRI,
      PSA_cost_drug,
      mcid
    )
    PSA_all_p <- as.data.frame(PSA_all_p)
    colnames(PSA_all_p) <- c(
      "PSA_gamma_survival_1",
      "PSA_gamma_survival_2",
      "PSA_gamma_survival_3",
      "PSA_meta_survival_54",
      "PSA_meta_survival_74",
      "PSA_meta_survival_99",
      "PSA_beta_1",
      "PSA_beta_2",
      'PSA_VDG1_sen',
      'PSA_VDG2_sen',
      'PSA_VDG3_sen',
      'PSA_VDG4_sen',
      "PSA_MRI_cdr",
      "PSA_US_cdr",
      "PSA_log_norm_mean",
      "PSA_log_norm_sd",
      "PSA_eff_ana",
      "PSA_eff_tam",
      "PSA_dropout_ana",
      "PSA_dropout_tam",
      "PSA_uptake_1",
      "PSA_uptake_2",
      "PSA_cost_strat",
      "PSA_costvar",
      "PSA_util_1to3",
      "PSA_util_4",
      "PSA_costscreen",
      "PSA_cost_follow_up",
      "PSA_cost_biop",
      "PSA_cost_US",
      "PSA_cost_MRI",
      "PSA_cost_drug",
      "mcid"
    )

    #Bind individual level parameters and monte carlo draws
    masterframe <- data.frame(matrix(
      nrow = inum * mcruns,
      ncol = length(risksample[1, ]) + length(PSA_all_p[1, ])
    ))
    masterframe[, 1:16] <- risksample
    masterframe[, 17:49] <- PSA_all_p
    colnames(masterframe)[1:16] <- colnames(risksample)
    colnames(masterframe)[17:49] <- colnames(PSA_all_p)

    #Split the dataframe into chunks for easier computation
    masterframe$split <- (rep(
      1:chunks,
      times = round(length(masterframe$VBD) / chunks)
    ))
    masterframe <- masterframe %>% filter(masterframe$life_expectancy >= 50)

    negsample <- masterframe %>% filter(masterframe$cancer == 0)
    save(
      negsample,
      file = paste("Risksamplewithmisclass/negsample.Rdata", sep = "")
    )
    risksplit <- split(masterframe, masterframe$split)

    #Clean up redundant inputs
    rm(masterframe, risksample, PSA_all_p, risk_mat)
    gc()

    #Save risk sample in chunks
    for (i in 1:chunks) {
      cancer_col <- paste("X", i, ".cancer", sep = "") %>% as.name()
      splitsample <- risksplit[i] %>%
        as.data.frame() %>%
        filter(!!cancer_col == 1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
      save(
        splitsample,
        file = paste("Risksamplewithmisclass/possample_", i, ".Rdata", sep = "")
      )
    }
  } else {
    risksample$split <- (rep(
      1:chunks,
      times = round(length(risksample$VBD) / chunks)
    ))
    risksample <- risksample %>% filter(risksample$life_expectancy >= 50)

    negsample <- risksample %>% filter(risksample$cancer == 0)
    save(
      negsample,
      file = paste("Risksamplewithmisclass/negsample.Rdata", sep = "")
    )
    risksplit <- split(risksample, risksample$split)

    #Save risk sample in chunks
    for (i in 1:chunks) {
      cancer_col <- paste("X", i, ".cancer", sep = "") %>% as.name()
      splitsample <- risksplit[i] %>%
        as.data.frame() %>%
        filter(!!cancer_col == 1) # We give this the same name as the merged sample to avoid extraneous if statements in the simulation script
      save(
        splitsample,
        file = paste("Risksamplewithmisclass/possample_", i, ".Rdata", sep = "")
      )
    }
  }
}

# Work out new version of IncidenceMortality adjusted for effect of drug on
# hazard ratios for a single individual.
get_drug_adj_IM <- function(
  ind_from_risksample,
  uptake,
  persistence,
  risk_red
) {
  drug_IM <- Incidence_Mortality
  # If takes drug assign time taking, otherwise set to zero
  if (
    dqrunif(1, 0, 1) <
      uptake[
        ind_from_risksample$risk_group,
        ind_from_risksample$starting_menses_status
      ]
  ) {
    time_taking_drug <- min(
      rexp(
        1,
        rate = persistence[
          ind_from_risksample$risk_group,
          ind_from_risksample$starting_menses_status
        ]
      ),
      course_length
    )
    hazard_ratio <- (1 -
      risk_red[
        ind_from_risksample$risk_group,
        ind_from_risksample$starting_menses_status
      ]) *
      time_taking_drug *
      log(1 / completion_prob[ind_from_risksample$starting_menses_status])
    new_weibull_scale <- inc_scale / hazard_ratio

    drug_IM$BC_age <- dweibull(
      Incidence_Mortality$age,
      shape = inc_shape,
      scale = new_weibull_scale
    )
  } else {
    time_taking_drug <- 0
  }

  return(list(time_taking_drug, drug_IM))
}

# The following function is used during PSA to get new risk reduction and
# dropout rate matrices for the preventative drug courses based on log hazard
# ratios drawn during the sample generation.
redraw_drug_pars <- function(risksample) {
  # New efficacies
  ana_eff <- exp(risksample$PSA_eff_ana)
  tam_eff <- exp(risksample$PSA_eff_tam)

  full_course_len <- 5

  # Assume constant drop out rate with 77% of individuals reaching 5yr mark
  tam_completion_prob <- .77

  # New dropout estimates
  ana_dropout_rate <- risksample$PSA_dropout_ana
  tam_dropout_rate <- risksample$PSA_dropout_tam

  # Estimate Anastrozole completion probability from Tamoxifen estimate and hazard ratios:
  ana_completion_prob <- exp(ana_dropout_rate - tam_dropout_rate) *
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
    risk_red <- matrix(
      c(ana_eff, tam_eff, ana_eff, tam_eff),
      nrow = 2,
      ncol = 2
    )

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
  return(list(risk_red, uptake, persistence))
}
