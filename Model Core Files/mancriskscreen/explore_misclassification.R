# In this script we show how to get beta distributed risk distributions for
# "true" risk based on comparing observed and estimated risk from Tyrer-Cuzick,
# and use Gaussian processes to estimate corrected life risks from these
# corrected 10-year risks.

EXPLORE_GPS <- FALSE # Set to true to look at a few alternative targets for fitting in the Gaussian process section

PRINT_RMSES <- FALSE # Set to true to print each RMSE calculated during the beta distribution calibration stage - this can be useful for doing diagnostics

SAVE_OUTPUTS <- FALSE # Set to true to overwrite existing sample with misclassification

library(dplyr)
library(GauPro)
library(ggplot2)
library(magrittr)
library(pdftools)
library(purrr)
library(tidyr)

# Load in risk table from Brentnall 2018 supplement:
supp_data <- pdf_data("coi180011supp1_prod.pdf")[[5]][-c(1:15), ] %>%
  mutate(
    x = round(x / 3), #reduce resolution to minimise inconsistent coordinates
    y = round(y / 3)
  ) %>%
  arrange(y, x) %>% #sort in reading order
  mutate(
    group = cumsum(
      !sapply(1:length(x), FUN = function(i) {
        # Gather nearby strings into same column (cutoff of 3 chosen by trial and error)
        if (i == 1) {
          return(TRUE)
        } else {
          return(
            (y[i - 1] == y[i]) & ((x[i] - x[i - 1]) <= 3)
          )
        }
      })
    )
  ) %>% #identify text with spaces and paste
  group_by(group) %>%
  summarise(x = first(x), y = first(y), text = paste(text, collapse = " ")) %>%
  filter(!grepl('95%', text)) %>% # Remove column headings for 95% confidence intervals
  filter(!grepl('%)', text)) %>% # Remove percentage conversions for numbers
  group_by(y) %>%
  mutate(colno = row_number()) %>% #add column numbers for table data
  ungroup() %>%
  dplyr::select(text, colno, y) %>%
  pivot_wider(names_from = colno, values_from = text) %>% #pivot into table format
  dplyr::select(-y) %>%
  set_colnames(c(
    "Risk_class",
    "N",
    "FU",
    "O",
    "E",
    "O_E",
    "O_E 95% CI",
    "IR",
    "IRR",
    "IRR 95% CI"
  )) %>%
  dplyr::select(1:10) %>%
  slice(-c(1, 44, 45, 46))
model_idxs <- grep("Ty", supp_data$Risk_class)
supp_data$Model_age <- sapply(1:nrow(supp_data), FUN = function(i) {
  supp_data[model_idxs[ceiling(i / 7)], ] %>% # Add column specifying age group and risk class
    as.character() %>%
    discard(is.na) %>%
    paste(collapse = " ")
})
supp_data <- supp_data[-c(model_idxs), ]

supp_data[, c("N", "E", "O", "O_E")] <- supp_data[, c("N", "E", "O", "O_E")] %>%
  sapply(as.numeric)

supp_data <- supp_data %>%
  filter(grepl("density", Model_age))

#### Do least-squares beta regression systematically for each age/model combination

n_classes <- supp_data$Model_age %>%
  unique() %>%
  length()
beta_fits_df <- data.frame(
  "obs_mean" = numeric(length = n_classes),
  "alpha_e" = numeric(length = n_classes),
  "beta_e" = numeric(length = n_classes),
  "mean_e" = numeric(length = n_classes),
  "alpha_o" = numeric(length = n_classes),
  "beta_o" = numeric(length = n_classes),
  "mean_o" = numeric(length = n_classes),
  row.names = unique(supp_data$Model_age)
)

risk_bds <- c(.02, .03, .05, .08, 1.)
for (m in unique(supp_data$Model_age)) {
  print(m)
  supp_data_m <- supp_data[which(supp_data$Model_age == m), ]
  data_quantiles <- cumsum(supp_data_m$N[2:6]) / supp_data_m$N[1]
  cohort_cumsum <- c(0, cumsum(supp_data_m$N[2:6]))

  cases_df <- supp_data_m %>%
    select(Risk_class, N, O, E)

  weight_calculator <- function(alpha, beta) {
    sapply(1:1000, FUN = function(i) {
      risk_draw <- rbeta(supp_data_m$N[1], alpha, beta) %>% sort()
      sapply(2:6, FUN = function(j) {
        (rbernoulli(
          n = supp_data_m$N[j],
          p = risk_draw[(cohort_cumsum[j - 1] + 1):cohort_cumsum[j]]
        ) %>%
          which() %>%
          length())
      })
    })
  }

  # First do expected
  tenyr_risk <- c(supp_data_m$E[2:5] / supp_data_m$N[2:5], 1.)
  exp_mean <- supp_data_m$E[which(supp_data_m$Risk_class == "Total")] /
    supp_data_m$N[which(supp_data_m$Risk_class == "Total")]

  alpha_0 <- supp_data_m$E[which(supp_data_m$Risk_class == "Total")]
  beta_0 <- supp_data_m$N[which(supp_data_m$Risk_class == "Total")] - alpha_0

  get_cases_match_rmse <- function(alpha) {
    beta <- (alpha / exp_mean) - alpha
    exp_weights <- sapply(1:100, FUN = function(i) {
      risk_draw <- rbeta(supp_data_m$N[1], alpha, beta) %>% sort()
      sapply(2:6, FUN = function(j) {
        (rbernoulli(
          n = supp_data_m$N[j],
          p = risk_draw[(cohort_cumsum[j - 1] + 1):cohort_cumsum[j]]
        ) %>%
          which() %>%
          length())
      })
    }) %>%
      rowMeans()
    rmse <- ((1 / length(exp_weights)) *
      (exp_weights - supp_data_m$E[2:6])^2) %>%
      sum() %>%
      sqrt()
    if (PRINT_RMSES) {
      print(rmse)
    }
    return(rmse)
  }
  par_estim_e <- optimize(get_cases_match_rmse, c(0, alpha_0))
  alpha_e <- par_estim_e$minimum
  beta_e <- (par_estim_e$minimum / exp_mean) - par_estim_e$minimum

  exp_weights <- weight_calculator(alpha_e, beta_e)
  cases_df$E_case_match <- c(sum(rowMeans(exp_weights)), rowMeans(exp_weights))

  # Now do observed cases
  tenyr_risk <- c(supp_data_m$O[2:5] / supp_data_m$N[2:5], 1)
  obs_mean <- supp_data_m$O[which(supp_data_m$Risk_class == "Total")] /
    supp_data_m$N[which(supp_data_m$Risk_class == "Total")]
  get_beta_rmse_case_match_E <- function(alpha) {
    beta <- (alpha / obs_mean) - alpha
    obs_weights <- sapply(1:1000, FUN = function(i) {
      risk_draw <- rbeta(supp_data_m$N[1], alpha, beta) %>% sort()
      sapply(2:6, FUN = function(j) {
        (rbernoulli(
          n = supp_data_m$N[j],
          p = risk_draw[(cohort_cumsum[j - 1] + 1):cohort_cumsum[j]]
        ) %>%
          which() %>%
          length())
      })
    }) %>%
      rowMeans()
    rmse <- ((1 / length(obs_weights)) *
      (obs_weights - supp_data_m$O[2:6])^2) %>%
      sum() %>%
      sqrt()
    if (PRINT_RMSES) {
      print(rmse)
    }
    return(rmse)
  }
  alpha_0 <- (supp_data_m$O[which(supp_data_m$Risk_class == "Total")])
  par_estim_o <- optimize(get_beta_rmse_case_match_E, c(0, alpha_0))
  alpha_o <- par_estim_o$minimum
  beta_o <- (par_estim_o$minimum / obs_mean) - par_estim_o$minimum

  obs_weights <- weight_calculator(alpha_o, beta_o)
  cases_df$O_case_match <- c(sum(rowMeans(obs_weights)), rowMeans(obs_weights))

  cases_df <- cases_df %>% relocate(O, .before = O_case_match)

  if (m == unique(supp_data$Model_age)[1]) {
    all_cases_df <- cases_df
  } else {
    all_cases_df <- cbind(all_cases_df, cases_df)
  }

  beta_fits_df[m, ] <- c(
    alpha_0,
    alpha_e,
    beta_e,
    alpha_e / (alpha_e + beta_e),
    alpha_o,
    beta_o,
    alpha_o / (alpha_o + beta_o)
  )
}

if (SAVE_OUTPUTS) {
  save(beta_fits_df, file = "fitted_risk_distributions_update.csv")
}


#### Next problem: inferring life risk from ten year risk

risk_mat <- read.csv("mancriskscreen/Data/synthetic_risk_data.csv") %>%
  dplyr::select(-"X")

risk_mat$syn.X10yr <- .01 * risk_mat$syn.X10yr # Rescale to numbers rather than percentages
risk_mat$syn.life <- .01 * risk_mat$syn.life

# Inspect visually:

p <- ggplot(risk_mat, aes(x = syn.X10yr, y = syn.life)) +
  geom_point() +
  xlab("Ten year estimated risk") +
  ylab("Lifetime estimated risk")
p

# Try splitting up by age:
risk_mat <- risk_mat %>%
  mutate(
    age_grp = case_when(
      syn.Age < 50 ~ "<50",
      syn.Age >= 50 & syn.Age < 60 ~ "50-60",
      syn.Age >= 60 ~ "60+"
    )
  )
p <- ggplot(risk_mat, aes(x = syn.X10yr, y = syn.life, colour = age_grp)) +
  geom_point(alpha = .2) +
  xlab("Ten year estimated risk") +
  ylab("Lifetime estimated risk") +
  scale_color_discrete()
p

# Should see obvious difference in relationship, so need to include age in regression model

risk_mat$bin.10yr <- risk_mat$syn.X10yr %>%
  cut(breaks = seq(from = 0, to = 1, by = .001))

# Needs to be age stratified, so try first on 50-60 group (main group of interest)

smpl <- risk_mat[which(risk_mat$age_grp == "50-60"), ] %>%
  slice_sample(by = bin.10yr, n = 1)

# Simple version:
x <- smpl$syn.X10yr
y <- smpl$syn.life
kern <- Matern52$new(0)
gp <- GauPro(x, y, kernel = kern)
ypred <- gp$predict(seq(0., .3, .01), se = T) %>%
  mutate(upper = mean + 2 * se, lower = mean - 2 * se, xval = seq(0., .3, .01))
p <- ggplot(smpl, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(data = ypred, aes(x = xval, y = mean), colour = "red") +
  geom_line(data = ypred, aes(x = xval, y = upper), colour = "blue") +
  geom_line(data = ypred, aes(x = xval, y = lower), colour = "blue") +
  labs(title = "GP fit for 50-60y subgroup")
print(p)

# Last point is a bit of an outlier, try retraining without:
smpl_order <- smpl$syn.X10yr %>% order() %>% head(-1) # Order samples by x value and drop last one
x_short <- x[smpl_order]
y_short <- y[smpl_order]
kern <- Matern52$new(0)
gp_short <- GauPro(x_short, y_short, kernel = kern)
ypred <- gp_short$predict(seq(0., .3, .01), se = T) %>%
  mutate(upper = mean + 2 * se, lower = mean - 2 * se, xval = seq(0., .3, .01))
p <- ggplot(smpl, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(data = ypred, aes(x = xval, y = mean), colour = "red") +
  geom_line(data = ypred, aes(x = xval, y = upper), colour = "blue") +
  geom_line(data = ypred, aes(x = xval, y = lower), colour = "blue") +
  labs(title = "GP fit for 50-60y subgroup with final datapoint removed")
print(p)

# Rename to reflect age group
gp_50_60 <- gp_short

# Now do other age bands:

smpl <- risk_mat[which(risk_mat$age_grp == "<50"), ] %>%
  slice_sample(by = bin.10yr, n = 1)

x <- smpl$syn.X10yr
y <- smpl$syn.life
kern <- Matern52$new(0)
gp_u50 <- GauPro(x, y, kernel = kern)
ypred <- gp_u50$predict(seq(0., .3, .01), se = T) %>%
  mutate(upper = mean + 2 * se, lower = mean - 2 * se, xval = seq(0., .3, .01))
p <- ggplot(smpl, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(data = ypred, aes(x = xval, y = mean), colour = "red") +
  geom_line(data = ypred, aes(x = xval, y = upper), colour = "blue") +
  geom_line(data = ypred, aes(x = xval, y = lower), colour = "blue") +
  labs(title = "GP fit for <50y subgroup")
print(p)

# Actually looks like it's points with >25% 10yr risk causing trouble, so this
# might be a better rule to use if we ever do this with new risk data

smpl <- smpl %>% filter(syn.X10yr <= 25)
x <- smpl$syn.X10yr
y <- smpl$syn.life
kern <- Matern52$new(0)
gp_u50 <- GauPro(x, y, kernel = kern)
ypred <- gp_u50$predict(seq(0., .3, .01), se = T) %>%
  mutate(upper = mean + 2 * se, lower = mean - 2 * se, xval = seq(0., .3, .01))
p <- ggplot(smpl, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(data = ypred, aes(x = xval, y = mean), colour = "red") +
  geom_line(data = ypred, aes(x = xval, y = upper), colour = "blue") +
  geom_line(data = ypred, aes(x = xval, y = lower), colour = "blue") +
  labs(title = "GP fit for <50y subgroup with all 10yr risks >25% removed")
print(p)

gp_list <- c(gp_u50, gp_50_60)

#### In this optional section we look at transformations of the data to improve the GP fit.
# Note that the code below this section will still use gp_short, the Gaussian process
# fit to the raw risk data with the final datapoint removed, to generate synthetic
# lifetime risks.

if (EXPLORE_GPS) {
  # Version that estimates log(y-x)

  x <- smpl$syn.X10yr
  y <- (smpl$syn.life - smpl$syn.X10yr) %>% log()
  gp <- GauPro(x, y, kernel = Exponential$new(0))

  # Try plotting risk prediction:

  ypred <- function(z) {
    z + (gp$predict(z) %>% exp())
  }
  yse <- function(z) {
    exp(gp$predict(z) + 0.5 * gp$predict(z, se = T)$se^2) *
      sqrt(exp(gp$predict(z, se = T)$se^2) - 1)
  }
  plot(x, smpl$syn.life)
  curve(ypred(x), add = T, col = 2)
  curve(ypred(x) + 2 * yse(x), add = T, col = 4)
  curve(ypred(x) - 2 * yse(x), add = T, col = 4)

  # Alternative: estimate log(dy/dx)
  smpl_order <- smpl$syn.X10yr %>% order()
  x <- smpl$syn.X10yr[smpl_order][2:length(smpl_order)]
  y <- diff(smpl$syn.life[smpl_order]) / diff(smpl$syn.X10yr[smpl_order])

  x <- x[which(y > 0)]
  y <- y[which(y > 0)]
  gp <- GauPro(x, log(y), kernel = Exponential$new(0))
  plot(x, y)
  yse <- function(z) {
    exp(gp$predict(z) + 0.5 * gp$predict(z, se = T)$se^2) *
      sqrt(exp(gp$predict(z, se = T)$se^2) - 1)
  }
  curve(exp(gp$predict(x)), add = T, col = 2)
  curve(exp(gp$predict(x)) + yse(x), add = T, col = 4)
  curve(exp(gp$predict(x) + 1.1 * min(y)) - yse(x), add = T, col = 4)
}

#### Now add "true" risks to risk_mat, assuming we're working with prediction
# with density.
#
bf_dens_only <- beta_fits_df[grep("density", rownames(beta_fits_df)), ]
supp_data_dens_only <- supp_data %>%
  filter(grepl("density", Model_age)) %>%
  filter(grepl("Total", Risk_class))
cutoffs <- c(0, 50, 60)

# Function to correct 10 year risk by mapping between clibrated beta distributions
true_risk_from_syn <- function(syn_risk, age) {
  age_class <- which(cutoffs < age) %>% max()

  beta_pars_e <- c(
    bf_dens_only$alpha_e[age_class],
    bf_dens_only$beta_e[age_class]
  )
  beta_pars_o <- c(
    bf_dens_only$alpha_o[age_class],
    bf_dens_only$beta_o[age_class]
  )

  # Get individual's position in estimated risk distribution
  quantile <- pbeta(syn_risk, beta_pars_e[1], beta_pars_e[2])

  # Get risk at same position in observed risk distribution:
  true_risk <- qbeta(quantile, beta_pars_o[1], beta_pars_o[2])
  return(true_risk)
}

true_risk_from_syn_v <- Vectorize(true_risk_from_syn)

# Create copy of risk_mat with misclassification
risk_mat_with_mc <- risk_mat %>%
  mutate(syn.X10yr.true = true_risk_from_syn_v(syn.X10yr, syn.Age))
# The next plot should indicate that there's one problematic case where 10 year
# risk gets mapped to 1:
ggplot(risk_mat_with_mc, aes(syn.X10yr, syn.X10yr.true, col = age_grp)) +
  geom_point() +
  geom_abline() +
  labs(title = "10 year risk correction")

# Quick exploration of mapping between risks:

beta_match_df <- data.frame(X = seq(0, 1, .0001)) %>%
  mutate(
    Ye_u50 = pbeta(X, bf_dens_only$alpha_e[1], bf_dens_only$beta_e[1]),
    Yo_u50 = pbeta(X, bf_dens_only$alpha_o[1], bf_dens_only$beta_o[1]),
    Ye_5060 = pbeta(X, bf_dens_only$alpha_e[2], bf_dens_only$beta_e[2]),
    Yo_5060 = pbeta(X, bf_dens_only$alpha_o[2], bf_dens_only$beta_o[2])
  )
beta_match_plt <- ggplot(beta_match_df) +
  geom_line(aes(x = Ye_u50, y = Yo_u50), colour = 'red') +
  geom_line(aes(x = Ye_5060, y = Yo_5060), colour = 'blue') +
  labs(x = "Predicted risk", y = "Corrected risk")
print(beta_match_plt)

# Should see that problematic cases are all ones with high risk
problem_cases <- risk_mat_with_mc %>% filter(syn.X10yr.true == 1)
print(problem_cases$syn.X10yr)

# Try switching to linear interpolation for these given plots look close to linear:
linint_df <- risk_mat_with_mc %>%
  filter(age_grp == "50-60") %>%
  dplyr::select(syn.X10yr, syn.X10yr.true) %>%
  filter(!syn.X10yr.true == 1)
lmodel <- lm(linint_df$syn.X10yr.true ~ linint_df$syn.X10yr)

# Try replacing problem cases with linear extrapolated values:
problem_locs <- which(risk_mat_with_mc$syn.X10yr.true == 1)
risk_mat_with_mc$syn.X10yr.true[problem_locs] <- problem_cases$syn.X10yr *
  lmodel$coefficients[2]

# Should look more reasonable now
ggplot(risk_mat_with_mc, aes(syn.X10yr, syn.X10yr.true, col = age_grp)) +
  geom_point() +
  geom_abline() +
  labs(
    title = "10 year risk correction with outlier handled by linear extrapolation"
  )

# Now add life risks:

get_liferisk <- function(age_grp_idx, syn.X10yr) {
  gp_i <- gp_list[[age_grp_idx]]
  return(gp_i$predict(syn.X10yr))
}

get_liferisk <- Vectorize(get_liferisk)

risk_mat_with_mc <- risk_mat_with_mc %>%
  mutate(age_grp_idx = ifelse(age_grp == "<50", 1, 2)) %>%
  mutate(syn.life.true = get_liferisk(age_grp_idx, syn.X10yr.true))

# Check life risk plots:
ggplot(risk_mat_with_mc, aes(syn.life, syn.life.true, col = age_grp)) +
  geom_point() +
  geom_abline() +
  labs("Lifetime risk correction")
risk_mat_with_mc <- risk_mat_with_mc %>%
  dplyr::select(-age_grp) %>%
  dplyr::select(-age_grp_idx) %>%
  dplyr::select(-bin.10yr)

# Rescale to percentages for consistency before saving
risk_mat_with_mc[, c(2, 3, 5, 6)] <- 100 * risk_mat_with_mc[, c(2, 3, 5, 6)]

if (SAVE_OUTPUTS) {
  write.csv(risk_mat_with_mc, "synthetic_risk_data_with_misclassification.csv")
}
