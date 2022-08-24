install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("mgcv")
install.packages("BCEA")
install.packages("devtools")
install.packages("here")


# Libraries -----------------------------------------------------------------------------------
library("mgcv")
library("BCEA")
library(tidyverse)
library(here)
library(magrittr)
library(devtools)

# Load data -----------------------------------------------------------------------------------
alt_names<-c("noscreening","procas","tertiles","3yr","2yr","5yr",
             "10yr","lowrisk5yr","lowrisk6yr","fullstrat")
subD <- "GitHub/MANC-RISK-SCREEN/Model Core Files/GAM inputs"
load(here(subD, "PSA_values.RData"))

PSA_all <- tibble(strat = factor(alt_names)) %>%
  pmap_dfr(~ read.csv(here(subD, paste0(.x, "_PSA.csv"))) %>%
             set_names(c("id","QALY","cost","screens","cancer","screen_detected","LY")) %>%
             bind_cols(alternative = .x,
                       PSA_all_p))

# CE results ----------------------------------------------------------------------------------
## CE plane with raw results from monte carlo --------------------------------------------------
incCU_PSA <- PSA_all %>%
  filter(alternative != "noscreening") %>%
  dplyr::select(id, alternative, cost, QALY) %>%
  left_join(PSA_all %>%
              filter(alternative == "noscreening") %>%
              dplyr::select(id,
                            baseCost = cost,
                            baseQALY = QALY),
            by = "id") %>%
  mutate(delta_cost = cost - baseCost,
         delta_QALY = QALY - baseQALY)

incCU_PSA %>%
  ggplot() +
  geom_point(aes(x=delta_QALY, y = delta_cost, colour = alternative))

## Incremental cost effectiveness at mean ------------------------------------------------------
incCU <- PSA_all %>%
  group_by(StratName = alternative) %>%
  summarise(Costs = mean(cost),
            QALYs = mean(QALY)) %>%
  mutate(ID = row_number())

source_gist("https://gist.github.com/gbrlrgrs/5b905bdb5ce597c01f3551746154b45a", encoding = "UTF-8")
incCU %>%
  fnIncCU()

#  Fit GAMs ------------------------------------------------------------------------------------
## QALYs ---------------------------------------------------------------------------------------
modQ <- bam(data = PSA_all,
            formula = QALY ~
              s(PSA_gamma_survival_1, by = alternative, bs = "tp") +
              s(PSA_gamma_survival_2, by = alternative, bs = "tp") +
              s(PSA_gamma_survival_3, by = alternative, bs = "tp") +
              s(PSA_meta_survival_54, by = alternative, bs = "tp") +
              s(PSA_meta_survival_74, by = alternative, bs = "tp") +
              s(PSA_meta_survival_99, by = alternative, bs = "tp") +
              s(PSA_beta_1, by = alternative, bs = "tp") +
              s(PSA_beta_2, by = alternative, bs = "tp") +
              s(PSA_VDG1_sen, by = alternative, bs = "tp") +
              s(PSA_VDG2_sen, by = alternative, bs = "tp") +
              s(PSA_VDG3_sen, by = alternative, bs = "tp") +
              s(PSA_VDG4_sen, by = alternative, bs = "tp") +
              s(PSA_MRI_cdr, by = alternative, bs = "tp") +
              s(PSA_US_cdr, by = alternative, bs = "tp") +
              s(PSA_log_norm_mean, by = alternative, bs = "tp") +
              s(PSA_log_norm_sd, by = alternative, bs = "tp") +
              s(PSA_util_1to3, by = alternative, bs = "tp") +
              s(PSA_util_4, by = alternative, bs = "tp") +
              alternative)
summary(modQ)

## Costs ---------------------------------------------------------------------------------------
modC <- bam(data = PSA_all,
            formula = cost ~ 
              s(PSA_cost_strat, by = alternative, bs = "tp") +
              s(PSA_costvar, by = alternative, bs = "tp") +
              s(PSA_gamma_survival_1, by = alternative, bs = "tp") +
              s(PSA_gamma_survival_2, by = alternative, bs = "tp") +
              s(PSA_gamma_survival_3, by = alternative, bs = "tp") +
              s(PSA_meta_survival_54, by = alternative, bs = "tp") +
              s(PSA_meta_survival_74, by = alternative, bs = "tp") +
              s(PSA_meta_survival_99, by = alternative, bs = "tp") +
              s(PSA_beta_1, by = alternative, bs = "tp") +
              s(PSA_beta_2, by = alternative, bs = "tp") +
              s(PSA_VDG1_sen, by = alternative, bs = "tp") +
              s(PSA_VDG2_sen, by = alternative, bs = "tp") +
              s(PSA_VDG3_sen, by = alternative, bs = "tp") +
              s(PSA_VDG4_sen, by = alternative, bs = "tp") +
              s(PSA_MRI_cdr, by = alternative, bs = "tp") +
              s(PSA_US_cdr, by = alternative, bs = "tp") +
              s(PSA_log_norm_mean, by = alternative, bs = "tp") +
              s(PSA_log_norm_sd, by = alternative, bs = "tp") +
              alternative)
summary(modC)

# Model fit ---------------------------------------------------------------------------
## Expected costs and QALYs fitted at the mean of sampled parameter values --------------------
ftd <- PSA_all_p %>%
  summarise(across(everything(), mean)) %>%
  crossing(tibble(alternative = alt_names))

ftdCE <- ftd %>%
  bind_cols(ftdQ = modQ %>%
              predict(newdata = ftd),
            ftdC = modC %>%
              predict(newdata = ftd))
incCU %>%
  left_join(ftdCE, by = c("StratName" = "alternative")) %>%
  select(ID, StratName, Costs, QALYs, ftdC, ftdQ)
