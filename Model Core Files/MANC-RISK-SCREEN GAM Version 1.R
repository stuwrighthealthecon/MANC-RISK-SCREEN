install.packages("dplyer")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("mgcv")
install.packages("BCEA")

#Libraries
library("dplyr")
library("tidyr")
library("ggplot2")
library("mgcv")
library("BCEA")

load("PSA_values.RData")

#Factor variable to indicate which alternative the results belong to
alt_num <- 10 #number of alternatives
alternative <- c(rep(1,1000))
for (i in 2:alt_num){
  alternative <- c(alternative,c(rep(i,1000)))
}
alt_names<-c("noscreening","procas","tertiles","3yr","2yr","5yr",
                   "10yr","lowrisk5yr","lowrisk6yr","fullstrat")
alternative <- as.data.frame(factor(alternative, labels = alt_names))
colnames(alternative) <- c("alternative")

#Load the results data and store in a list
#noscreening
PSA_results <- read.csv(file = "noscreening_PSA.csv") 
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- PSA_results

#Procas
PSA_results <- read.csv(file = "procas_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#tertiles
PSA_results <- read.csv(file = "tertiles_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#3yr
PSA_results <- read.csv(file = "3yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#2yr
PSA_results <- read.csv(file = "2yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#5yr
PSA_results <- read.csv(file = "5yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results)

#10yr
PSA_results <- read.csv(file = "10yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#lowrisk5yr
PSA_results <- read.csv(file = "lowrisk5yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#lowrisk6yr
PSA_results <- read.csv(file = "lowrisk6yr_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#fullstrat
PSA_results <- read.csv(file = "fullstrat_PSA.csv")
colnames(PSA_results) <- c("id","QALY","cost","screens","cancer","screen_detected","LY")
PSA_results <- bind_cols(PSA_all_p,PSA_results)
PSA_results_all <- bind_rows(PSA_results_all, PSA_results) 

#Data frames with results and parameters
PSA_all <- bind_cols(PSA_results_all, alternative)

#CE plane with raw results from monte carlo

#Make wide, QALY and cost by alternative
base_data <- (PSA_all) %>%
  dplyr::select(one_of(c("id","alternative","QALY","cost"))) %>%
  pivot_wider(names_from = alternative, values_from = c(QALY, cost))
#Get incremental QALY and cost for each alternative with no screening as the reference
base_data <- base_data %>%
  mutate(delta_QALY_noscreening = QALY_noscreening - QALY_noscreening) %>%
  mutate(delta_QALY_procas = QALY_procas - QALY_noscreening) %>%
  mutate(delta_QALY_tertiles = QALY_tertiles - QALY_noscreening) %>%
  mutate(delta_QALY_3yr = QALY_3yr - QALY_noscreening) %>%
  mutate(delta_QALY_2yr = QALY_2yr - QALY_noscreening) %>%
  mutate(delta_QALY_5yr = QALY_5yr - QALY_noscreening) %>%
  mutate(delta_QALY_10yr = QALY_10yr - QALY_noscreening) %>%
  mutate(delta_QALY_lowrisk5yr = QALY_lowrisk5yr - QALY_noscreening) %>%
  mutate(delta_QALY_lowrisk6yr = QALY_lowrisk6yr - QALY_noscreening) %>%
  mutate(delta_QALY_fullstrat = QALY_fullstrat - QALY_noscreening) %>%
  
  mutate(delta_cost_noscreening = cost_noscreening - cost_noscreening) %>%
  mutate(delta_cost_procas = cost_procas - cost_noscreening) %>%
  mutate(delta_cost_tertiles = cost_tertiles - cost_noscreening) %>%
  mutate(delta_cost_3yr = cost_3yr - cost_noscreening) %>%
  mutate(delta_cost_2yr = cost_2yr - cost_noscreening) %>%
  mutate(delta_cost_5yr = cost_5yr - cost_noscreening) %>%
  mutate(delta_cost_10yr = cost_10yr - cost_noscreening) %>%
  mutate(delta_cost_lowrisk5yr = cost_lowrisk5yr - cost_noscreening) %>%
  mutate(delta_cost_lowrisk6yr = cost_lowrisk6yr - cost_noscreening) %>%
  mutate(delta_cost_fullstrat = cost_fullstrat - cost_noscreening)

#Fitting the GAM
#colnames(PSA_all_p) <- var_names
#Matrix for fitted values
#Fitted_mat <-matrix(0,nrow=1000,ncol=alt_num*2)
#Data frame to store fitted values, parameter values, QALY and cost
fitted <- dplyr::select(PSA_all, one_of(c("id","alternative","QALY","cost"))) %>%
  mutate(RID = row_number())

#QALY model
rm(PSA_results_all)
rm(PSA_all_p)

model1 <- gam(data = PSA_all, formula = QALY ~ s(PSA_gamma_survival_1, by = alternative, bs = "tp")+ s(PSA_gamma_survival_2, by = alternative, bs = "tp")+s(PSA_gamma_survival_3, by = alternative, bs = "tp")+s(PSA_meta_survival_54, by = alternative, bs = "tp")+ s(PSA_meta_survival_74, by = alternative, bs = "tp") + s(PSA_meta_survival_99, by = alternative, bs = "tp") + s(PSA_beta_1, by = alternative, bs = "tp") + s(PSA_beta_2, by = alternative, bs = "tp") + s(PSA_VDG1_sen, by = alternative, bs = "tp") + s(PSA_VDG2_sen, by = alternative, bs = "tp") + s(PSA_VDG3_sen, by = alternative, bs = "tp") + s(PSA_VDG4_sen, by = alternative, bs = "tp") +s(PSA_MRI_cdr, by = alternative, bs = "tp")+ s(PSA_US_cdr, by = alternative, bs = "tp")+ s(PSA_log_norm_mean, by = alternative, bs = "tp") + s(PSA_log_norm_sd, by = alternative, bs = "tp") + s(PSA_util_1to3, by = alternative, bs = "tp") + s(PSA_util_4, by = alternative, bs = "tp") + alternative)
summary(model1)

#Cost model
model2 <- gam(data = PSA_all, formula = cost ~ s(PSA_cost_strat, by = alternative, bs = "tp") + s(PSA_costvar, by = alternative, bs = "tp") + alternative)
                #s(PSA_gamma_survival_1, by = alternative, bs = "tp") + s(PSA_gamma_survival_2, by = alternative, bs = "tp") + s(PSA_gamma_survival_3, by = alternative, bs = "tp")+ s(PSA_meta_survival_54, by = alternative, bs = "tp") + s(PSA_meta_survival_74, by = alternative, bs = "tp") + s(PSA_meta_survival_99, by = alternative, bs = "tp") + s(PSA_beta_1, by = alternative, bs = "tp") + s(PSA_beta_2, by = alternative, bs = "tp") + s(PSA_VDG1_sen, by = alternative, bs = "tp") + s(PSA_VDG2_sen, by = alternative, bs = "tp") + s(PSA_VDG3_sen, by = alternative, bs = "tp") + s(PSA_VDG4_sen, by = alternative, bs = "tp") +s(PSA_MRI_cdr, by = alternative, bs = "tp")+ s(PSA_US_cdr, by = alternative, bs = "tp")+ s(PSA_log_norm_mean, by = alternative, bs = "tp") + s(PSA_log_norm_sd, by = alternative, bs = "tp") + s(PSA_cost_strat, by = alternative, bs = "tp") + s(PSA_costvar, by = alternative, bs = "tp") + alternative)
summary(model2)
