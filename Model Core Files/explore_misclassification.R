# In this script we show how to get beta distributed risk distributions for
# "true" risk based on comp[aring observed and estimated risk from Tyrer-Cuzick

EXPLORE_GPS <- FALSE # Set to true to look at a few alternative targets for fitting in the Gaussian process section

library(dplyr)
library(ggplot2)
library(magrittr)
library(pdftools)
library(purrr)
library(tidyr)

# Load in risk table from Brentnall 2018 supplement:
supp_data <- pdf_data("coi180011supp1_prod.pdf")[[5]][-c(1:15), ] %>% mutate(x = round(x/3),        #reduce resolution to minimise inconsistent coordinates
                                                                             y = round(y/3)) %>% 
  arrange(y, x) %>%                        #sort in reading order
  mutate(group = cumsum(!sapply(1:length(x), FUN=function(i){ # Gather nearby strings into same column (cutoff of 3 chosen by trial and error)
    if (i==1){
      return(TRUE)
    }else{
      return(
        (y[i-1]==y[i])&((x[i]-x[i-1])<=3))
    }
  }
  )
  )
  )%>%  #identify text with spaces and paste
  group_by(group) %>%
  summarise(x = first(x),
            y = first(y),
            text = paste(text, collapse = " ")) %>% 
  filter(!grepl('95%', text)) %>% # Remove column headings for 95% confidence intervals
  filter(!grepl('%)', text)) %>% # Remove percentage conversions for numbers
  group_by(y) %>% 
  mutate(colno = row_number()) %>%         #add column numbers for table data 
  ungroup() %>% 
  dplyr::select(text, colno, y) %>% 
  pivot_wider(names_from = colno, values_from = text) %>% #pivot into table format
  dplyr::select(-y) %>%
  set_colnames(c("Risk_class",
                 "N",
                 "FU",
                 "O",
                 "E",
                 "O_E",
                 "O_E 95% CI",
                 "IR",
                 "IRR",
                 "IRR 95% CI")) %>%
  dplyr::select(1:10) %>%
  slice(-c(1, 44, 45, 46))
model_idxs <- grep("Ty", supp_data$Risk_class)
supp_data$Model_age <- sapply(1:nrow(supp_data), FUN = function(i){supp_data[model_idxs[ceiling(i/7)],] %>% # Add column specifying age group and risk class
    as.character() %>%
    discard(is.na) %>%
    paste(collapse=" ")})
supp_data <- supp_data[-c(model_idxs), ]

supp_data[, c("N", "E", "O", "O_E")] <- supp_data[, c("N", "E", "O", "O_E")] %>% sapply(as.numeric)

# Calculate parameters for an empirical distribution based on numbers positive and negative:
alpha_vals <- supp_data$O[which(supp_data$Risk_class=="Total")] / supp_data$N[which(supp_data$Risk_class=="Total")]
beta_vals <- 1 - alpha_vals

# Now look at how well this fits...

# Get cdf of risk from data, compare to cdf from beta fit
supp_data_tc50 <- supp_data[which(supp_data$Model_age=="Tyrer-Cuzick (<50y)"), ]
data_quantiles <- cumsum(supp_data_tc50$N[2:6])/ supp_data_tc50$N[1]
heights <- supp_data_tc50$O[2:6] / supp_data_tc50$N[2:6]
heights[5] <- 1
beta_cdf <- pbeta(heights, alpha_vals[1], beta_vals[1])
cdf_df <- data.frame("data_quantiles"=data_quantiles,
                          "heights"=heights,
                          "beta_cdf"=beta_cdf)
h <- ggplot(cdf_df, aes(x=heights, y=data_quantiles)) +
  geom_line(colour="red") +
  geom_point(aes(x=heights, y=beta_cdf)) +
  scale_x_continuous(transform = "log10")

# Should see from plot that fit is fairly poor!
print(h)

# Try finding optimal parameters using least squares:
get_beta_rmse <- function(pars){
  alpha <- pars[1]
  beta <- pars[2]
  cdf <- pbeta(heights, alpha, beta)
  rmse <- (1/length(data_quantiles)) * (data_quantiles - cdf)^2 %>% sum() %>% sqrt()
  return(rmse)
}

par_estim <- optim(c(alpha_vals[1], beta_vals[1]), get_beta_rmse)

beta_cdf <- pbeta(heights, par_estim$par[1], par_estim$par[2])
cdf_df <- data.frame("data_quantiles"=data_quantiles,
                          "heights"=heights,
                          "beta_cdf"=beta_cdf)
h <- ggplot(cdf_df, aes(x=heights, y=data_quantiles)) +
  geom_line(colour="red") +
  geom_point(aes(x=heights, y=beta_cdf)) +
  scale_x_continuous(transform = "log10")

# Fit now looks much better!
print(h)

# This fit does not preserve true mean:
estim_mean <- par_estim$par[1] / (par_estim$par[1] + par_estim$par[2])
true_mean <- supp_data_tc50$O[1] / supp_data_tc50$N[1]

# Alternative approach: force fits to match mean
get_beta_rmse_with_mean <- function(alpha){
  beta <- alpha / true_mean - alpha
  cdf <- pbeta(heights, alpha, beta)
  rmse <- (1/length(data_quantiles)) * (data_quantiles - cdf)^2 %>% sum() %>% sqrt()
  return(rmse)
}
par_estim <- optim(alpha_vals[1], get_beta_rmse_with_mean)
alpha_estim <- par_estim$par[1]
beta_estim <- alpha_estim / true_mean - alpha_estim
estim_mean <- alpha_estim / (alpha_estim + beta_estim)

beta_cdf <- pbeta(heights, alpha_estim, beta_estim)
cdf_df <- data.frame("data_quantiles"=data_quantiles,
                          "heights"=heights,
                          "beta_cdf"=beta_cdf)
# Should see much worse eyeball fit:
h <- ggplot(cdf_df, aes(x=heights, y=data_quantiles)) +
  geom_line(colour="red") +
  geom_point(aes(x=heights, y=beta_cdf)) +
  scale_x_continuous(transform = "log10")
print(h)

#### Do least-squares beta regression systematically for each age/model combination

n_classes <- supp_data$Model_age %>% unique() %>% length()
beta_fits_df <- data.frame("obs_mean"=numeric(length=n_classes),
                           "alpha_e"=numeric(length=n_classes),
                           "beta_e"=numeric(length=n_classes),
                           "mean_e"=numeric(length=n_classes),
                           "err_e"=numeric(length=n_classes),
                           "alpha_o"=numeric(length=n_classes),
                           "beta_o"=numeric(length=n_classes),
                           "mean_o"=numeric(length=n_classes),
                           "err_o"=numeric(length=n_classes),
                           row.names = unique(supp_data$Model_age))

for (m in unique(supp_data$Model_age)){
  supp_data_m <- supp_data[which(supp_data$Model_age==m), ]
  data_quantiles <- cumsum(supp_data_m$N[2:6])/ supp_data_m$N[1]
  
  # First do expected
  heights <- supp_data_m$E[2:6] / supp_data_m$N[2:6]
  heights[5] <- 1
  get_beta_rmse <- function(pars){
    alpha <- pars[1]
    beta <- pars[2]
    cdf <- pbeta(heights, alpha, beta)
    rmse <- (1/length(data_quantiles)) * (data_quantiles - cdf)^2 %>% sum() %>% sqrt()
    return(rmse)
  }
  alpha_0 <- supp_data_m$E[which(supp_data_m$Risk_class=="Total")]
  beta_0 <- supp_data_m$N[which(supp_data_m$Risk_class=="Total")] - supp_data_m$E[which(supp_data_m$Risk_class=="Total")]
  par_estim_e <- optim(c(alpha_0, beta_0), get_beta_rmse)
  
  beta_cdf <- pbeta(heights, par_estim_e$par[1], par_estim_e$par[2])
  cdf_df <- data.frame("data_quantiles"=data_quantiles,
                       "heights"=heights,
                       "beta_cdf"=beta_cdf)
  
  h <- ggplot(cdf_df, aes(x=heights, y=data_quantiles)) +
    geom_line(colour="red") +
    geom_point(aes(x=heights, y=beta_cdf)) +
    scale_x_continuous(transform = "log10") +
    ggtitle(paste(m, ", expected", sep=""))
  print(h)
  
  # Now do estimates
  heights <- supp_data_m$O[2:6] / supp_data_m$N[2:6]
  heights[5] <- 1
  get_beta_rmse <- function(pars){
    alpha <- pars[1]
    beta <- pars[2]
    cdf <- pbeta(heights, alpha, beta)
    rmse <- (1/length(data_quantiles)) * (data_quantiles - cdf)^2 %>% sum() %>% sqrt()
    return(rmse)
  }
  alpha_0 <- supp_data_m$O[which(supp_data_m$Risk_class=="Total")]
  beta_0 <- supp_data_m$N[which(supp_data_m$Risk_class=="Total")] - supp_data_m$O[which(supp_data_m$Risk_class=="Total")]
  par_estim_o <- optim(c(alpha_0, beta_0), get_beta_rmse)
  
  beta_cdf <- pbeta(heights, par_estim_o$par[1], par_estim_o$par[2])
  cdf_df <- data.frame("data_quantiles"=data_quantiles,
                       "heights"=heights,
                       "beta_cdf"=beta_cdf)
  
  h <- ggplot(cdf_df, aes(x=heights, y=data_quantiles)) +
    geom_line(colour="red") +
    geom_point(aes(x=heights, y=beta_cdf)) +
    scale_x_continuous(transform = "log10") +
    ggtitle(paste(m, ", observed", sep=""))
  print(h)
  
  beta_fits_df[m, ] <- c(alpha_0,
                         par_estim_e$par[1],
                         par_estim_e$par[2],
                         par_estim_e$par[1] / (par_estim_e$par[1] + par_estim_e$par[2]),
                         par_estim_e$value,
                         par_estim_o$par[1],
                         par_estim_o$par[2],
                         par_estim_o$par[1] / (par_estim_e$par[1] + par_estim_e$par[2]),
                         par_estim_o$value)
}

save(beta_fits_df, file="fitted_risk_distributions.csv")

#### Next problem: inferring life risk from ten year risk

risk_mat<-read.csv("synthetic_risk_data.csv") %>% dplyr::select(-"X")

# Inspect visually:

p <- ggplot(risk_mat, aes(x=syn.X10yr, y=syn.life)) +
  geom_point() +
  xlab("Ten year estimated risk") +
  ylab("Lifetime estimated risk")
p

# Try splitting up by age:
risk_mat <- risk_mat %>% mutate(age_grp=case_when(syn.Age<50 ~ "<50",
                                                  syn.Age>=50 & syn.Age<60 ~ "50-60",
                                                  syn.Age>=60 ~ "60+"))
p <- ggplot(risk_mat, aes(x=syn.X10yr, y=syn.life, colour=age_grp)) +
  geom_point(alpha=.2) +
  xlab("Ten year estimated risk") +
  ylab("Lifetime estimated risk") +
  scale_color_discrete()
p

# Should see obvious difference in relationship, so need to include age in regression model
library(GauPro)

risk_mat$bin.10yr <- risk_mat$syn.X10yr %>% cut(breaks=seq(from=0, to=100, by=.1))

smpl <- risk_mat[which(risk_mat$age_grp=="50-60"),] %>%
  slice_sample(by=bin.10yr, n=1)

# Simple version:
x <- smpl$syn.X10yr
y <- smpl$syn.life
kern <- Matern52$new(0)
gp <- GauPro(x, y, kernel=kern)
ypred <- gp$predict(risk_mat$syn.X10yr, se=T)
y_upper <- ypred$mean + 2*ypred$se
y_lower <- ypred$mean - 2*ypred$se
p <- ggplot(risk_mat, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(aes(x=risk_mat$syn.X10yr, y=ypred$mean), colour="red") +
  geom_line(aes(x=risk_mat$syn.X10yr, y=y_upper), colour="blue") +
  geom_line(aes(x=risk_mat$syn.X10yr, y=y_lower), colour="blue")
print(p)

# Last point is a bit of an outlier, try retraining without:
smpl_order <- smpl$syn.X10yr %>% order() %>% head(-1) # Order samples by x value and drop last one
x_short <- x[smpl_order]
y_short <- y[smpl_order]
kern <- Matern52$new(0)
gp_short <- GauPro(x_short, y_short, kernel=kern)
ypred <- gp_short$predict(risk_mat$syn.X10yr, se=T)
y_upper <- ypred$mean + 2*ypred$se
y_lower <- ypred$mean - 2*ypred$se
p <- ggplot(risk_mat, aes(syn.X10yr, syn.life)) +
  geom_point() +
  geom_line(aes(x=risk_mat$syn.X10yr, y=ypred$mean), colour="red") +
  geom_line(aes(x=risk_mat$syn.X10yr, y=y_upper), colour="blue") +
  geom_line(aes(x=risk_mat$syn.X10yr, y=y_lower), colour="blue")
print(p)

#### In this optional section we look at transformations of the data to improve the GP fit.
# Note that the code below this section will still use gp_short, the Gaussian process
# fit to the raw risk data with the final datapoint removed, to generate synthetic
# lifetime risks.

if (EXPLORE_GPS){

# Version that estimates log(y-x)

x <- smpl$syn.X10yr
y <- (smpl$syn.life - smpl$syn.X10yr) %>% log()
gp <- GauPro(x, y, kernel=Exponential$new(0))

# Try plotting risk prediction:

ypred <- function(z){z + (gp$predict(z) %>% exp())}
yse <- function(z){exp(gp$predict(z) + 0.5*gp$predict(z, se=T)$se^2) * sqrt(exp(gp$predict(z, se=T)$se^2) - 1)}
plot(x, smpl$syn.life)
curve(ypred(x), add=T, col=2)
curve(ypred(x)+2*yse(x), add=T, col=4)
curve(ypred(x)-2*yse(x), add=T, col=4)

# Alternative: estimate log(dy/dx)
smpl_order <- smpl$syn.X10yr %>% order()
x <- smpl$syn.X10yr[smpl_order][2:length(smpl_order)]
y <- diff(smpl$syn.life[smpl_order]) / diff(smpl$syn.X10yr[smpl_order])

x <- x[which(y>0)]
y <- y[which(y>0)]
gp <- GauPro(x, log(y), kernel=Exponential$new(0))
plot(x, y)
yse <- function(z){exp(gp$predict(z) + 0.5*gp$predict(z, se=T)$se^2) * sqrt(exp(gp$predict(z, se=T)$se^2) - 1)}
curve(exp(gp$predict(x)), add=T, col=2)
curve(exp(gp$predict(x))+yse(x), add=T, col=4)
curve(exp(gp$predict(x)+ 1.1*min(y))-yse(x), add=T, col=4)
}

#### Now add "true" risks to risk_mat, assuming we're working with prediction
# with density.

bf_dens_only <- beta_fits_df[grep("density", rownames(beta_fits_df)),]
supp_data_dens_only <- supp_data %>%
                       filter(grepl("density", Model_age)) %>%
                       filter(grepl("Total", Risk_class))
cutoffs <- c(0, 50, 60)

# Function for synthetic to true risk mapping
# Problem: this can generate negatives
true_risk_from_syn <- function(syn_risk, age){
  age_class <- which(cutoffs<age)  %>% max()
  beta_pars <- c(bf_dens_only$alpha_o[age_class],
                 bf_dens_only$beta_o[age_class])
  true_risk <- 100*(syn_risk/100 * as.numeric(supp_data_dens_only$O_E[age_class]) -
                bf_dens_only$obs_mean[age_class] + rbeta(1, beta_pars[1], beta_pars[2]) )
  return(true_risk)
}

true_risk_from_syn_v <- Vectorize(true_risk_from_syn)

# Alternative approach: fit a maximum likelihood beta distribution to synthetic
# risk and map between quantiles

# TODO: STRATIFY BY AGE
beta_logpmf <- function(alpha, risk_data){
  if((alpha[[1]]>=0)&(alpha[[2]]>=0)){
    lpmf<- dbeta(.01*risk_data$syn.X10yr, alpha[[1]], alpha[[2]]) %>% log() %>% sum()
    return(-lpmf)
  }else{return(100000)}
}

syn_beta_pars <- sapply(c("<50", "50-60"), FUN=function(str){
  pars <- optim(c(1., 100.), beta_logpmf, risk_data=risk_mat[grep(str, risk_mat$age_grp),])$par
  return(pars)
  }
  )
true_risk_from_syn <- function(syn_risk, age){
  # print(syn_risk)
  age_class <- which(cutoffs<age)  %>% max()
  
  beta_pars_e <- c(bf_dens_only$alpha_e[age_class],
                   bf_dens_only$beta_e[age_class])
  # Get cdf evaluated at synthetic risk
  cdf <- pbeta(.01*syn_risk,
                    beta_pars_e[1],
                    beta_pars_e[2])
  # print(cdf)
  beta_pars_o <- c(bf_dens_only$alpha_o[age_class],
                 bf_dens_only$beta_o[age_class])
  # print(c(syn_beta_pars[1, age_class],
  #         syn_beta_pars[2, age_class]))
  # print(beta_pars)
  # Get quantile associated with same position along true risk distribution:
  true_risk <- 100*qbeta(cdf, beta_pars_o[1], beta_pars_o[2])
  # print(pbeta(.01*true_risk, beta_pars[1], beta_pars[2]))
  return(true_risk)
}

true_risk_from_syn_v <- Vectorize(true_risk_from_syn)

# Create copy of risk_mat with misclassification
risk_mat_with_mc <- risk_mat %>%
                    mutate(syn.X10yr.true = true_risk_from_syn_v(syn.X10yr,
                                                                 syn.Age)) %>%
                    mutate(syn.life.true.dist = gp_short$predict(syn.life,
                                                                 se=T)) %>%
                    mutate(syn.life.true = rnorm(nrow(risk_mat), mean=syn.life.true.dist$mean,
                                                 sd = syn.life.true.dist$se)) %>%
                    select(-syn.life.true.dist)

# Plot to check:

ggplot(risk_mat_with_mc, aes(syn.X10yr, syn.X10yr.true, col=age_grp)) + geom_point() + geom_abline()
