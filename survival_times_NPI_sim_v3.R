#Simulate survival by NPI Fong 2015 estimates (women aged 50-65)
#Needs to know NPI_Cat to generate a survival time from current age (currently 10-yr with cancer survival is irrespective of age), mort_age and age are need for those suviving beyond 10 years.

#Further assumption to guard against (reverse)lead-time bias is that cancer-specific survival is calculated from the age the cancer would have been clinically detected. Assumes no mortality effect of treatment.

Ca_survival_time <- function(NPI_cat, mort_age,age, CD_age){

gamma_NPI <- c(gamma_survival_1,gamma_survival_2,gamma_survival_3)
metastatic_survival <- c(meta_survival_49, meta_survival_69, meta_survival_99)

if (NPI_cat< 4){
  survival_time <- -(log(x = runif(1,0,1))/gamma_NPI[NPI_cat]) #inverse of cdf when rate is gamma_NPI[x]

#adjust for additional mortality at ages above 65
if (CD_age > 65){
  survival_time <- -(log(x = runif(1,0,1))/((Incidence_Mortality[min((floor(CD_age)+1),100),5]/Incidence_Mortality[66,5])*gamma_NPI[NPI_cat]))
}
#data are for 10-year survival, after 10 years assume that pop mortality rates apply
if(survival_time > 10){
  mort_age <- qweibull(p = runif(n = 1,min = pweibull(q = CD_age+10,shape = mortality_wb_a,scale = mortality_wb_b),max = 1),shape = mortality_wb_a, scale = mortality_wb_b)
  if(mort_age > time_horizon){mort_age <- time_horizon}
  survival_time <- mort_age - age
  }
}

if (NPI_cat == 5){
  if (age < 50){age_cat_M <- 1}
  if (age >=50 && age <70){age_cat_M <- 2}
  if (age >= 70){age_cat_M <- 3}
  survival_time <- (log(runif(1,0,1))/metastatic_survival[age_cat_M])
  #check lifetime does not exceed horizon and set to less than 100 if it does
  if (CD_age+survival_time >=100){survival_time <- time_horizon - CD_age}
}
if (NPI_cat == 4){
  survival_time <- (mort_age-CD_age) # no effect on mortality
}

if(CD_age+survival_time > time_horizon){survival_time <- time_horizon-CD_age}
result <- CD_age+survival_time
  return(result)
}