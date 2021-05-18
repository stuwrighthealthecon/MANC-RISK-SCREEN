install.packages("doParallel")
library("doParallel")

#Set working directory
setwd(dir="")

#Register number of cores for foreach loop
registerDoParallel(cores=5)

#Set timer to record duration of simulation
ptm <- proc.time()

#Load files containing functions
source()
source()
source()
source()

#Set loop numbers
inum<-1000000
jnum<-1

#Set screening programme related parameters
#Turn supplemental Screening (MRI and US) on (1) or off (0)
supplemental_screening<-0

#Set the screening strategy: 1=PROCAS, 2=Optimal, 3=Risk tertiles, 4=3 yearly, 5=2 yearly, 6=5 yearly, 7=2 rounds at 50 and 60, Other num=no screening
screen_strategy<-

#Age of an individual at start of simulation
start_age<-38

#Set time horizon
time_horizon<-100

#Set health and cost discount rates
discount_health<-0.035
discount_cost<-0.035

######################################End Programme parameters###################################################

#############################################################################################################

#Set fixed parameters (non PSA varying)
#Set proportion of cancers in screening age range detected by screen
screen_detected<-0.5

#Set the mean and standard deviation of the doubling rate for tumours
screen_detection_m<-4.12
screen_detection_sd<-3.93

#Set parameters of a Weibull survival curve to represent all cause mortality
acmmortality_wb_a<-8.97
acmmortality_wb_b<-86.74

#Set parameters for all cause mortality following breast cancer
gamma_survival_3<-exp(-2.465) #exponential distribution scale parameter NPI 3
gamma_survival_2<-exp(-4.023) #exponential distribution scale parameter NPI 2
gamma_sruvival_1<-exp(-5.413) #exponential distribution scale parameter NPI 1



#Set loop to divide i loop into 10 sub-loops in case of simulation break
for (ii in 10) {
  
  
  
}

