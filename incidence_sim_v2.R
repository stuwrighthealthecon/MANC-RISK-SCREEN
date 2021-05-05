#Load necessary ONS data
Incidence_function <- function(){
  
#x<- cumsum(Incidence_Mortality[,2][21:101])
#Sample an incidence time (based on vector of probabilities of getting cancer at age t conditional on getting cancer and surviving to age t)
incidence_time_1 <- sample(x = Incidence_Mortality[,1][start_age:101],size = 1,prob = Incidence_Mortality[,2][start_age:101])
#add within year time (as uniform)
incidence_time <- incidence_time_1+ runif(1,0,1)
#back calculate the genesis time of the tumour

#First determine if screen detected or clinical detected in current data 

detect_mode <- 1 #clinical detected
if (incidence_time <= 70 & incidence_time >= 50 & runif(1,0,1)<prop_screen_detected){detect_mode <- 0} #screen detected

#size at detection - as number of tumour doublings in diameter from a 0.25mm diameter
  clin_detect_size_g <- rnorm(n = 1,mean = clin_detection_m,sd = clin_detection_sd)
  clin_detect_size_g[clin_detect_size_g < 4] <- 4 # to prevent unrealistic left tail
  clin_detect_size_g[clin_detect_size_g >= 9] <- 8.99 # to prevent unrealistic right tail
  clin_detect_size_g <- start_size*2^clin_detect_size_g
  ca_size_incidence <- clin_detect_size_g
  if(detect_mode == 0){
  screen_detect_size_g <- rnorm(n = 1,mean = screen_detection_m,sd = screen_detection_sd)
  screen_detect_size_g[screen_detect_size_g < 3.5] <- 3.5 # to prevent unrealistic left tail
  screen_detect_size_g[screen_detect_size_g >= 9] <- 8.99 # to prevent unrealistic right tail
  screen_doubles <- screen_detect_size_g #size as number of doublings
  screen_detect_size_g <- start_size*2^screen_detect_size_g 
  ca_size_incidence <- screen_detect_size_g
  #generate a valid potential clinical detection time - ie after screen detection
  if(clin_detect_size_g< screen_detect_size_g){
    
    clin_detect_size_g <- qnorm(runif(n = 1,min = pnorm(screen_doubles,clin_detection_m,clin_detection_sd),max = 1),mean = clin_detection_m,sd = clin_detection_sd)
    clin_detect_size_g[clin_detect_size_g >= 9] <- 8.999 # to prevent unrealistic right tail
    clin_detect_size_g <- start_size*2^clin_detect_size_g
  }
  }


#generate (non-cancer) mortality time conidtional on survivng to time t(peroid covered by BC survival function)
if (incidence_time<90){
mort_time <- sample(x = Incidence_Mortality[,1][incidence_time_1:101],size = 1,prob = Incidence_Mortality[,3][incidence_time_1:101])
}else{mort_time <- 101} # if get cancer at >90 then all time in simulation is covered 
# to replace mort_age
#How does this interact with the calculation of overall survival given breast cancer diagnosis?
#Breast Cancer Survial reverts to mort_age after 10 years

result<-c(incidence_time,detect_mode,ca_size_incidence, clin_detect_size_g,mort_time)
return(result)
}
