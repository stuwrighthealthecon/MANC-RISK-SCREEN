#Screen test results simulation (function)
#Updated 0915 with weedon-fekjaer 2008 estimates
#args (inputs) are tumour diameter, VDG, MRI_screening(0/1), US_screening(0/1) 
screening_result <- function(Ca_size,VDG,MRI_screening,US_screening){

  
#Caculate size/density specific sensitivity  
Sensitivity <- if(
    exp((Ca_size - beta2)/beta1)/(1+exp((Ca_size-beta2)/beta1))>0.95){0.95}else{exp((Ca_size - beta2)/beta1)/(1+exp((Ca_size-beta2)/beta1))} #use to set max sensitivity 0.95

dense_OR <- (Sen_VDG[VDG]/(1-Sen_VDG[VDG]))/(Sen_VDG_av/(1-Sen_VDG_av))
Sensitivity_dense <- ((Sensitivity/(1-Sensitivity))*dense_OR)/(1+((Sensitivity/(1-Sensitivity))*dense_OR))

Sensitivity <- Sensitivity_dense

rnd_1 <- runif(1,0,1) # random number used to compare to Sensitivity with and without supplemental screening
  if(rnd_1<Sensitivity){
    Screen_detected_ca <- 1
    Mammo_detected_ca <- 1 #keep track of which stage it is detected
    }else{Screen_detected_ca <-0
          Mammo_detected_ca <- 0} #test if cancer is detected by screening mammo
  #Additional screening depends on density/risk
  if(Screen_detected_ca == 0){
    if(MRI_screening == 1){
      MRI_supp_odds <- (Sensitivity/(1-Sensitivity))*((MRI_cdr+Mammo_cdr)/Mammo_cdr)
      MRI_supp_sens <- MRI_supp_odds/(MRI_supp_odds+1)
      if(rnd_1 < MRI_supp_sens){
        Screen_detected_ca <- 1
        MRI_detected_ca <- 1
      }else{MRI_detected_ca <- 0}}else{MRI_detected_ca <- 0}
  
    if(US_screening==1){
    US_supp_odds <- (Sensitivity/(1-Sensitivity))*((US_cdr+Mammo_cdr)/Mammo_cdr)
    US_supp_sens <- US_supp_odds/(US_supp_odds+1)
    if(rnd_1 < US_supp_sens){
      Screen_detected_ca <-1
      US_detected_ca <- 1
    }else{US_detected_ca <- 0}}else{US_detected_ca <- 0}}else{US_detected_ca <- 0
                                                              MRI_detected_ca <- 0}
    #Uses estimate based on increased cancer detection rate - check this then matches the number of interval cancers that are now found at screening

  
  result <- c(Screen_detected_ca,Mammo_detected_ca,MRI_detected_ca,US_detected_ca)
    return(result)
}