#Mammography with sensitivity conditional on tumour diameter parameters W-F
beta1 <- 1.47 #Increases result in longer time to reach maximum sensitivity
beta2 <- 6.51 #Increases results in longer time to reach maximum sensitivity
sensitivity_max<-0.95 #This is a parameter in the model but could be hard coded here

#Create vector of cancer sizes (sensitivity pretty much always maxes out at 20mm)
Ca_size<-c(seq(from=0,to=20,by=1))

#Create vector of sensitivities by size
sens_size<-exp((Ca_size - beta2)/beta1)/(1+exp((Ca_size-beta2)/beta1))

#Set max sensitivity to 0.95
sens_size<-replace(sens_size,sens_size>sensitivity_max,sensitivity_max)

#Plot results
plot(Ca_size,sens_size,type="l",ylim=c(0,1),ylab=("Sensitivity"),xlab=("Cancer size (mm)"))

