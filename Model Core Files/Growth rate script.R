#Set tumour growth rate parameters
log_norm_mean <- 1.07
log_norm_sd <- 1.31
max_size <- 128 #mm diameter
start_size <- 0.25 #starting size of tumours, diameter in mm
Vc = (4/3)*pi*(start_size/2)^3 #Volume at start
Vm = (4/3)*pi*(max_size/2)^3 #Max volume

#Get estimates for Mean growth rate and upper and lower confidence intervals
grow_rate<-rlnorm(runif(100000,0,1),meanlog=log_norm_mean,sdlog=sqrt(log_norm_sd))
rate_centiles<-quantile(grow_rate,c(0.025,0.5,0.975))

#Estimate time taken to grow to different sizes
upperbound<-seq(from=10,to=100,b=5)
lowerbound<-seq(from=5,to=95,by=5)

#Function estimates time take to grow from one size to another
#Creates a vector where each element is a cumulative sum of the time taken to grow 
#to that size
growth_time<-function(upperbound,lowerbound,grow_rate){
  t<-vector(length=length(upperbound))
  for (i in 1:length(upperbound)) {
   t[i]<-((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(upperbound[i]/2)^3))^0.25-1))/(0.25*grow_rate)) - 
  ((log((Vm/Vc)^0.25-1)-log((Vm/((4/3)*pi*(lowerbound[i]/2)^3))^0.25-1))/(0.25*grow_rate))}
  t<-cumsum(t)
return(t)
}

#Create a data.frame of time taken to grow to different sizes for lower, mean, 
#and upper estimates of the growth rate
grow_time_all<-data.frame(matrix(nrow=20,ncol=4))
colnames(grow_time_all)<-c("Size","Lower","Mean","Upper")
grow_time_all[1,]<-c(0,0,0,0)
grow_time_all$Size[2:20]<-upperbound

grow_time_all$Lower[2:20]<-growth_time(upperbound,lowerbound,grow_rate=rate_centiles[1])
grow_time_all$Mean[2:20]<-growth_time(upperbound,lowerbound,grow_rate=rate_centiles[2])
grow_time_all$Upper[2:20]<-growth_time(upperbound,lowerbound,grow_rate=rate_centiles[3])

#Plot results
plot(grow_time_all$Mean,grow_time_all$Size,
     type="l",xlab="Time (Years)",ylab="Size (mm)",xlim=c(0,20))
lines(grow_time_all$Lower,grow_time_all$Size,lty=2)
lines(grow_time_all$Upper,grow_time_all$Size,lty=2)
