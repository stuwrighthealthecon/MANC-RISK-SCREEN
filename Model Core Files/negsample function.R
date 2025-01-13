screen_strategy<-3


if(SEPARATE_SAMPLES){
  if(MISCLASS){
    load("Risksamplewithmisclass/negsample.Rdata")}else{
   load("Risksample/negsample.Rdata")
    }
  }

negsample<-data.frame("risk_group"=negsample$risk_group,
                      "MRI_screen"=negsample$MRI_screen,
                      "US_screen"=negsample$US_screen,
                      "risk_predicted"=negsample$risk_predicted,
                      "feedback"=negsample$feedback,
                      "interval_change"=negsample$interval_change,
                      "life_expectancy"=negsample$life_expectancy)

if(screen_strategy==3){
  screen_times <- low_risk_screentimes
}
if(screen_strategy==4){
  screen_times <- med_risk_screentimes
}
if(screen_strategy==5){
  screen_times <- seq(screen_startage, screen_startage+(5*4),5)
}
if(screen_strategy==6){
  screen_times <- seq(screen_startage, screen_startage+10,10)
}

for(i in 1:length(screen_times)){
  negsample[,7+i]<-numeric(length(negsample$negsample.risk_group))
}

negsample[,8]<-rbinom(length(negsample$risk_group),1,uptakefirstscreen)

for (i in 1:length(screen_times)-1){
negsample[,8+i]<-ifelse(rowSums(negsample[8:(7+i)])>=1,
                      rbinom(length(negsample$risk_group),1,uptakeotherscreen),
                      rbinom(length(negsample$risk_group),1,uptakenoscreen))
}
negsample$total_screens<-rowSums(negsample[8:length(negsample[1,])])

#################PEOPLE CAN HAVE SCREENING AFTER THEY DIE CURRENTLY##############

for (i in 1:length(screen_times)){
  negsample[,7+i]<-negsample[,7+i]*(cost_screen*((1/((1+discount_cost)^(screen_times[i]-screen_startage)))))
}
negsample$screencost<-rowSums(negsample[8:length(negsample[1,])])


results<-data.frame(rep(0,length=length(negsample$risk_group)),
                    negsample$screencost,
                    negsample$total_screens,
                    rep(0,length=length(negsample$risk_group)),
                    rep(0,length=length(negsample$risk_group)),
                    rep(0,length=length(negsample$risk_group)),
                    rep(screen_strategy,length=length(negsample$risk_group)),
                    rep(0,length=length(negsample$risk_group)),
                    negsample$life_expectancy-rep(start_age,length=length(negsample$risk_group)),
                    rep(0,length=length(negsample$risk_group)),
                    rep(0,length=length(negsample$risk_group)),
                    negsample$life_expectancy,
                    rep(0,length=length(negsample$risk_group)))
names(results) <- c('QALY',
                    'Cost',
                    'Screens',
                    "Cancer Diagnosed Age",
                    "Cancer",
                    "screen detected",
                    "alternative",
                    "Growth rate",
                    "Life Years",
                    "Stage",
                    "Cancer Size",
                    "Death Age",
                    "Cancer Screen Number")