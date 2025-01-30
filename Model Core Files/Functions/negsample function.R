##############################Function for estaimating outcomes for non-cancer############

negsamplefn<-function(screen_strategy,MISCLASS){

#Load appropriate data
  if(MISCLASS){
    load("Risksamplewithmisclass/negsample.Rdata")}else{
   load("Risksample/negsample.Rdata")
    }

if(screen_strategy==1 | screen_strategy==2 | screen_strategy==7 | screen_strategy==8 | screen_strategy==9){

#Retain key data
negsample<-data.frame("risk_group"=negsample$risk_group,
                      "MRI_screen"=negsample$MRI_screen,
                      "US_screen"=negsample$US_screen,
                      "risk_predicted"=negsample$risk_predicted,
                      "feedback"=negsample$feedback,
                      "interval_change"=negsample$interval_change,
                      "life_expectancy"=negsample$life_expectancy)

negsample$risk_group<-negsample$risk_group*negsample$interval_change

subsamples<-unique(negsample$risk_group)
mastersample<-negsample

for (i in 1:length(subsamples)){
  negsample<-mastersample %>% filter(mastersample$risk_group==subsamples[i])

#Set screen times
if(negsample$risk_group[1]==5){
  screen_times<-high_risk_screentimes
}
if(negsample$risk_group[1]==4){
  screen_times<-med_risk_screentimes
}
if(negsample$risk_group[1]==3 | negsample$risk_group[1]==2 | negsample$risk_group[1]==0){
  screen_times<-low_risk_screentimes
}
if(negsample$risk_group[1]==1 & (screen_strategy==9 | screen_strategy==7)){
  screen_times<-seq(screen_startage, screen_startage+(5*4),5)
}
  if(negsample$risk_group[1]==1 & screen_strategy==8){
    screen_times<-seq(screen_startage,screen_startage+(6*3),6)
  }

#Add blank columns for potential screen times
for(i in 1:length(screen_times)){
  negsample[,7+i]<-numeric(length(negsample$negsample.risk_group))
}

#Draw attendance at first screen
negsample[,8]<-rbinom(length(negsample$risk_group),1,uptakefirstscreen)

#Loop through remaining screens conditional on previous attendance
for (i in 1:length(screen_times)-1){
negsample[,8+i]<-ifelse(rowSums(negsample[8:(7+i)])>=1,
                      rbinom(length(negsample$risk_group),1,uptakeotherscreen),
                      rbinom(length(negsample$risk_group),1,uptakenoscreen))
}

#Remove screening attendance after death
for (i in 1:length(screen_times)){
  negsample[,7+i][negsample$life_expectancy<rep(screen_times[i],length(negsample$life_expectancy))]<-0
}

#Calculate screens attended
negsample$total_screens<-rowSums(negsample[8:length(negsample[1,])])

#Calculate screening cost
for (i in 1:length(screen_times)){
  negsample[,7+i]<-negsample[,7+i]*(cost_screen*((1/((1+discount_cost)^(screen_times[i]-screen_startage)))))
}

negsample <- negsample %>%
  mutate(first_case = {
    tmp <- select(., starts_with('V'))
    ifelse(rowSums(tmp) == 0, NA, max.col(tmp != 0, ties.method = 'first'))
  })

negsample$screencost<-rowSums(negsample[8:length(negsample[1,])])
negsample$riskcost<-rep(cost_strat,length=nrow(negsample))*
  ((1/((1+discount_cost)^((screen_times[negsample$first_case]-rep(screen_startage,length(nrow(negsample))))))))
negsample$screencost<-negsample$screencost+negsample$riskcost

#Create QALY vector
negsample$QALY<-rep(0,length=length(negsample$risk_group))

#Create utility weight lookup table
qalylookup<-data.frame("age"=seq(from=screen_startage,to=100,by=1),
                       "qalyweight"=rep(0,length=100-screen_startage+1))

#Fill in utility values for each age with discounting
for (i in 1:length(qalylookup$qalyweight)){
qalylookup$qalyweight[i]<-utility_ages[match((ceiling(((screen_startage-1)+i)/5)*5),utility_ages[,1]),2]*(1/(1+discount_health)^i)
}

#Calculate cumulative QALYs for round ages
qalylookup$qalyyear<-cumsum(qalylookup$qalyweight)

negsample$QALY<-(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),3])+
  ((negsample$life_expectancy-floor(negsample$life_expectancy))*(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),2]))

if(screen_strategy==0 | screen_strategy>9){
  negsample$screencost<-rep(0,length=nrows(negsample))
}

results<-data.frame(negsample$QALY,
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

save(results,file = paste(det_output_path,
                          "Determ_",
                          screen_strategy,
                          "_",
                          negsample$risk_group[1],
                          "negresults",
                          ".Rdata",
                          sep = ""))
}
}else{

#Set screen times
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

#Add blank columns for potential screen times
for(i in 1:length(screen_times)){
  negsample[,7+i]<-numeric(length(negsample$negsample.risk_group))
}

#Draw attendance at first screen
negsample[,8]<-rbinom(length(negsample$risk_group),1,uptakefirstscreen)

#Loop through remaining screens conditional on previous attendance
for (i in 1:length(screen_times)-1){
  negsample[,8+i]<-ifelse(rowSums(negsample[8:(7+i)])>=1,
                          rbinom(length(negsample$risk_group),1,uptakeotherscreen),
                          rbinom(length(negsample$risk_group),1,uptakenoscreen))
}

#Remove screening attendance after death
for (i in 1:length(screen_times)){
  negsample[,7+i][negsample$life_expectancy<rep(screen_times[i],length(negsample$life_expectancy))]<-0
}

#Calculate screens attended
negsample$total_screens<-rowSums(negsample[8:length(negsample[1,])])

#Calculate screening cost
for (i in 1:length(screen_times)){
  negsample[,7+i]<-negsample[,7+i]*(cost_screen*((1/((1+discount_cost)^(screen_times[i]-screen_startage)))))
}
negsample$screencost<-rowSums(negsample[8:length(negsample[1,])])

#Create QALY vector
negsample$QALY<-rep(0,length=length(negsample$risk_group))

#Create utility weight lookup table
qalylookup<-data.frame("age"=seq(from=screen_startage,to=100,by=1),
                       "qalyweight"=rep(0,length=100-screen_startage+1))

#Fill in utility values for each age with discounting
for (i in 1:length(qalylookup$qalyweight)){
  qalylookup$qalyweight[i]<-utility_ages[match((ceiling(((screen_startage-1)+i)/5)*5),utility_ages[,1]),2]*(1/(1+discount_health)^i)
}

#Calculate cumulative QALYs for round ages
qalylookup$qalyyear<-cumsum(qalylookup$qalyweight)

negsample$QALY<-(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),3])+
  ((negsample$life_expectancy-floor(negsample$life_expectancy))*(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),2]))

if(screen_strategy==0 | screen_strategy>9){
  negsample$screencost<-rep(0,length=nrows(negsample))
}

results<-data.frame(negsample$QALY,
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

save(results,file = paste(det_output_path,
                          "Determ_",
                          screen_strategy,
                          "_",
                          "negresults",
                          ".Rdata",
                          sep = ""))
}
}
