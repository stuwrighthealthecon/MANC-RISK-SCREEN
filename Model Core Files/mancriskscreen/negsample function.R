##############################Function for estaimating outcomes for non-cancer############

negsamplefn<-function(screen_strategy,MISCLASS,PSA){
  
#Load appropriate data
  if(MISCLASS){
    load("Risksamplewithmisclass/negsample.Rdata")}else{
   load("Risksample/negsample.Rdata")
    }
  if(MISCLASS){
    #Assign women to risk groups based on 10yr risk if using risk-stratified approach  
    if(screen_strategy==1 | screen_strategy==9) {
      negsample$risk_group<-1+findInterval(negsample$tenyrrisk_est,risk_cutoffs_procas)
    } else
      if(screen_strategy==2) {
        negsample$risk_group<-1+findInterval(negsample$tenyrrisk_est,risk_cutoffs_tert)
        negsample$risk_group<-negsample$risk_group+2
      } else
        if(screen_strategy==7 | screen_strategy==8) {
          negsample$risk_group<-ifelse(negsample$tenyrrisk_est<low_risk_cut,1,2)
        } }
  
  if(MISCLASS==FALSE){
          if(screen_strategy==1 | screen_strategy==9) {
            negsample$risk_group<-1+findInterval(negsample$tenyrrisk,risk_cutoffs_procas)
          } else
            if(screen_strategy==2) {
              negsample$risk_group<-1+findInterval(negsample$tenyrrisk,risk_cutoffs_tert)
              negsample$risk_group<-negsample$risk_group+2
            } else
              if(screen_strategy==7 | screen_strategy==8) {
                negsample$risk_group<-ifelse(negsample$tenyrrisk<low_risk_cut,1,2) 
              }}
if(screen_strategy==1 | screen_strategy==2 | screen_strategy==7 | screen_strategy==8 | screen_strategy==9){

#Retain key data including PSA values if PSA used
  if(PSA==1){
    savePSA<-as.data.frame(negsample[c("risk_group","interval_change","PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                                       "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                                       "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                                       'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                                       "PSA_log_norm_mean","PSA_log_norm_sd","PSA_eff_ana","PSA_eff_tam",
                                       "PSA_dropout_ana","PSA_dropout_tam","PSA_uptake_1","PSA_uptake_2",
                                       "PSA_cost_strat","PSA_costvar",
                                       "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                                       "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","PSA_cost_drug","mcid")])
  }
  
negsample<-data.frame("risk_group"=negsample$risk_group,
                      "MRI_screen"=negsample$MRI_screen,
                      "US_screen"=negsample$US_screen,
                      "risk_predicted"=negsample$risk_predicted,
                      "feedback"=negsample$feedback,
                      "interval_change"=negsample$interval_change,
                      "life_expectancy"=negsample$life_expectancy,
                      "cost_screen"=cost_screen_base,
                      "cost_strat"=cost_strat,
                      "cost_follow_up"=cost_follow_up_base,
                      "cost_biop"=cost_biop_base)

#Overwrite deterministic input values for PSA
if(PSA==1){
  negsample$cost_screen<-(1+savePSA$PSA_costscreen)*cost_screen_base
  negsample$cost_strat<-savePSA$PSA_cost_strat
  negsample$cost_follow_up<-(1+savePSA$PSA_cost_follow_up)*cost_follow_up_base
  negsample$cost_biop<-(1+savePSA$PSA_cost_biop)*cost_biop_base
}

negsample$risk_group<-negsample$risk_group*negsample$interval_change

if(PSA==1){
savePSA$risk_group<-savePSA$risk_group*savePSA$interval_change
}

#Create an index of risk-gorups to cycle over for stratified screening programmes
subsamples<-unique(negsample$risk_group)
mastersample<-negsample
rm(negsample)

#Iterate over risk-groups, loading in all individuals from each group
for (ii in 1:length(subsamples)){
  negsample<-filter(mastersample,risk_group==subsamples[ii])
  if(PSA==1){
    #Load in PSA values for the group
subPSA<-filter(savePSA,risk_group==subsamples[ii])
subPSA<-subPSA[,-c(1:2)]}
  
#Set screen times
screen_times<-low_risk_screentimes
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
  negsample[,11+i]<-numeric(length(negsample$negsample.risk_group))
}

#Draw attendance at first screen
negsample[,12]<-rbinom(length(negsample$risk_group),1,uptakefirstscreen)

#Loop through remaining screens conditional on previous attendance
for (i in 1:(length(screen_times)-1)){
negsample[,12+i]<-ifelse(rowSums(negsample[12:(11+i)])>=1,
                      rbinom(length(negsample$risk_group),1,uptakeotherscreen),
                      rbinom(length(negsample$risk_group),1,uptakenoscreen))
}

#Remove screening attendance after death
for (i in 1:length(screen_times)){
  negsample[,11+i][negsample$life_expectancy<rep(screen_times[i],length(negsample$life_expectancy))]<-0
}

#Calculate screens attended
negsample$total_screens<-rowSums(negsample[12:length(negsample[1,])])

#Calculate screening cost
for (i in 1:length(screen_times)){
  negsample[,11+i]<-negsample[,11+i]*(negsample$cost_screen+
                                      (recall_rate*negsample$cost_follow_up)+
                                      (recall_rate*biopsy_rate*negsample$cost_biop)+
                                        (negsample$MRI_screen*cost_MRI)+
                                        (negsample$US_screen*cost_US)*
                                      ((1/((1+discount_cost)^(screen_times[i]-screen_startage)))))
}

#Find first screening event to add risk prediciton cost
negsample <- negsample %>%
  mutate(first_case = {
    tmp <- dplyr::select(negsample,starts_with('V'))
    ifelse(rowSums(tmp) == 0, NA, max.col(tmp != 0, ties.method = 'first'))
  })

#Add up screening costs and risk prediciton costs
negsample$screencost<-rowSums(negsample[12:length(negsample[1,])])
negsample$riskcost<-rep(negsample$cost_strat,length=nrow(negsample))*
  ((1/((1+discount_cost)^((screen_times[negsample$first_case]-rep(screen_startage,length(nrow(negsample))))))))
negsample$screencost<-negsample$screencost+negsample$riskcost

#Create QALY vector
negsample$QALY<-rep(0,length=length(negsample$risk_group))

#Create QALY weight lookup table
qalylookup<-data.frame("age"=seq(from=screen_startage,to=100,by=1),
                       "qalyweight"=rep(0,length=100-screen_startage+1))

#Fill in utility values for each age with discounting
for (i in 1:length(qalylookup$qalyweight)){
  qalylookup$qalyweight[i]<-utility_ages[match((ceiling((screen_startage-1)+i)),utility_ages[,1]),2]*(1/(1+discount_health)^(i))
}

#Calculate cumulative QALYs for round ages
qalylookup$qalyyear<-cumsum(qalylookup$qalyweight)

#Adjust QALYs for partial years at end of life
negsample$QALY<-(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),3])+
  ((negsample$life_expectancy-floor(negsample$life_expectancy))*((qalylookup[match(ceiling(negsample$life_expectancy),qalylookup[,1]),2])-
                                                                   (qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),2]))
  )

#Create results table
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

#Add 0 cost if na found in table due to no screening events
results$Cost[is.na(results$Cost)]<-0

#Bind PSA values to results if PSA 
if(PSA==1){
results<-cbind(results,subPSA)
}

#Save outputs
save(results,file = paste(ifelse(PSA==0,det_output_path,psa_output_path),
                          screen_strategy,
                          "_",
                          ii,
                          "_",
                          "negresults",
                          ".Rdata",
                          sep = ""))
}
}

  if(screen_strategy>2 & screen_strategy<7){
    
    #Save PSA values
    if(PSA==1){
      savePSA<-as.data.frame(negsample[c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                                         "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                                         "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                                         'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                                         "PSA_log_norm_mean","PSA_log_norm_sd","PSA_eff_ana","PSA_eff_tam",
                                         "PSA_dropout_ana","PSA_dropout_tam","PSA_uptake_1","PSA_uptake_2",
                                         "PSA_cost_strat","PSA_costvar",
                                         "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                                         "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","PSA_cost_drug","mcid")])
    }
    
    #Retain key values in negsample
    negsample<-data.frame("risk_group"=negsample$risk_group,
                          "MRI_screen"=negsample$MRI_screen,
                          "US_screen"=negsample$US_screen,
                          "risk_predicted"=negsample$risk_predicted,
                          "feedback"=negsample$feedback,
                          "interval_change"=negsample$interval_change,
                          "life_expectancy"=negsample$life_expectancy,
                          "cost_screen"=cost_screen_base,
                          "cost_strat"=cost_strat,
                          "cost_follow_up"=cost_follow_up_base,
                          "cost_biop"=cost_biop_base)
    
    #if PSA then add in PSA input values
    if(PSA==1){
      negsample$cost_screen<-(1+savePSA$PSA_costscreen)*cost_screen_base
      negsample$cost_strat<-savePSA$PSA_cost_strat
      negsample$cost_follow_up<-(1+savePSA$PSA_cost_follow_up)*cost_follow_up_base
      negsample$cost_biop<-(1+savePSA$PSA_cost_biop)*cost_biop_base
    }
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
  negsample[,11+i]<-numeric(length(negsample$negsample.risk_group))
}

#Draw attendance at first screen
negsample[,12]<-rbinom(length(negsample$risk_group),1,uptakefirstscreen)

#Loop through remaining screens conditional on previous attendance
for (i in 1:(length(screen_times)-1)){
  negsample[,12+i]<-ifelse(rowSums(negsample[12:(11+i)])>=1,
                          rbinom(length(negsample$risk_group),1,uptakeotherscreen),
                          rbinom(length(negsample$risk_group),1,uptakenoscreen))
}

#Remove screening attendance after death
for (i in 1:length(screen_times)){
  negsample[,11+i][negsample$life_expectancy<rep(screen_times[i],length(negsample$life_expectancy))]<-0
}

#Calculate screens attended
negsample$total_screens<-rowSums(negsample[12:length(negsample[1,])])

#Calculate screening cost
for (i in 1:length(screen_times)){
  negsample[,11+i]<-negsample[,11+i]*(negsample$cost_screen+
                                        (recall_rate*negsample$cost_follow_up)+
                                        (recall_rate*biopsy_rate*negsample$cost_biop)+
                                        (negsample$MRI_screen*cost_MRI)+
                                        (negsample$US_screen*cost_US)*
                                        ((1/((1+discount_cost)^(screen_times[i]-screen_startage)))))
}
negsample$screencost<-rowSums(negsample[12:length(negsample[1,])])

#Create QALY vector
negsample$QALY<-rep(0,length=length(negsample$risk_group))

#Create utility weight lookup table
qalylookup<-data.frame("age"=seq(from=screen_startage,to=100,by=1),
                       "qalyweight"=rep(0,length=100-screen_startage+1))

#Fill in utility values for each age with discounting
for (i in 1:length(qalylookup$qalyweight)){
  qalylookup$qalyweight[i]<-utility_ages[match((ceiling((screen_startage-1)+i)),utility_ages[,1]),2]*(1/(1+discount_health)^(i))
}

#Calculate cumulative QALYs for round ages
qalylookup$qalyyear<-cumsum(qalylookup$qalyweight)

#Adjust QALY weights for partial year in last year of life
negsample$QALY<-(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),3])+
  ((negsample$life_expectancy-floor(negsample$life_expectancy))*((qalylookup[match(ceiling(negsample$life_expectancy),qalylookup[,1]),2])-
                                                                   (qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),2]))
  )

#Create results frame
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

#Replace na with 0 cost in event of no screens
results$Cost[is.na(results$Cost)]<-0

#If PSA bind results and PSA values
if(PSA==1){
  results<-cbind(results,savePSA)
}

#Save outputs
save(results,file = paste(ifelse(PSA==0,det_output_path,psa_output_path),
                          screen_strategy,
                          "_",
                          ii,
                          "_",
                          "negresults",
                          ".Rdata",
                          sep = ""))

  }
  if (screen_strategy==0 | screen_strategy>9){
    
    #Save PSA values
    if(PSA==1){
      savePSA<-as.data.frame(negsample[c("PSA_gamma_survival_1","PSA_gamma_survival_2","PSA_gamma_survival_3",
                                         "PSA_meta_survival_54","PSA_meta_survival_74","PSA_meta_survival_99",
                                         "PSA_beta_1","PSA_beta_2",'PSA_VDG1_sen','PSA_VDG2_sen',
                                         'PSA_VDG3_sen', 'PSA_VDG4_sen',"PSA_MRI_cdr","PSA_US_cdr",
                                         "PSA_log_norm_mean","PSA_log_norm_sd","PSA_eff_ana","PSA_eff_tam",
                                         "PSA_dropout_ana","PSA_dropout_tam","PSA_uptake_1","PSA_uptake_2",
                                         "PSA_cost_strat","PSA_costvar",
                                         "PSA_util_1to3","PSA_util_4","PSA_costscreen","PSA_cost_follow_up",
                                         "PSA_cost_biop","PSA_cost_US","PSA_cost_MRI","PSA_cost_drug","mcid")])
    }
    
#Retain key data from negsample
    negsample<-data.frame("risk_group"=negsample$risk_group,
                          "MRI_screen"=negsample$MRI_screen,
                          "US_screen"=negsample$US_screen,
                          "risk_predicted"=negsample$risk_predicted,
                          "feedback"=negsample$feedback,
                          "interval_change"=negsample$interval_change,
                          "life_expectancy"=negsample$life_expectancy,
                          "cost_screen"=cost_screen_base,
                          "cost_strat"=cost_strat,
                          "cost_follow_up"=cost_follow_up_base,
                          "cost_biop"=cost_biop_base)
    
    #If PSA then add in PSAV input values
    if(PSA==1){
      negsample$cost_screen<-(1+savePSA$PSA_costscreen)*cost_screen_base
      negsample$cost_strat<-savePSA$PSA_cost_strat
      negsample$cost_follow_up<-(1+savePSA$PSA_cost_follow_up)*cost_follow_up_base
      negsample$cost_biop<-(1+savePSA$PSA_cost_biop)*cost_biop_base
    }
    
    #Create vector of 0 screen costs as no screening in strategy
    negsample$screencost<-rep(0,length=nrow(negsample))
    negsample$total_screens<-rep(0,length=nrow(negsample))
    
    #Create QALY vector
    negsample$QALY<-rep(0,length=length(negsample$risk_group))
    
    #Create utility weight lookup table
    qalylookup<-data.frame("age"=seq(from=screen_startage,to=100,by=1),
                           "qalyweight"=rep(0,length=100-screen_startage+1))
    
    #Fill in utility values for each age with discounting
    for (i in 1:length(qalylookup$qalyweight)){
      qalylookup$qalyweight[i]<-utility_ages[match((ceiling((screen_startage-1)+i)),utility_ages[,1]),2]*(1/(1+discount_health)^(i))
    }
    
    #Calculate cumulative QALYs for round ages
    qalylookup$qalyyear<-cumsum(qalylookup$qalyweight)
    
    #Adjust QALYs for partial year in last year of life
    negsample$QALY<-(qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),3])+
      ((negsample$life_expectancy-floor(negsample$life_expectancy))*((qalylookup[match(ceiling(negsample$life_expectancy),qalylookup[,1]),2])-
        (qalylookup[match(floor(negsample$life_expectancy),qalylookup[,1]),2]))
      )

    #Save results
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

#Replace na with 0 where no screening events
results$Cost[is.na(results$Cost)]<-0

#If PSA bind results and PSA values
if(PSA==1){
  results<-cbind(results,savePSA)
}

#Save outputs
save(results,file = paste(ifelse(PSA==0,det_output_path,psa_output_path),
                          screen_strategy,
                          "_",
                          ii,
                          "_",
                          "negresults",
                          ".Rdata",
                          sep = ""))
}
}