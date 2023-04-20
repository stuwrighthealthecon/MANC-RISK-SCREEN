library("mgcv")

#Load the predictive model
#This is quite large at the moment (1.8GB)
modQ<-readRDS("QALYmodelslim.RDS")
modC<-readRDS("costmodelslim.RDS")

#Input some values, these would be the values that decision makers can change
#For some we'd change it to a more logical value and then convert it to the data we need
#For example instead of gamma survival we'd have an input for 5 year cancer survival which we'd then convert in the R code
input_vector<-c("PSA_gamma_survival_1"= -5.462,
                "PSA_gamma_survival_2"= -3.814,
                "PSA_gamma_survival_3"= -2.723,
                "PSA_meta_survival_54"= -1.787,
                "PSA_meta_survival_74"= -1.388,
                "PSA_meta_survival_99"= -1.011,
                "PSA_beta_1"= 1.47,
                "PSA_beta_2"= 6.51,
                "PSA_VDG1_sen"= 0.85,
                "PSA_VDG2_sen"= 0.776,
                "PSA_VDG3_sen"= 0.695,
                "PSA_VDG4_sen"= 0.61,
                "PSA_log_norm_mean" = 1.07,
                "PSA_log_norm_sd" = 1.32,
                "PSA_cost_strat"= 8.12,
                "PSA_costvar" = 0,
                "PSA_util_1to3"= 0.82,
                "PSA_util_4"= 0.75,
                "PSA_costscreen" = 0.01,
                "PSA_cost_follow_up" = 0.01,
                "PSA_cost_biop"= 0.01,
                "PSA_cost_US" = 0.01,
                "PSA_cost_MRI" = 0.01
                )

#Make a row for each alternative (screening strategy)
#I'm sure there's a better way to do this!
input_vector_strategies<-data.frame(matrix(nrow=6,ncol=23))
input_vector_strategies[1,]<-as.numeric(input_vector)
input_vector_strategies[2,]<-as.numeric(input_vector)
input_vector_strategies[3,]<-as.numeric(input_vector)
input_vector_strategies[4,]<-as.numeric(input_vector)
input_vector_strategies[5,]<-as.numeric(input_vector)
input_vector_strategies[6,]<-as.numeric(input_vector)

input_vector_strategies[,24]<-as.factor(c("noscreening","tertiles","2yr","3yr","procas","fullstrat"))
names(input_vector_strategies)<-c(names(input_vector),"alternative")

#Predict the QALYs and costs for each strategy
#This is not giving sensible values at the moment so I need to fix
#Should be ok for now
predict.bam(modQ,input_vector_strategies)
predict.bam(modC,input_vector_strategies)


