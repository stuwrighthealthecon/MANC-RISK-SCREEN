#matrix with proporiton in each ?NPI group for each size category
#Updated 1010 include DCIS
NPI_by_size <- function(Ca_size,screen_detected_ca){
  NPI_cat <- 0
  if (runif(1,0,1)<0.21 && screen_detected_ca == 1){
    NPI_cat <-4
  }
#first determine if advanced cancer or not based on metastatic prob by size (categorical)
    if(Ca_size<= 25){m_size <- 25}else{m_size <- ceiling((Ca_size-25)/10)*10+25}
    if (m_size > 85){m_size <- 85}
    if(runif(1,0,1) < metastatic_prob[match(m_size, metastatic_prob[,1]),2] && NPI_cat == 0){NPI_cat <- 5} 
  
    #sample from categories 1,2 & 3 with probability of each based on the correct row of NPIbysize matrix 
    # Ca_size is continuous, need to match to closest larger value in NPI_by_size column 1
  if(NPI_cat == 0){
    v <- c(0.025, 5, 10, 15, 20, 30, 128) #category cut-points from Kolias 1999
  size_cat <- findInterval(Ca_size,v)
  NPI_cat <- sample(x=c(1,2,3),size = 1,prob = c(NPI_by_size_mat[size_cat,])) #1 best 3 worst prognosis
  }
#return the NPI category
result <- NPI_cat

  return(result)
}