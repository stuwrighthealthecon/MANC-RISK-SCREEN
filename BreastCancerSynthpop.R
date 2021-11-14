library(tidyverse)
library(magrittr)
library(synthpop)

## this code imagines your original data are in a tibble called tblIn
## I've used totally made-up variable names

## define predictor.matrix
## I've defined this example simply, with a list of predictor variables that are used to predict 
## all synthesised variables
## obviously, you could be more judicious about that
vars  <- c("gender", "age", "BRCA2", "density", "10yrRisk") # variables to synthesise
preds <- c("gender", "age", "BRCA2", "density") # predictor variables

## construct predictor matrix
## each row represents a variable to be synthesised; a 1 in any cell means the column variable will 
## be used as a predictor
mm <- matrix(0L,
             nrow = NCOL(tblIn),
             ncol = NCOL(tblIn),
             dimnames = rep(list(names(tblIn)), 2))
for (i in 1:length(vars)) {
  for (j in 1:length(preds)) {
    mm[vars[i],preds[j]] <- 1
  }
}

## define methods
meth <- rep("cart", NCOL(tblIn)) %>% # method for all vars starts off as CART
  set_names(names(tblIn))

## I had a time-to-event variable in my dataset; I don't think you have that complication
## but, if so, here's the syntax for defining one:
# meth["study_days"] <- "survctree" # different method for TTE variable
# ev <- list("study_days" = "OtherDeathStat") # identify event indicator to go with TTE variable

nSynth <- 2000 # number of 'people' to generate

## synthesis
synth <- tblIn %>%
  syn(predictor.matrix = mm,
      visit.sequence   = vars,
      method           = meth,
      k                = nSynth,
      # event            = ev, # reinsert this line if you have a TTE var as described above
      drop.not.used    = T,
      drop.pred.only   = T,
      minnumlevels     = 5)

## check results
tblSynth <- as_tibble(synth$syn)
summary(tblSynth)
compare(synth, tblIn)

## output results
write.csv(tblSynth, file = "synthpop.csv", row.names = FALSE)