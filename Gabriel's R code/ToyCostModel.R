library(tidyverse)

# read in and refactor data from Laudicella et al. (2016) https://doi.org/10.1038/bjc.2016.77
tbl <- tribble(~Yr, ~Early_18.64, ~Late_18.64, ~Diff1, ~Early_65plus, ~Late_65plus, ~Diff2,
               0, 464, 607, 143, 1086, 1324, 238,
               1, 10746, 13315, 2569, 7597, 8804, 1207,
               2, 3357, 5785, 2429, 2529, 3650, 1121,
               3, 1953, 3782, 1829, 2156, 3170, 1014,
               4, 1627, 2932, 1305, 2230, 2924, 693,
               5, 1617, 2841, 1225, 2077, 2957, 880,
               6, 1547, 2645, 1099, 2174, 2783, 609,
               7, 1394, 2618, 1225, 2063, 2903, 840,
               8, 1376, 2559, 1183, 2134, 2454, 320,
               9, 1279, 1848, 569, 2204, 2932, 728) %>%
  select(-Diff1, -Diff2) %>%
  pivot_longer(cols      = contains("6"),
               names_to  = c("Stage", "Age"),
               names_sep = "_",
               values_to = "Cost") %>%
  group_by(Stage, Age) %>%
  mutate(DCost      = Cost - first(Cost),
         DCost.i    = DCost * 1.219312579, # NHSCII inflator for 2010/11-->2020/21
         disc       = 1/1.035^(Yr-0.5),
         DCost.i.d  = DCost.i * disc,
         CDCost.i.d = cumsum(DCost.i.d),
         Yr1        = as.factor(Yr==1),
         Yr2        = as.factor(Yr==2),
         Yr3        = as.factor(Yr==3)) %>%
  filter(Yr > 0) %>%
  arrange(Stage, Age, Yr)

# log-linear model
mod <- lm(data = tbl,
          formula = log(DCost) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age)
mod %>% AIC()

# prediction matrix
tblNewDat <- crossing(Yr=1:50, Stage=c("Early", "Late"), Age=c("18.64", "65plus")) %>%
  mutate(Yr1 = as.factor(Yr==1),
         Yr2 = as.factor(Yr==2),
         Yr3 = as.factor(Yr==3))

# generate predictions
tblNewDat %>%
  bind_cols(pred = mod %>% predict(newdata = tblNewDat)) %>%
  mutate(DCost.p = exp(pred)) -> tblPred

# make lookup table
tblLookup <- tblPred %>%
  filter(Yr==1) %>%
  mutate(across(c(Yr, pred, DCost.p), ~0)) %>%
  bind_rows(tblPred) %>%
  group_by(Stage, Age) %>%
  mutate(DCost.p.i    = DCost.p * 1.219312579, # NHSCII inflator for 2010/11-->2020/21
         disc         = 1/1.035^(Yr-0.5),
         DCost.p.i.d  = DCost.p.i * disc,
         CDCost.p.i.d = cumsum(DCost.p.i.d),
         StageEarly   = Stage=="Early",
         AgeYoung     = Age=="18.64") %>%
  arrange(Stage, Age, Yr) %>%
  ungroup()

# visualise fitted model
tblPred %>%
  ggplot() +
  geom_line(aes(x=Yr, y=DCost.p, colour = interaction(Stage, Age)), size = 1, lty = 1) +
  geom_point(data = tbl, aes(x=Yr, y=DCost, colour = interaction(Stage, Age)), size = 3)

# visualise fitted model (cumulative, inflated, discounted costs)
tblLookup %>%
  filter(Yr>0) %>%
  ggplot() +
  geom_line(aes(x=Yr, y=CDCost.p.i.d, colour = interaction(Stage, Age)), size = 1, lty = 1) +
  geom_point(data = tbl, aes(x=Yr, y=CDCost.i.d, colour = interaction(Stage, Age)), size = 3)

# function to generate expected discounted lifetime costs, given age & stage categories and life expectancy
fnLookup_BCCosts <- function(blnStageEarly, blnAgeYoung, LE) {
  as.numeric(tblLookup$CDCost.p.i.d[tblLookup$StageEarly==blnStageEarly & tblLookup$AgeYoung==blnAgeYoung & tblLookup$Yr==floor(LE)]) +
    as.numeric(tblLookup$DCost.p.i.d[tblLookup$StageEarly==blnStageEarly & tblLookup$AgeYoung==blnAgeYoung & tblLookup$Yr==ceiling(LE)]) * (LE - floor(LE))
}
# example usage
fnLookup_BCCosts(blnStageEarly = T, blnAgeYoung = T, LE = 0.9999)


## ==== alternative approaches explored but turned out slower ==== 

modC <- lm(data = tbl,
          formula = (CDCost.i.d) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age)
modC %>% AIC()

tblNewDat <- crossing(Yr=1:50, Stage=c("Early", "Late"), Age=c("18.64", "65plus")) %>%
  mutate(Yr1 = as.factor(Yr==1),
         Yr2 = as.factor(Yr==2),
         Yr3 = as.factor(Yr==3))
tblNewDat %>%
  bind_cols(pred = modC %>% predict(newdata = tblNewDat)) %>%
  mutate(CDCost.i.d.p = (pred)) -> tblPredC

tblPredC %>%
  ggplot() +
  geom_line(aes(x=Yr, y=CDCost.i.d.p, colour = interaction(Stage, Age)), size = 1, lty = 1) +
  geom_point(data = tbl, aes(x=Yr, y=CDCost.i.d, colour = interaction(Stage, Age)), size = 3)


tblPredC %>%
  arrange(Stage, Age, Yr)

tblLookup %>%
  ggplot() +
  geom_line(aes(x=Yr, y=CDCost.p.i.d, colour = interaction(Stage, Age)), size = 1, lty = 1) +
  geom_line(data = tblPredC, aes(x=Yr, y=CDCost.i.d.p, colour = interaction(Stage, Age)), size = 1, lty = 2) +
  geom_point(data = tbl, aes(x=Yr, y=CDCost.i.d, colour = interaction(Stage, Age)), size = 3)

fnLookupTidy <- function(iStage, iAge, iLE) {
  tblLookup %>%
    filter(Stage==iStage, Age==iAge, Yr==iLE) %>%
    select(CDCost.p.i.d) %>%
    as.numeric()
}
fnLookupTidy("Early", "18.64", 12)

fnLookupBase <- function(iStage, iAge, iLE) {
  as.numeric(tblLookup$CDCost.p.i.d[tblLookup$Stage==iStage & tblLookup$Age==iAge & tblLookup$Yr==iLE])
}
fnLookupBase("Early", "18.64", 12)

library(data.table)
data_table <- data.table(tblLookup)
fnLookupDT <- function(iStage, iAge, iLE) {
  data_table[Stage==iStage & Age==iAge & Yr==iLE]$CDCost.p.i.d
}
fnLookupDT("Early", "18.64", 12)

setkey(data_table, Stage, Age, Yr)
fnLookupDTi <- function(iStage, iAge, iLE) {
  data_table[.(iStage, iAge, iLE), nomatch = 0L]$CDCost.p.i.d
}
fnLookupDTi("Early", "18.64", 12)

fnModPred <- function(iStage, iAge, iLE) {
  modC %>% predict(newdata = list(Stage=iStage, Age=iAge, Yr=iLE, Yr1 = as.factor(iLE==1),
                                  Yr2 = as.factor(iLE==2),
                                  Yr3 = as.factor(iLE==3))) %>%
    as.numeric()
}

fnModPred("Early", "18.64", 12)

fnModPred2 <- function(iStage, iAge, iLE) {
  as.numeric(predict(modC, newdata = list(Stage=iStage, Age=iAge, Yr=iLE, Yr1 = as.factor(iLE==1),
                                  Yr2 = as.factor(iLE==2),
                                  Yr3 = as.factor(iLE==3))))
}
fnModPred2("Early", "18.64", 12)

library(microbenchmark)
tblTst <- tibble(iLE    = sample(x = 1:50, size = 100, replace = TRUE),
                 iStage = sample(x = c("Early", "Late"), size = 100, replace = TRUE),
                 iAge   = sample(x = c("18.64", "65plus"), size = 100, replace = TRUE))
mbm <- microbenchmark(tblTst %>% pmap_dbl(fnLookupTidy),
                      tblTst %>% pmap_dbl(fnLookupBase),
                      tblTst %>% pmap_dbl(fnLookupDT),
                      tblTst %>% pmap_dbl(fnLookupDTi),
                      tblTst %>% pmap_dbl(fnModPred),
                      tblTst %>% pmap_dbl(fnModPred2))
mbm
autoplot(mbm)
