library(tidyverse)
library(magrittr)
library(readr)
library(MASS, exclude = 'select')
library(EnvStats)

fnLogit <- function(x) {log(x / (1-x))}
fnALogit <- function(x) {exp(x) / (1+exp(x))}
fnBoxCox <- function(y) {
  xx <- EnvStats::boxcox(y, objective.name = "PPCC", optimize = TRUE)
  return(xx$lambda)
}

tblSynth <- read_csv("synthetic_risk_data.csv")

tblSynth %<>%
  rename_with(~str_remove(., "syn.")) %>%
  rename(ID = 1,
         VBD = VBD.) %>%
  mutate(across(-c(ID, Age), ~./100))

tblS <- tblSynth %>%
  pivot_longer(cols = -c(ID, Age),
               names_to = "var",
               values_to = "raw") %>%
  mutate(log = log(raw),
         logit = fnLogit(raw)) %>%
  group_by(var) %>%
  mutate(lambda = fnBoxCox(raw),
         BC = (raw^lambda-1)/lambda)

tblBC <- tblS %>%
  mutate(var = paste0("BC_", var)) %>%
  distinct(var, lambda)

tblS %<>%
  pivot_wider(id_cols = c(ID, Age),
              names_from = "var",
              values_from = c(raw, log, logit, BC))

m <- tblS %>%
  select(Age, contains("_")) %>%
  summarise(across(everything(), mean)) %>%
  as.numeric()

v <- tblS %>%
  select(Age, contains("_")) %>%
  cov() %>%
  as.matrix()

nSample <- NROW(tblSynth)
tblMV <- mvrnorm(n     = nSample,
                 mu    = m,
                 Sigma = v) %>%
  as_tibble() %>%
  pivot_longer(cols         = starts_with("BC_"),
               names_to     = "var",
               values_to    = "val") %>%
  left_join(tblBC) %>%
  mutate(val = exp(log(lambda*val+1)/lambda)) %>%
  select(-lambda) %>%
  pivot_wider(names_from  = "var",
              values_from = "val") %>%
  mutate(across(starts_with("log_"), exp),
         across(starts_with("logit_"), fnALogit))

summary(tblSynth)
summary(tblMV)

tblMV %<>%
  pivot_longer(cols      = -Age,
               names_to  = c("trans", "var"),
               names_sep = "_",
               values_to = "val") %>%
  mutate(trans = ifelse(trans=="raw", "norm", trans)) %>%
  bind_rows(tblS %>%
              select(Age, starts_with("raw_")) %>%
              pivot_longer(cols      = -Age,
                           names_to  = c("trans", "var"),
                           names_sep = "_",
                           values_to = "val") %>%
              select(Age, var, val) %>%
              mutate(trans = "orig"))

tblGG <- tblMV %>%
  mutate(bin    = cut(val,
                      breaks = 0:100 / 100,
                      labels = 0:99 / 100,
                      right = F,
                      include.lowest = T),
         binVal = (as.numeric(bin)-1) * 0.01) %>%
  select(-bin) %>%
  group_by(trans, var, binVal) %>%
  summarise(n = n()) %>%
  mutate(c = n/sum(n)) %>%
  filter(binVal<0.2) %>%
  ungroup()

tblGG %>%
  filter(trans != "orig") %>%
  ggplot() +
  geom_col(data     = tblGG %>% filter(trans == "orig"),
           aes(x=binVal, y=c),
           colour   = "black",
           fill     = NA) +
  geom_col(aes(x=binVal, y=c, fill=trans),
           width    = 0.7,
           position = position_dodge2(width = 1)) +
  facet_wrap(vars(var),
             nrow   = 3,
             scales = "free") +
  scale_x_binned(breaks = 0:20 / 100,
                 right  = F,
                 limits = c(0,0.2),
                 expand = c(0,0)) +
  theme_bw()

# tblMV %>%
#   filter(trans != "orig") %>%
#   ggplot() +
#   geom_freqpoly(data     = tblMV %>% filter(trans == "orig"),
#            aes(x=val),
#            colour   = "black",
#            fill     = "red") +
#   geom_freqpoly(aes(x=val, colour=trans), binwidth=0.01) +
#   facet_wrap(vars(var),
#              nrow   = 3,
#              scales = "free")  +
#   scale_x_continuous(expand       = c(0, 0),
#                      limits       = c(-0.1,0.2))
#   
# 
# tblMV %>%
#   filter(trans != "orig") %>%
#   ggplot() +
#   geom_histogram(data     = tblMV %>%
#                    filter(trans == "orig"),
#                  aes(x=val, fill=trans),
#                  colour   = "black",
#                  fill     = NA,
#                  position = position_dodge2(width = 0.9, preserve = "single")) +
#   geom_histogram(aes(x=val, fill=trans),
#                  position = position_dodge2(width = 0.9, preserve = "single")) +
#   facet_wrap(vars(var),
#              nrow = 3,
#              scales = "free") +
#   theme_bw()