library("mgcv")
library("dampack")
library(gt)
library(tidyverse)
library(magrittr)

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
input_df<-data.frame(matrix(nrow=6,ncol=23))
input_df[1,]<-as.numeric(input_vector)
input_df[2,]<-as.numeric(input_vector)
input_df[3,]<-as.numeric(input_vector)
input_df[4,]<-as.numeric(input_vector)
input_df[5,]<-as.numeric(input_vector)
input_df[6,]<-as.numeric(input_vector)

strategies <- as.factor(levels(modQ$var.summary$alternative))
input_df[,24] <- data.frame(alternative = strategies)
names(input_df)<-c(names(input_vector),"alternative")

#Predict the QALYs and costs for each strategy
#This is not giving sensible values at the moment so I need to fix
#Should be ok for now
#Predict the QALYs and costs for each strategy
qaly <- mgcv::predict.bam(modQ,input_df)
cost <- mgcv::predict.bam(modC,input_df)

#Simple output table and graph [REQUIRES DAMPACK]
output_df <- data.frame(qualy = qaly, cost = cost)
rownames(output_df) <- strategies
output_df[,"incQALYS"]<-c(output_df$qualy-output_df$qualy[4])
output_df[,"incCost"]<-c(output_df$cost-output_df$cost[4])
output_df[,"ICER"]<-c(output_df$incCost/output_df$incQALYS)
output_df[,"NB20k"]<-c((output_df$incQALYS*20000)-output_df$incCost)
output_df[,"NB30k"]<-c((output_df$incQALYS*30000)-output_df$incCost)

icer_strat<-calculate_icers(cost=output_df$cost,
                            effect=output_df$qualy,
                            strategies = c(row.names(output_df)))
plot(icer_strat,currency="£",label="all")

#Complex but nicer table and graph [REQUIRES gt, tidyverse, magrittr]
fnIncCU <- function(Names, Costs, QALYs, blnNHB=TRUE, WTP=c(20000,30000), blnCUP=TRUE, costDP = 2, qalyDP = 5, icerDP = 2) {
  nn <- NROW(Names)
  IncCU <- tibble(StratName  = Names,
                  Cost       = Costs,
                  QALY       = QALYs,
                  dom        = FALSE,
                  extdom     = FALSE,
                  comp       = 1,
                  IncCost    = 0,
                  IncQALY    = 0,
                  ICER       = 0,
                  strICER    = "-",
                  strIncCost = "-",
                  strIncQALY = "-",
                  strCost    = format(round(Cost, costDP), nsmall = costDP, big.mark=","),
                  strQALY    = format(round(QALY, qalyDP), nsmall = qalyDP, big.mark=",")) %>%
    arrange(Cost) %>%
    mutate(ID = row_number())
  for (i in 2:nn) {
    comp <- max(which(!IncCU$dom[1:(i-1)] & !IncCU$extdom[1:(i-1)]))
    IncCU$comp[i]       <- comp
    IncCU$IncCost[i]    <- IncCU$Cost[i]-IncCU$Cost[comp]
    IncCU$IncQALY[i]    <- IncCU$QALY[i]-IncCU$QALY[comp]
    IncCU$ICER[i]       <- IncCU$IncCost[i] / IncCU$IncQALY[i]
    IncCU$dom[i]        <- max(IncCU$QALY[1:(i-1)]) > IncCU$QALY[i]
    IncCU$extdom[i]     <- if (i==nn) FALSE else any(IncCU$QALY[(i+1):nn]>IncCU$QALY[i] & (IncCU$Cost[(i+1):nn]-IncCU$Cost[comp])/(IncCU$QALY[(i+1):nn]-IncCU$QALY[comp]) < IncCU$ICER[i])
    IncCU$strICER[i]    <- ifelse(IncCU$dom[i], "dominated", ifelse(IncCU$extdom[i], "extendedly dominated", format(round(IncCU$ICER[i], icerDP), nsmall = icerDP, big.mark=",")))
    IncCU$strIncCost[i] <- format(round(IncCU$IncCost[i], costDP), nsmall = costDP, big.mark=",")
    IncCU$strIncQALY[i] <- format(round(IncCU$IncQALY[i], qalyDP), nsmall = qalyDP, big.mark=",")
  }
  if (blnNHB) {
    IncCU %<>%
      rowwise() %>%
      mutate(NHB = map2(Cost, QALY, ~ (QALY - Cost / WTP) %>%
                          setNames(paste0("NHB", WTP)))) %>%
      unnest_wider(NHB) %>%
      mutate(across(starts_with("NHB"), ~ rank(desc(.x)), .names = "rank{.col}"))
  }
  
  if (blnCUP) {
    or.x <- IncCU$QALY[which(IncCU$Cost==min(IncCU$Cost))]
    or.y <- min(IncCU$Cost)
    ax.x <- fnAxisScale(min(IncCU$QALY), max(IncCU$QALY))
    ax.y <- fnAxisScale(or.y, max(IncCU$Cost))
    tl.x <- (ax.x$dMax - ax.x$dMin) / 50
    tl.y <- (ax.y$dMax - ax.y$dMin) / 50
    nt.x <- ((ax.x$dMax - ax.x$dMin) / ax.x$dMajor + 1) %>% round(0)
    nt.y <- ((ax.y$dMax - ax.y$dMin) / ax.y$dMajor + 1) %>% round(0)
    tblXT <- tibble(x  = ax.x$dMin + (1:nt.x-1)*ax.x$dMajor,
                    y1 = or.y - tl.y/2,
                    y2 = or.y + tl.y/2)
    tblYT <- tibble(y  = ax.y$dMin + (1:nt.y-1)*ax.y$dMajor,
                    x1 = or.x - tl.x/2,
                    x2 = or.x + tl.x/2)
    
    plt <- IncCU %>%
      ggplot(aes(x=QALY, y=Cost)) +
      geom_abline(slope     = WTP[1],
                  intercept = seq(ax.y$dMin - ax.y$dMajor * 100,
                                  ax.y$dMax + ax.y$dMajor * 100,
                                  ax.y$dMajor) - or.x * WTP[1],
                  colour    = "white",
                  lty       = 5,
                  linewidth      = 1) +
      geom_hline(yintercept = or.y,
                 linewidth       = 1.25) +
      geom_segment(data = tblXT,
                   aes(x    = x,
                       xend = x,
                       y    = y1,
                       yend = y2),
                   linewidth = 1.25) +
      geom_vline(xintercept = or.x, linewidth = 1.25) +
      geom_segment(data = tblYT, 
                   aes(y    = y,
                       yend = y,
                       x    = x1,
                       xend = x2),
                   linewidth = 1.25) +
      geom_line(data = IncCU %>%
                  filter(str_detect(strICER, "dominated", negate = T)),
                aes(x = QALY,
                    y = Cost),
                colour = "red", linewidth = 1.25) +
      geom_point(shape = 19, colour = "darkgreen", size = 10) +
      geom_text(aes(label = ID),
                colour = "white", size = 12 * 1/72 * 25.4) +
      scale_x_continuous(limits = c(ax.x$dMin, ax.x$dMax*1.000000001),
                         breaks = seq(ax.x$dMin, ax.x$dMax, ax.x$dMajor),
                         labels = seq(ax.x$dMin, ax.x$dMax, ax.x$dMajor),
                         name   = "QALYs") +
      scale_y_continuous(limits = c(ax.y$dMin, ax.y$dMax*1.000000001),
                         breaks = seq(ax.y$dMin, ax.y$dMax, ax.y$dMajor),
                         labels = scales::number_format(prefix = "£",
                                                        big.mark = ","),
                         name   = "Costs\n") +
      theme(panel.background = element_rect(fill = "grey75"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 12),
            axis.title       = element_text(size = 12, face = "bold"),
            axis.ticks       = element_blank())
  } else plt <- NULL
  
  prettyIncCU <- IncCU %>%
    select(ID,
           StratName,
           strCost,
           strQALY,
           strIncCost,
           strIncQALY,
           strICER,
           contains("NHB")) %>%
    gt() %>%
    cols_label(
      StratName  = "Strategy",
      strCost    = "Costs\n(£)",
      strQALY    = "Effects\n(QALYs)",
      strIncCost = "Costs\n(£)",
      strIncQALY = "Effects\n(QALYs)",
      strICER    = "ICER\n(£/QALY)") %>%
    cols_label_with(
      starts_with("NHB"), ~ paste0("£", as.numeric(str_sub(.x, 4))/1000, "K/QALY")) %>%
    cols_label_with(
      starts_with("rankNHB"), ~ paste0("£", as.numeric(str_sub(.x, 8))/1000, "K/QALY")) %>%
    tab_spanner(label = "Incremental",
                columns = c(strIncCost, strIncQALY, strICER)) %>%
    tab_spanner(label = "Net health benefit",
                columns = starts_with("NHB")) %>%
    tab_spanner(label = "Rank",
                columns = starts_with("rankNHB")) %>%
    cols_align(align = "center",
               columns = -c(StratName)) %>%
    opt_horizontal_padding(scale = 2) %>%
    opt_table_font(font = "Arial")
  
  return(list(IncCU       = IncCU,
              prettyIncCU = prettyIncCU,
              CUplane     = plt))
}
fnAxisScale <- function(dMin, dMax) {
  ## adapted from some VBA of Jon Peltier's
  ## https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/  
  if (dMax == dMin) {
    dTemp = dMax
    dMax = dMax * 1.01
    dMin = dMin * 0.99
  }
  
  if (dMax < dMin) {
    dTemp = dMax
    dMax = dMin
    dMin = dTemp
  }
  
  if (dMax != 0) {
    if (dMax > 0) {
      dMax = dMax + (dMax - dMin) * 0.01
    }
    else {
      dMax = min(dMax + (dMax - dMin) * 0.01, 0)
    }
  }
  
  if (dMin != 0) {
    if (dMin > 0) {
      dMin = max(dMin - (dMax - dMin) * 0.01, 0)
    }
    else {
      dMin = dMin - (dMax - dMin) * 0.01
    }
  }
  
  if (dMax == 0 & dMin == 0) dMax = 1
  
  dPower = log(dMax - dMin) / log(10)
  dScale = 10 ^ (dPower - trunc(dPower))
  
  dScale <- case_when(
    between(dScale, 0, 2.5) ~ 0.2,
    between(dScale, 2.5, 5) ~ 0.5,
    between(dScale, 5, 7.5) ~ 1,
    TRUE                    ~ 2,
  )
  dSmall <- case_when(
    between(dScale, 0, 2.5) ~ 0.05,
    between(dScale, 2.5, 5) ~ 0.1,
    between(dScale, 5, 7.5) ~ 0.2,
    TRUE                    ~ 0.5,
  )
  
  dScale = dScale * 10 ^ trunc(dPower)
  dSmall = dSmall * 10 ^ trunc(dPower)
  
  list (dMin = dScale * trunc(dMin / dScale),
        dMax = dScale * (trunc(dMax / dScale) + 1),
        dMajor = dScale,
        dMinor = dSmall)
}

IncCU<-fnIncCU(Names=strategies,Costs=output_df$cost,QALYs=output_df$qualy)

#Formatted output table
IncCU$prettyIncCU

#Formatted cost-effectiveness plane
IncCU$CUplane
