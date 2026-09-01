library(tidyverse)

#Load data
filenames <- list.files(det_output_path, full.names = TRUE)
alldata <- lapply(filenames, function(x) {
  get(load(x, .GlobalEnv))
})
detresults <- as.data.frame(do.call("rbind", alldata))

#Create data.frame of screen and clinically detected cancers
screendet <- detresults %>% filter(detresults$`screen detected` == 1)
clindet <- detresults %>% filter(detresults$`screen detected` == 0 & detresults$Cancer==1)

#Calculate Distribution of Stages of Screen Detected Cancers
screen_stage_props<-data.frame("Names"=c("I","II","III","IV","DCIS"),
                               "Predicted"=c(tabulate(screendet$Stage)/sum(tabulate(screendet$Stage))),
                               "Observed"=c(0.4020,0.2310,0.1310,0.0240,0.2120))

#Caluclate Distribution of Stages of all cancers
all_stage_props<-data.frame("Names"=c("I","II","III","IV","DCIS"),
                            "Predicted"=c(tabulate(detresults$Stage)/sum(tabulate(detresults$Stage))),
                            "Observed"=c(0.3940,0.3530,0.0800,0.0440,0.1290))

#Calculate mean cancer sizes
cancer_sizes<-data.frame("Detection"=c("Screen Detected","Clinically Detected"),
                         "Size"=c(mean(screendet$`Cancer Size`),mean(clindet$`Cancer Size`)))

#Plot cancer detection sizes
plot(density(clindet$`Cancer Size`),col="black",main="Size of Cancers Detected",xlab="Diameter (mm)",ylim=c(0,0.08))
lines(density(screendet$`Cancer Size`),col="blue")
legend("right",legend=c("Clinically Detected","Screen Detected"),fill=c("black","blue"))


#Calculate calculate age band incidence rates
#Drop rows with no death age
n_missing <- sum(is.na(detresults$`Death Age`))
if (n_missing > 0) {
  message(sprintf("Dropping %d rows with missing Death Age", n_missing))
  detresults <- detresults[!is.na(detresults$`Death Age`), ]
}

#Define age bands
breaks <- c(seq(0, 90, by = 5), Inf)
labels <- c(
  sprintf("%02d to %02d", seq(0, 85, by = 5), seq(4, 89, by = 5)),
  "90+"
)

# Person-years and cases per band 
# A woman followed from age 0 until her Death Age contributes
# max(0, min(DeathAge, upper) - lower) person-years to band [lower, upper)
person_years_in_band <- function(death_age, lower, upper) {
  pmax(0, pmin(death_age, upper) - lower)
}

incidence <- data.frame(
  `Age Range` = labels,
  Cases = NA_real_,
  Person_Years = NA_real_,
  check.names = FALSE
)

for (i in seq_along(labels)) {
  lower <- breaks[i]
  upper <- breaks[i + 1]
  
  incidence$Person_Years[i] <- sum(person_years_in_band(detresults$`Death Age`, lower, upper))
  
  incidence$Cases[i] <- sum(
    detresults$Cancer == 1 &
      detresults$`Cancer Diagnosed Age` >= lower &
      detresults$`Cancer Diagnosed Age` < upper,
    na.rm = TRUE
  )
}

# Incidence rate per 100,000 person-years
incidence$Rate_per_100000 <- round(incidence$Cases / incidence$Person_Years * 100000, 1)
incidence$Observed<-c(0,0,0,0.1,1.6,11.5,31.2,65.8,124.6,214.8,279.8,285.5,337.9,
                      412.3,372.7,403,430.4,447.7,448.4)
plot(incidence$Rate_per_100000,type="l",ylim=c(0,450))
lines(incidence$Observed)

incidence_long <- pivot_longer(
  incidence,
  cols = c(Rate_per_100000, Observed),
  names_to = "Source",
  values_to = "Rate"
)

ggplot(incidence_long, aes(x = `Age Range`, y = Rate, color = Source, group = Source)) +
  geom_line() +
  geom_point() +
  labs(x = "Age Range", y = "Incidence per 100,000 person-years",
       title = "Simulated vs Observed Cancer Incidence by Age Band") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(incidence, row.names = FALSE)

