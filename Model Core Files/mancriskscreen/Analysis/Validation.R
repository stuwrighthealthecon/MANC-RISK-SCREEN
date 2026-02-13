filenames <- list.files(det_output_path, full.names = TRUE)
alldata <- lapply(filenames, function(x) {
  get(load(x, .GlobalEnv))
})
alldata <- do.call("rbind", alldata)

#Remove women with cancer diagnosed before age 50
alldata <- alldata %>% filter(alldata[, 4] > 50 | alldata[, 4] == 0)
detresults <- alldata

screendet <- detresults %>% filter(detresults$`screen detected` == 1)
clindet <- detresults %>% filter(detresults$`screen detected` == 0)

tabulate(screendet$Stage)
tabulate(detresults$Stage)

mean(clindet$`Cancer Size`)
quantile(clindet$`Cancer Size`, probs = c(0.25, 0.75))
plot(density(clindet$`Cancer Size`))
mean(screendet$`Cancer Size`)
quantile(screendet$`Cancer Size`, probs = c(0.25, 0.75))
plot(density(screendet$`Cancer Size`))


cancersbyscreen <- tabulate(screendet$`Cancer Screen Number`)
median(cancersbyscreen)
