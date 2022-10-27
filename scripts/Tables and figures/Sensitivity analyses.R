# "PHIRE Data Analysis - SENSITIVITY ANALYSES
# Ruwan Thilakaratne
# 6/7/21

#packages
memory.limit(size=100000)
library(lattice)
library(gridExtra)
library(grid)
library(haven)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(tableone)
library(splines)
library(tidycensus)
library(tidyverse)
library(lubridate)
library(AER)
library(knitr)
library(splines)
library(data.table)
library(speedglm)
library(sp)
library(leaflet)
library(raster)
library(rgdal)
library(geosphere)
library(rgeos)
library(wesanderson)
library(stats)
library(sf)
library(lme4)
library(weathermetrics)
library(lubridate)
library(effects)
library(sjPlot)
library(arm)
library(broom.mixed)
library(dlnm)

unname <- function(x) {
  names(x) <- NULL  
  x
}

## Reading in data and creating covariates
data6 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")
data6 <- data.frame(data6)
#-----------------------------------------------
# STRATIFIED MODELS - EMERGENCY DEPARTMENT VISITS
#------------------------------------------------

# note that each outcome has its own lag length chosen based on best BIC from previous model building

# zips <- sample(unique(data5$PATZIP),20)
# test <- data6 %>% filter(PATZIP %in% zips)

dof <- c(36,45,54)

#asthma -  lag 0-2
ed_asthma_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=2, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(ed_asthma ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_ed_asthma_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_asthma_sens[i,1] <- "Asthma ED"
  ed_asthma_sens[i,2] <- dof[i]/9
  ed_asthma_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_asthma_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_asthma_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(ed_asthma_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(ed_asthma_sens,"output/Risk estimate models/Sens analysis - spline dof/ed_asthma_sens.rds")

#copd - lag 0-1
ed_copd_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(ed_copd ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_ed_copd_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_copd_sens[i,1] <- "copd ED"
  ed_copd_sens[i,2] <- dof[i]/9
  ed_copd_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_copd_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_copd_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(ed_copd_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(ed_copd_sens,"output/Risk estimate models/Sens analysis - spline dof/ed_copd_sens.rds")

#Resp - lag 0-4
ed_allcause_resp_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=4, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(ed_allcause_resp ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_ed_allcause_resp_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_allcause_resp_sens[i,1] <- "Respiratory ED"
  ed_allcause_resp_sens[i,2] <- dof[i]/9
  ed_allcause_resp_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_allcause_resp_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_allcause_resp_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(ed_allcause_resp_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(ed_allcause_resp_sens,"output/Risk estimate models/Sens analysis - spline dof/ed_allcause_resp_sens.rds")

# CV - lag0 
ed_allcause_cv_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
for (i in 1:length(dof)){
  model <- glmer(ed_allcause_cv ~ pm25_residual_cz + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_ed_allcause_cv_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  ed_allcause_cv_sens[i,1] <- "Cardiovascular ED"
  ed_allcause_cv_sens[i,2] <- dof[i]/9
  ed_allcause_cv_sens[i,3] <- round((exp((unname(fixef(model))*10)[2])-1)*100,2)
  ed_allcause_cv_sens[i,4] <- round((exp((unname(fixef(model))*10)[2] - 1.96*(unname(se.fixef(model))*10)[2])-1)*100,2)
  ed_allcause_cv_sens[i,5] <- round((exp((unname(fixef(model))*10)[2] + 1.96*(unname(se.fixef(model))*10)[2])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(ed_allcause_cv_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(ed_allcause_cv_sens,"output/Risk estimate models/Sens analysis - spline dof/ed_allcause_cv_sens.rds")

ed_allcause_cv_sens <- readRDS()

#######################
# HOSPITAL ADMISSIONS #
#######################
#asthma -  lag 0-1
pdd_asthma_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(pdd_asthma ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_pdd_asthma_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_asthma_sens[i,1] <- "Asthma HA"
  pdd_asthma_sens[i,2] <- dof[i]/9
  pdd_asthma_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_asthma_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_asthma_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(pdd_asthma_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(pdd_asthma_sens,"output/Risk estimate models/Sens analysis - spline dof/pdd_asthma_sens.rds")

#copd - lag 0-1
pdd_copd_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(pdd_copd ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_pdd_copd_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_copd_sens[i,1] <- "COPD HA"
  pdd_copd_sens[i,2] <- dof[i]/9
  pdd_copd_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_copd_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_copd_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(pdd_copd_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(pdd_copd_sens,"output/Risk estimate models/Sens analysis - spline dof/pdd_copd_sens.rds")

#Resp - lag 0-1
pdd_allcause_resp_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(pdd_allcause_resp ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_pdd_allcause_resp_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_allcause_resp_sens[i,1] <- "Respiratory HA"
  pdd_allcause_resp_sens[i,2] <- dof[i]/9
  pdd_allcause_resp_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_allcause_resp_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_allcause_resp_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(pdd_allcause_resp_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(pdd_allcause_resp_sens,"output/Risk estimate models/Sens analysis - spline dof/pdd_allcause_resp_sens.rds")

#CV - lag 0-1
pdd_allcause_cv_sens <- data.frame(matrix(NA,nrow=3,ncol=5))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in 1:length(dof)){
  model <- glmer(pdd_allcause_cv ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,dof[i]) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/Sens analysis - spline dof/MODEL_pdd_allcause_cv_spline_dof_",dof[i]/9,"_per_year.rds",sep=""))
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_allcause_cv_sens[i,1] <- "CV HA"
  pdd_allcause_cv_sens[i,2] <- dof[i]/9
  pdd_allcause_cv_sens[i,3] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_allcause_cv_sens[i,4] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_allcause_cv_sens[i,5] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",i,sep=" "))
}
names(pdd_allcause_cv_sens) <- c("Outcome","Time spline df/year","RR","LB","UB")
saveRDS(pdd_allcause_cv_sens,"output/Risk estimate models/Sens analysis - spline dof/pdd_allcause_cv_sens.rds")



# COMBINING AND SAVING
ed_asthma <- readRDS("output/Risk estimate models/Sens analysis - spline dof/ed_asthma_sens.rds")
ed_copd <- readRDS("output/Risk estimate models/Sens analysis - spline dof/ed_copd_sens.rds")
ed_allcause_cv <- readRDS("output/Risk estimate models/Sens analysis - spline dof/ed_allcause_cv_sens.rds")
ed_allcause_resp <- readRDS("output/Risk estimate models/Sens analysis - spline dof/ed_allcause_resp_sens.rds")
pdd_asthma <- readRDS("output/Risk estimate models/Sens analysis - spline dof/pdd_asthma_sens.rds")
pdd_copd <- readRDS("output/Risk estimate models/Sens analysis - spline dof/pdd_copd_sens.rds")
pdd_allcause_cv <- readRDS("output/Risk estimate models/Sens analysis - spline dof/pdd_allcause_cv_sens.rds")
pdd_allcause_resp <- readRDS("output/Risk estimate models/Sens analysis - spline dof/pdd_allcause_resp_sens.rds")
all <- rbind(ed_asthma,
             ed_copd,
             ed_allcause_cv,
             ed_allcause_resp,
             pdd_asthma,
             pdd_copd,
             pdd_allcause_cv,
             pdd_allcause_resp)
write.csv(all,"output/Risk estimate models/Sens analysis - spline dof/all_sensitivity_estimates.csv")
