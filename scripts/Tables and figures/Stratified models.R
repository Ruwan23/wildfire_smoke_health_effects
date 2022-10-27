# "PHIRE Data Analysis - Stratified risk estimates"
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
data5 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")
data5 <- data.frame(data5)
# removing columns that won't be analyzed
data6 <- data5 %>%
          dplyr::select(-c("ed_allcause_cv_0_1","ed_allcause_cv_0_4","ed_allcause_cv_5_18","ed_allcause_cv_19_64",
                    "ed_copd_0_1","ed_copd_0_4","ed_copd_5_18","ed_copd_19_64",
                    "ed_asthma_45_64",
                    "ed_allcause_resp_45_64",
                    "pdd_allcause_cv_0_1","pdd_allcause_cv_0_4","pdd_allcause_cv_5_18","pdd_allcause_cv_19_64",
                    "pdd_copd_0_1","pdd_copd_0_4","pdd_copd_5_18","pdd_copd_19_64",
                    "pdd_asthma_45_64",
                    "pdd_allcause_resp_45_64"))



#------------------------------------------------
# STRATIFIED MODELS - EMERGENCY DEPARTMENT VISITS
#------------------------------------------------

# note that each outcome has its own lag length chosen based on best BIC from previous model building

# zips <- sample(unique(data5$PATZIP),20)
# test <- data6 %>% filter(PATZIP %in% zips)

#asthma -  lag 0-2
start <- which(names(data6)=="ed_asthma")
end <- which(names(data6)=="ed_asthma_other")
ed_asthma_estimates <- data.frame(matrix(NA,nrow=15,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=2, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_asthma_estimates[(i-start)+1,1] <- out
  ed_asthma_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_asthma_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_asthma_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(ed_asthma_estimates) <- c("Group","RR","LB","UB")
saveRDS(ed_asthma_estimates,"output/Risk estimate models/ed_asthma_estimates.rds")

#copd - lag 0-1
start <- which(names(data6)=="ed_copd")
end <- which(names(data6)=="ed_copd_other")
ed_copd_estimates <- data.frame(matrix(NA,nrow=13,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_copd_estimates[(i-start)+1,1] <- out
  ed_copd_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_copd_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_copd_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(ed_copd_estimates) <- c("Group","RR","LB","UB")
saveRDS(ed_copd_estimates,"output/Risk estimate models/ed_copd_estimates.rds")


# Respiratory - lag 0-4
start <- which(names(data6)=="ed_allcause_resp")
end <- which(names(data6)=="ed_allcause_resp_other")
ed_allcause_resp_estimates <- data.frame(matrix(NA,nrow=15,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=4, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  ed_allcause_resp_estimates[(i-start)+1,1] <- out
  ed_allcause_resp_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  ed_allcause_resp_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  ed_allcause_resp_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(ed_allcause_resp_estimates) <- c("Group","RR","LB","UB")
saveRDS(ed_allcause_resp_estimates,"output/Risk estimate models/ed_allcause_resp_estimates.rds")


# CV - lag0 
start <- which(names(data6)=="ed_allcause_cv")
end <- which(names(data6)=="ed_allcause_cv_other")
ed_allcause_cv_estimates <- data.frame(matrix(NA,nrow=13,ncol=4))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ pm25_residual_cz + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  ed_allcause_cv_estimates[(i-start)+1,1] <- out
  ed_allcause_cv_estimates[(i-start)+1,2] <- round((exp((unname(fixef(model))*10)[2])-1)*100,2)
  ed_allcause_cv_estimates[(i-start)+1,3] <- round((exp((unname(fixef(model))*10)[2] - 1.96*(unname(se.fixef(model))*10)[2])-1)*100,2)
  ed_allcause_cv_estimates[(i-start)+1,4] <- round((exp((unname(fixef(model))*10)[2] + 1.96*(unname(se.fixef(model))*10)[2])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(ed_allcause_cv_estimates) <- c("Group","RR","LB","UB")
saveRDS(ed_allcause_cv_estimates,"output/Risk estimate models/ed_allcause_cv_estimates.rds")

rm(list=ls())

# creating figure
ed_asthma_estimates <- readRDS("output/Risk estimate models/ed_asthma_estimates.rds")
ed_copd_estimates <- readRDS("output/Risk estimate models/ed_copd_estimates.rds")
ed_allcause_resp_estimates <- readRDS("output/Risk estimate models/ed_allcause_resp_estimates.rds")
ed_allcause_cv_estimates <- readRDS("output/Risk estimate models/ed_allcause_cv_estimates.rds")
results_ed <- rbind(ed_asthma_estimates,ed_copd_estimates,ed_allcause_resp_estimates,ed_allcause_cv_estimates)
results_ed$Outcome <- ifelse(grepl("asthma",results_ed$Group),"Asthma (Lag 0-2)",
                           ifelse(grepl("copd",results_ed$Group),"COPD (Lag 0-1)",
                                  ifelse(grepl("allcause_resp",results_ed$Group),"Respiratory* (Lag 0-4)","Cardiovascular** (Lag 0)"
                                  )))
results_ed$Subgroup <- ifelse(grepl("0_1",results_ed$Group),"Age 0-1",
ifelse(grepl("0_4",results_ed$Group),"Age 0-4",
ifelse(grepl("5_18",results_ed$Group),"Age 5-18",
ifelse(grepl("19_44",results_ed$Group),"Age 19-44",
ifelse(grepl("45_64",results_ed$Group),"Age 45-64",
ifelse(grepl("19_64",results_ed$Group),"Age 19-64",
ifelse(grepl("65_plus",results_ed$Group),"Age 65+",
ifelse(grepl("female",results_ed$Group),"Female",
ifelse(grepl("male",results_ed$Group),"Male",
ifelse(grepl("white",results_ed$Group),"NH White",
ifelse(grepl("black",results_ed$Group),"NH Black",
ifelse(grepl("hispanic",results_ed$Group),"Hispanic",
ifelse(grepl("asian",results_ed$Group),"NH Asian",
ifelse(grepl("nam",results_ed$Group),"NH Nat Am",
ifelse(grepl("other",results_ed$Group),"NH Other",
ifelse(grepl("public",results_ed$Group),"Public",
ifelse(grepl("private",results_ed$Group),"Private","Overall"
)))))))))))))))))

results_ed$Subgroup <- factor(results_ed$Subgroup,levels = c("Overall","Age 0-1","Age 0-4","Age 5-18","Age 19-44","Age 45-64","Age 19-64","Age 65+",
                                                       "Female","Male","NH White","NH Black","Hispanic","NH Asian",
                                                       "NH Nat Am","NH Other","Public","Private"))

saveRDS(results_ed,"output/Risk estimate models/ed_all_results.rds")


#------------------------------------------------
# STRATIFIED MODELS - HOSPITALIZATIONS
#------------------------------------------------

# note that each outcome has its own lag length chosen based on best BIC from previous model building

# zips <- sample(unique(data5$PATZIP),20)
# test <- data6 %>% filter(PATZIP %in% zips)


#asthma -  lag 0-1
start <- which(names(data6)=="pdd_asthma")
end <- which(names(data6)=="pdd_asthma_other")
pdd_asthma_estimates <- data.frame(matrix(NA,nrow=15,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_asthma_estimates[(i-start)+1,1] <- out
  pdd_asthma_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_asthma_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_asthma_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(pdd_asthma_estimates) <- c("Group","RR","LB","UB")
saveRDS(pdd_asthma_estimates,"output/Risk estimate models/pdd_asthma_estimates.rds")

#copd - lag 0-1
start <- which(names(data6)=="pdd_copd")
end <- which(names(data6)=="pdd_copd_other")
pdd_copd_estimates <- data.frame(matrix(NA,nrow=13,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_copd_estimates[(i-start)+1,1] <- out
  pdd_copd_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_copd_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_copd_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(pdd_copd_estimates) <- c("Group","RR","LB","UB")
saveRDS(pdd_copd_estimates,"output/Risk estimate models/pdd_copd_estimates.rds")


# Respiratory - lag 0-1
start <- which(names(data6)=="pdd_allcause_resp")
end <- which(names(data6)=="pdd_allcause_resp_other")
pdd_allcause_resp_estimates <- data.frame(matrix(NA,nrow=15,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_allcause_resp_estimates[(i-start)+1,1] <- out
  pdd_allcause_resp_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_allcause_resp_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_allcause_resp_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(pdd_allcause_resp_estimates) <- c("Group","RR","LB","UB")
saveRDS(pdd_allcause_resp_estimates,"output/Risk estimate models/pdd_allcause_resp_estimates.rds")


# CV - lag 0-1
start <- which(names(data6)=="pdd_allcause_cv")
end <- which(names(data6)=="pdd_allcause_cv_other")
pdd_allcause_cv_estimates <- data.frame(matrix(NA,nrow=13,ncol=4))
cb.pm <- crossbasis(data6$pm25_residual_cz, lag=1, argvar=list(fun="lin"),
                    arglag=list(fun="integer"))
for (i in start:end){
  out <- names(data6)[i]
  names(data6)[i] <- "a"
  model <- glmer(a ~ cb.pm + 
                   ns(hi,3) + dow + ns(dayseq,36) + 
                   (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
  saveRDS(model,paste("output/Risk estimate models/MODEL_",out,".rds",sep=""))
  names(data6)[i] <- out
  pred1.pm <- crosspred(cb.pm, model, cumul=TRUE)
  pdd_allcause_cv_estimates[(i-start)+1,1] <- out
  pdd_allcause_cv_estimates[(i-start)+1,2] <- round((unname(pred1.pm$allRRfit["10"])-1)*100,2)
  pdd_allcause_cv_estimates[(i-start)+1,3] <- round((unname(pred1.pm$allRRlow["10"])-1)*100,2)
  pdd_allcause_cv_estimates[(i-start)+1,4] <- round((unname(pred1.pm$allRRhigh["10"])-1)*100,2)
  print(paste("Completed",out,sep=" "))
}
names(pdd_allcause_cv_estimates) <- c("Group","RR","LB","UB")
saveRDS(pdd_allcause_cv_estimates,"output/Risk estimate models/pdd_allcause_cv_estimates.rds")


rm(list=ls())

# creating figure
pdd_asthma_estimates <- readRDS("output/Risk estimate models/pdd_asthma_estimates.rds")
pdd_copd_estimates <- readRDS("output/Risk estimate models/pdd_copd_estimates.rds")
pdd_allcause_resp_estimates <- readRDS("output/Risk estimate models/pdd_allcause_resp_estimates.rds")
pdd_allcause_cv_estimates <- readRDS("output/Risk estimate models/pdd_allcause_cv_estimates.rds")
results_pdd <- rbind(pdd_asthma_estimates,pdd_copd_estimates,pdd_allcause_resp_estimates,pdd_allcause_cv_estimates)
results_pdd$Outcome <- ifelse(grepl("asthma",results_pdd$Group),"Asthma (Lag 0-2)",
                             ifelse(grepl("copd",results_pdd$Group),"COPD (Lag 0-1)",
                                    ifelse(grepl("allcause_resp",results_pdd$Group),"Respiratory* (Lag 0-4)","Cardiovascular** (Lag 0)"
                                    )))
results_pdd$Subgroup <- ifelse(grepl("0_1",results_pdd$Group),"Age 0-1",
                              ifelse(grepl("0_4",results_pdd$Group),"Age 0-4",
                                     ifelse(grepl("5_18",results_pdd$Group),"Age 5-18",
                                            ifelse(grepl("19_44",results_pdd$Group),"Age 19-44",
                                                   ifelse(grepl("45_64",results_pdd$Group),"Age 45-64",
                                                          ifelse(grepl("19_64",results_pdd$Group),"Age 19-64",
                                                                 ifelse(grepl("65_plus",results_pdd$Group),"Age 65+",
                                                                        ifelse(grepl("female",results_pdd$Group),"Female",
                                                                               ifelse(grepl("male",results_pdd$Group),"Male",
                                                                                      ifelse(grepl("white",results_pdd$Group),"NH White",
                                                                                             ifelse(grepl("black",results_pdd$Group),"NH Black",
                                                                                                    ifelse(grepl("hispanic",results_pdd$Group),"Hispanic",
                                                                                                           ifelse(grepl("asian",results_pdd$Group),"NH Asian",
                                                                                                                  ifelse(grepl("nam",results_pdd$Group),"NH Nat Am",
                                                                                                                         ifelse(grepl("other",results_pdd$Group),"NH Other",
                                                                                                                                ifelse(grepl("public",results_pdd$Group),"Public",
                                                                                                                                       ifelse(grepl("private",results_pdd$Group),"Private","Overall"
                                                                                                                                       )))))))))))))))))

results_pdd$Subgroup <- factor(results_pdd$Subgroup,levels = c("Overall","Age 0-1","Age 0-4","Age 5-18","Age 19-44","Age 45-64","Age 19-64","Age 65+",
                                                             "Female","Male","NH White","NH Black","Hispanic","NH Asian",
                                                             "NH Nat Am","NH Other","Public","Private"))

saveRDS(results_pdd,"output/Risk estimate models/pdd_all_results.rds")
 


