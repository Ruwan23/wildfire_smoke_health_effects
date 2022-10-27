# RANDOM EFFECTS META ANALYSIS
#packages
memory.limit(size=100000)
library(dlnm)
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
library(biglm)
library(speedglm)
library(sp)
library(leaflet)
library(raster)
library(rgdal)
library(geosphere)
library(rgeos)
library(wesanderson)
library(stats)
library(devtools)
library(sf)
library(lme4)
library(glmmADMB)
library(weathermetrics)
library(lubridate)
library(effects)
library(sjPlot)
library(broom.mixed)
library(rmeta)

# Reading in dataset
data4 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")
data5 <- data4 %>%
  group_by(PATZIP) %>%
  mutate(pm25_residual_cz_lag1 = lag(pm25_residual_cz, n = 1), 
         pm25_residual_cz_lag2 = lag(pm25_residual_cz, n = 2),
         pm25_residual_cz_lag3 = lag(pm25_residual_cz, n = 3),
         pm25_residual_cz_lag4 = lag(pm25_residual_cz, n = 4),
         pm25_residual_cz_lag5 = lag(pm25_residual_cz, n = 5)
  )

# first 6 obs
head(data5)

# estimate RR and SE at each zip
results <- as.data.frame(matrix(NA,length(unique(data5$zip3)),3))

# run models
zip3_list <- unique(data5$zip3)
for (i in 1:length(zip3_list)){
  tic <- Sys.time()
  subset <- data5 %>% filter(zip3==zip3_list[i])
  cb02.pm <- crossbasis(subset$pm25_residual_cz, lag=2, argvar=list(fun="lin"),arglag=list(fun="integer"),group=subset$PATZIP)
  model <- glm(ed_asthma ~ cb02.pm + ns(hi,3) + dow + ns(dayseq,36) + factor(PATZIP),data = subset,family=poisson(link="log"))
  pred1.pm <- crosspred(cb02.pm, model, at=0:10, cumul=TRUE)
  log_rr <- log(unname(pred1.pm$allRRfit["10"]))
  se <- (log(unname(pred1.pm$allRRhigh["10"]))-log_rr)/1.96
  results[i,1] <- zip3_list[i]
  results[i,2] <- log_rr
  results[i,3] <- se
  toc <- Sys.time()
  print(i)
  print(toc-tic)
  
}
names(results) <- c("ZIP3","log_rr","se")

# calculate meta analytic effect
meta_re <- meta.summaries(d=results$log_rr, se=results$se, method="random")
est_re <- c(exp(meta_re$summary),exp(meta_re$summary-1.96*meta_re$se.summary),exp(meta_re$summary+1.96*meta_re$se.summary))
names(est_re) <- c("RR","LB","UB")
est_re


#BOOTSTRAP
bs <- read.csv("data/bootstrap.csv")
std <- function(x) sd(x)/sqrt(length(x))
se <- std(bs$log_rr)
pe <- mean(bs$log_rr)
est <- c(round(exp(pe),4),round(exp(pe-1.96*se),4),round(exp(pe+1.96*se),4))
names(est) <- c("RR", "95% LB", "95% UB")

 
