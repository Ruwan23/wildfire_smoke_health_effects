# Improving run time for random intercept models
# Ruwan Thilakaratne

#packages
memory.limit(size=100000)
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

data3 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")

# testing impact of varying nAGQ to speed up fitting time
# nAGQ=1 is default Laplace approx
# nAGQ=0 is PIRLS estimation (faster)

# random subset of 100 zip codes
sub_size <- 10
zips <- unique(data3$PATZIP)
samp <- sample(zips,size=sub_size)
data_sub <- data3 %>% filter(PATZIP %in% samp)

# nAGQ = 0 (faster)
tic <- Sys.time()
model_sub_ri_nacq0 <- glmer(ed_asthma ~  pm25 + (1|zip3),data = data_sub,family=poisson(link="log"),nAGQ=0)
toc <- Sys.time()
toc-tic
summary(model_sub_ri_nacq0)

# nAGQ = 1 (default, slower)
tic <- Sys.time()
model_sub_ri_nacq1 <- glmer(ed_asthma ~  pm25 + (1|zip3),data = data_sub,family=poisson(link="log"),nAGQ=1)
toc <- Sys.time()
toc-tic
summary(model_sub_ri_nacq1)

test1 <- data.frame(nagq0_coef = unname(summary(model_sub_ri_nacq0)$coefficients[,1]),
                    nagq1_coef = unname(summary(model_sub_ri_nacq1)$coefficient[,1]),
                    nagq0_se = unname(summary(model_sub_ri_nacq0)$coefficients[,2]),
                    nagq1_se = unname(summary(model_sub_ri_nacq1)$coefficient[,2]),
                    coef = (((unname(summary(model_sub_ri_nacq0)$coefficients[,1])) - (unname(summary(model_sub_ri_nacq1)$coefficient[,1])))/(unname(summary(model_sub_ri_nacq1)$coefficients[,1])))*100,
                    se = (((unname(summary(model_sub_ri_nacq0)$coefficients[,2])) - (unname(summary(model_sub_ri_nacq1)$coefficient[,2])))/(unname(summary(model_sub_ri_nacq1)$coefficients[,2])))*100)
test1 <- test1[2,] # removing grand intercept comparison and leaving pM2.5 coef comparison
names(test1) <- c("nAGQ 0 PM2.5 coef","nAGQ 1 PM2.5 coef","nAGQ 0 PM2.5 se","nAGQ 1 PM2.5 se","Percent change in coef","Percent change in SE")
write.csv(test1,"output/Percent difference in coefficient and standard error between nAGQ=1 vs. nAGQ=0.csv")

















