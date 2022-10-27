# Ruwan Thilakaratne
# Calculating Moran's I for models 4a, 6a, 9, 10

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
library(spdep)

# Reading in dataset
data4 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")

# Data header
head(data4)

#random sample of zips
# zip_count=10
# zips <- sample(unique(data4$PATZIP),zip_count)
# data4 <- data4 %>% filter(PATZIP %in% zips)

# SPATIAL AUTOCORRELATION
ca_zip3 <- seq(900,961,1)
ca_zip_map <- st_read("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/CA zips/tl_2019_us_zcta510/tl_2019_us_zcta510.shp")
ca_zip_map <- ca_zip_map %>% dplyr::filter(substr(GEOID10,1,3) %in% ca_zip3)
ca_zip_map <- ca_zip_map %>% rename("PATZIP" = "GEOID10")
ca_zip_map$PATZIP <- as.numeric(ca_zip_map$PATZIP)

mi_results <- data.frame()

#---------------------------------------------
# Model 4a

#prediction from model 4a
model_4a <- speedglm(ed_asthma ~ pm25 + ns(hi,3) + dow + zip3 + ns(dayseq,36) + offset(log(popsize)),data = data4,family=poisson(link="log"))
prediction <- predict(model_4a, newdata=data4)
data4$pred_m4a <- prediction
data4$resid_m4a <- data4$ed_asthma - data4$pred_m4a
zip_mean_residual <- data4 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_m4a))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- left_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_4a <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("Model 4a",unname(round(moran_4a$estimate[1],8)),unname(round(moran_4a$statistic,8))))

#---------------------------------------
#prediction from model 6a
model_6a <- speedglm(ed_asthma ~ pm25_residual_zip3 + ns(hi,3) + dow + ns(dayseq,36) + zip3 + offset(log(popsize)),data = data4,family=poisson(link="log"))
prediction <- predict(model_6a, newdata=data4)
data4$pred_m6a <- prediction
data4$resid_m6a <- data4$ed_asthma - data4$pred_m6a
zip_mean_residual <- data4 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_m6a))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- left_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_6a <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("Model 6a",unname(round(moran_6a$estimate[1],8)),unname(round(moran_6a$statistic,8))))

#---------------------------------------
#prediction from model 9
model_9 <- speedglm(ed_asthma ~ pm25_residual_cz + ns(hi,3) + dow + zip3 + ns(dayseq,36) + offset(log(popsize)),data = data4,family=poisson(link="log"))
prediction <- predict(model_9, newdata=data4)
data4$pred_m9 <- prediction
data4$resid_m9 <- data4$ed_asthma - data4$pred_m9
zip_mean_residual <- data4 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_m9))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- left_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_9 <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("Model 9",unname(round(moran_9$estimate[1],8)),unname(round(moran_9$statistic,8))))


#---------------------------------------
#prediction from model 10
model_10 <- speedglm(ed_asthma ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + hpi + offset(log(popsize)),data = data4,family=poisson(link="log"))
prediction <- predict(model_10, newdata=data4)
data4$pred_m10 <- prediction
data4$resid_m10 <- data4$ed_asthma - data4$pred_m10
zip_mean_residual <- data4 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_m10))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- left_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_10 <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("Model 10","PM2.5 CZ residuals + HPI(linear)",unname(round(moran_10$estimate[1],8))))

#labeling results
names(mi_results) <- c("Model","Moran's I","Moran's I SDs")

#----------------------------------------------------------------------
# MORAN'S I FOR NEW ZIP CODE RANDOM INTERCEPT MODEL + COMPARISON TO OLD
#----------------------------------------------------------------------

mi_results <- data.frame()

#---------------------------------------
#prediction from new model - zip code random intercept
model_dlm_10 <- glmer(ed_asthma ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
prediction <- predict(model_dlm_10, newdata=data4)
data4$pred_dlm10 <- prediction
data4$resid_dlm10 <- data4$ed_asthma - data4$pred_dlm10
zip_mean_residual <- data4 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_dlm10))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- left_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_10 <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("DLM ","PM2.5 CZ residuals - lag0, zip RI",unname(round(moran_10$estimate[1],8))))


#---------------------------------------
#prediction from new model, lag0-5
cb05.pm <- crossbasis(data4$pm25_residual_cz, lag=5, argvar=list(fun="lin"),arglag=list(fun="integer"),group=data4$PATZIP)
model_dlm_10 <- glmer(ed_asthma ~ cb05.pm + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
prediction <- predict(model_dlm_10, newdata=data4)
data4$pred_dlm10 <- prediction
data4$resid_dlm10 <- data4$ed_asthma - data4$pred_dlm10
data5 <- data4[complete.cases(data4),]
zip_mean_residual <- data5 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_dlm10))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- inner_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_10 <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("DLM ","PM2.5 CZ residuals - lag0-5, zip RI",unname(round(moran_10$estimate[1],8))))

#---------------------------------------
#prediction from old model, lag0-5
cb05.pm <- crossbasis(data4$pm25_residual_cz, lag=5, argvar=list(fun="lin"),arglag=list(fun="integer"),group=data4$PATZIP)
model_dlm_10 <- speedglm(ed_asthma ~ cb05.pm + ns(hi,3) + dow + ns(dayseq,36) + zip3 + offset(log(popsize)),data = data4,family=poisson(link="log"))
prediction <- predict(model_dlm_10, newdata=data4)
data4$pred_dlm10 <- prediction
data4$resid_dlm10 <- data4$ed_asthma - data4$pred_dlm10
data5 <- data4[complete.cases(data4),]
zip_mean_residual <- data5 %>% group_by(PATZIP) %>% summarize(mean_residual = mean(resid_dlm10))
zip_mean_residual <- data.frame(zip_mean_residual)

#merging with polygons
ca_zip_map_resid <- inner_join(ca_zip_map,zip_mean_residual,by=c("PATZIP"))
ca_zip_map_resid <- ca_zip_map_resid %>% filter(PATZIP %in% unique(data4$PATZIP))

#creating neighborhood matrix
resid_nb <- poly2nb(ca_zip_map_resid)

# set weights for matrix
resid_w <- nb2listw(resid_nb,zero.policy=TRUE)

# calculate Moran's I
moran_10 <- moran.test(ca_zip_map_resid$mean_residual,listw=resid_w,zero.policy = TRUE)
mi_results <- rbind(mi_results,c("DLM ","PM2.5 CZ residuals - lag0-5, zip offset (old model)",unname(round(moran_10$estimate[1],8))))

#labeling results
names(mi_results) <- c("Model","Moran's I","Moran's I SDs")
