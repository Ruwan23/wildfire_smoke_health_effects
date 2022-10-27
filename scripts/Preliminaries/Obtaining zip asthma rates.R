# Model component extraction for health burden assessment
# Ruwan Thilakaratne

#packages
memory.limit(size=100000)
library(dlnm)
library(forecast)
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
library(readxl)
library(weathermetrics)



#-------------------------------------
# BUILDING MODEL 
#-------------------------------------

# Reading in dataset
data4 <- readRDS("S:/Wildfire-Hoshiko/Wildfire - RThilakaratne/PHIRE_Health_Analysis/data/Asthma ED data with covariates.rds")

# Data header
head(data4)

# zips <- sample(unique(data4$PATZIP),10)
# df3 <- data4 %>% filter(PATZIP %in% zips)

data5 <- data4 %>%
  group_by(PATZIP) %>%
  mutate(pm25_residual_cz_lag1 = lag(pm25_residual_cz, n = 1), 
         pm25_residual_cz_lag2 = lag(pm25_residual_cz, n = 2),
         pm25_residual_cz_lag3 = lag(pm25_residual_cz, n = 3),
         pm25_residual_cz_lag4 = lag(pm25_residual_cz, n = 4),
         pm25_residual_cz_lag5 = lag(pm25_residual_cz, n = 5)
  )

## ASTHMA MODEL
# NOTE! 5-day DLM was used for comparing observed vs. expected rates. 2-day DLM was chosen as final model because of minimization of BIC

model_dlm_02 <- glmer(asthma_count ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + 
                        ns(hi,3) + factor(dow) + ns(dayseq,36) + (1|PATZIP),data = data5,family=poisson(link="log"),nAGQ = 0)
cb03.pm <- crossbasis(data4$pm25_residual_cz, lag=2, argvar=list(fun="lin"),arglag=list(df=2),group=data4$PATZIP)
model_dlm_02_nonlinear <- glmer(asthma_count ~ cb03.pm + 
                        ns(hi,3) + factor(dow) + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
pred1.pm <- crosspred(cb03.pm, model_dlm_02_nonlinear, at=0:10, bylag=0.2, cumul=TRUE)

#-------------------------------------
# CREATING PREDICTIONS FOR VALIDATION WITH RAW RATES
#-------------------------------------

# Zip-level rates for comparison (whole study period, including PM values, predict all and average by zip code)
data5$pred <- predict(model_dlm_02,newdata=data5)
data6 <- data5 %>% filter(!is.na(pred)) %>% mutate(exp_pred = exp(pred))
zip_avg <- data6 %>% 
            group_by(PATZIP) %>%
            summarize(mean_exp_pred = mean(exp_pred))
zip_avg$annual_rate <- (zip_avg$mean_exp_pred)*365
zip_popsize <- data6 %>% group_by(PATZIP) %>% summarize(popsize = mean(popsize))
zip_avg_merge <- inner_join(zip_avg,zip_popsize,by=c("PATZIP"))
write.csv(zip_avg_merge,"S:/Wildfire-Hoshiko/Wildfire-BButler1/RT data/predicted annual asthma visit rates.csv")

#merging in 2010 census zip code population size data
pop <- data.frame(get_decennial(geography = "zcta",
                                variables = "P001001",
                                year = 2010))
pop <- pop %>% rename(PATZIP=GEOID,
                      popsize=value)
pop$popsize <- as.numeric(pop$popsize)
pop$PATZIP <- as.numeric(pop$PATZIP)
pop <- pop %>% dplyr::select(PATZIP,popsize)
ca_zip3 <- seq(900,961,1)
pop_ca <- pop %>% filter(substr(PATZIP,1,3) %in% ca_zip3)
write.csv(pop_ca,"S:/Wildfire-Hoshiko/Wildfire-BButler1/RT data/2010 census population CA zips.csv")

#-------------------------------------
# CONSTRUCTING DATASET FOR BASELINE RATES (2016)
#-------------------------------------

#census_api_key("e04b85d7a670f07ea35f483ab3f48275136e69ce",install=TRUE)
# ^^^ to understand why I wrote this, see here: https://walker-data.com/tidycensus/articles/basic-usage.html

#reading in daily asthma ED visit data
data <- fread("data/Zip-level daily asthma ED visits, 2008-2016.csv")

#merging in 2010 census zip code population size data
pop <- data.frame(get_decennial(geography = "zcta",
                                variables = "P001001",
                                year = 2010))
pop <- pop %>% rename(PATZIP=GEOID,
                      popsize=value)
pop$popsize <- as.numeric(pop$popsize)
pop$PATZIP <- as.numeric(pop$PATZIP)
pop <- pop %>% dplyr::select(PATZIP,popsize)
ca_zip3 <- seq(900,961,1)
pop_ca <- pop %>% filter(substr(PATZIP,1,3) %in% ca_zip3)

#merging in census data
data_census <- left_join(data,pop,by=c("PATZIP"))

#reading in meteorological data
noaa <- readRDS("S:/Wildfire-Hoshiko/Wildfire - RThilakaratne/PHIRE_Health_Analysis/data/Data For Analysis/CA_NOAA.rds")
head(noaa)
noaa <- noaa %>% rename(serv_dt = Date,
                        PATZIP = ZIP,
                        POP.Wght.noaa = POP.Wght,
                        N.tracts.noaa = N.tracts)
noaa$PATZIP <- as.numeric(noaa$PATZIP)
noaa_2008_2009 <- readRDS("S:/Wildfire-Hoshiko/Wildfire - RThilakaratne/PHIRE_Health_Analysis/data/Data For Analysis/ca_noaa_2008_2009.rds")
noaa <- rbind(noaa,noaa_2008_2009)
data_census$serv_dt <- as.Date(data_census$serv_dt)
data_census_noaa <- inner_join(data_census,noaa,by=c("PATZIP","serv_dt"))

# Creating heat index and DOW
data_census_noaa$rel_hum <- ifelse(data_census_noaa$RH.Mean>100,100,data_census_noaa$RH.Mean)
data_census_noaa_hi <- data_census_noaa %>%
  mutate(dow = lubridate::wday(serv_dt, label = TRUE, abbr = TRUE),
         hi = heat.index(t = TEMP.Mean, rh = rel_hum,temperature.metric="celsius"))

# Adding date indicator for spline
serv_dt <- seq.Date(as.Date("2008/1/1"),as.Date("2016/12/31"),"days")
dayseq <- seq(1,length(serv_dt))
dayseq_for_merge <- data.frame(serv_dt,dayseq)
df <- inner_join(data_census_noaa_hi,dayseq_for_merge,by=c("serv_dt"))

# Adding empty PM columns
df_2016 <- as.data.frame(df %>% 
          filter(serv_dt>=as.Date("2016/01/01")) %>%
            mutate(pm25_residual_cz = 0,
                    pm25_residual_cz_lag1 = 0,
                    pm25_residual_cz_lag2 = 0))

#--------------------------------
# 2016 BASELINE RATE ESTIMATES
#--------------------------------

# Zip-level rates for health burden formula (2016)
# creating prediction dataset - limiting to 2016 and setting all PM predictors to 0
df_2016 <- data5 %>% 
            filter(serv_dt>=as.Date("2016/01/01")) %>%
            mutate(pm25_residual_cz = 0,
                   pm25_residual_cz_lag1 = 0,
                   pm25_residual_cz_lag2 = 0)
df_2016$pred <- predict(model_dlm_02,newdata=df_2016)
df_2016$exp_pred <- exp(df_2016$pred)
zip_avg_2016 <- df_2016 %>% 
  group_by(PATZIP) %>%
  summarize(mean_exp_pred = mean(exp_pred),
            popsize = mean(popsize))
zip_to_zcta <- read.csv("data/zip_to_zcta_2018.csv") %>%
                rename(PATZIP=ZIP_CODE)
zip_avg_2016_merged <- inner_join(zip_to_zcta,zip_avg_2016,by=c("PATZIP"))
zip_avg_2016_col <- zip_avg_2016_merged %>% group_by(ZCTA) %>% summarize(mean_exp_pred = mean(mean_exp_pred),
                                                                         popsize = mean(popsize)) # taking mean but really all values are the same within a ZCTA
rate_sum <- sum(zip_avg_2016_col$mean_exp_pred)
pop_sum <- sum(zip_avg_2016_col$popsize)
annual_rate <- (rate_sum/pop_sum)*365
beta <- sum(unname(fixef(model_dlm_02))[2:4])
outcome <- "Asthma ED"
demo_subgroup <- "All"
health_effects_comp <- data.frame(outcome = outcome,
                                  demo_subgroup = demo_subgroup,
                                  annual_rate = annual_rate,
                                  beta = beta)
write.csv(health_effects_comp,"S:/Wildfire-Hoshiko/Wildfire-BButler1/RT data/incidence rates 2016.csv")
write.csv(health_effects_comp,"output/Betas and incidence rates for health burden analysis/incidence rates 2016.csv")

# seeing where zcta and patzips differ
names(pop_ca) <- c("ZCTA","popsize")
merged <- inner_join(pop_ca, zip_to_zcta, by=c("ZCTA"))
merged$indicator <- 1
analysis_zips <- unique(data5$PATZIP)
merged2 <- merged %>% filter(PATZIP %in% analysis_zips)
merged_col <- merged2 %>% group_by(ZCTA) %>% summarize(avg_popsize = mean(popsize),
                                                      zips = sum(indicator))

# time detrend in outcome only
cb05.pm <- crossbasis(data5$pm25, lag=5, argvar=list(fun="lin"),arglag=list(fun="integer"),group=data5$PATZIP)
model_dlm_02_3 <- glmer(asthma_count ~ cb05.pm +
                          ns(hi,3) + factor(dow) + ns(dayseq,72) + (1|PATZIP),data = data5,family=poisson(link="log"),nAGQ = 0)
pred1.pm <- crosspred(cb05.pm, model_dlm_02_3, at=0:10, cumul=TRUE)
round(unname(pred1.pm$allRRfit["10"]),5)
round(unname(pred1.pm$allRRlow["10"]),5)
round(unname(pred1.pm$allRRhigh["10"]),5)














