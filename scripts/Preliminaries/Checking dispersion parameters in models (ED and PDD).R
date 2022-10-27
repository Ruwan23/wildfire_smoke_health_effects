# Figure - HA risk estimates

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
library(sf)
library(lme4)
library(weathermetrics)
library(lubridate)
library(effects)
library(sjPlot)
library(arm)
library(broom.mixed)
library(dlnm)
library(flextable)


# checking for overdispersion (disp values = 1 if properly dispersed; <1 okay)

# ED
# asthma
model_ed_asthma <- readRDS("output/Risk estimate models/MODEL_ed_asthma.rds")
disp_ed_asthma <- deviance(model_ed_asthma)/df.residual(model_ed_asthma)
disp_ed_asthma
rm(model_ed_asthma)
# COPD
model_ed_copd <- readRDS("output/Risk estimate models/MODEL_ed_copd.rds")
disp_ed_copd <- deviance(model_ed_copd)/df.residual(model_ed_copd)
disp_ed_copd
rm(model_ed_copd)
# all-resp
model_ed_allcause_resp <- readRDS("output/Risk estimate models/MODEL_ed_allcause_resp.rds")
disp_ed_allcause_resp <- deviance(model_ed_allcause_resp)/df.residual(model_ed_allcause_resp)
disp_ed_allcause_resp
rm(model_ed_allcause_resp)
# all-CV
model_ed_allcause_cv <- readRDS("output/Risk estimate models/MODEL_ed_allcause_cv.rds")
disp_ed_allcause_cv <- deviance(model_ed_allcause_cv)/df.residual(model_ed_allcause_cv)
disp_ed_allcause_cv
rm(model_ed_allcause_cv)

# PDD
# asthma
model_pdd_asthma <- readRDS("output/Risk estimate models/MODEL_pdd_asthma.rds")
disp_pdd_asthma <- deviance(model_pdd_asthma)/df.residual(model_pdd_asthma)
disp_pdd_asthma
rm(model_pdd_asthma)
# COPD
model_pdd_copd <- readRDS("output/Risk estimate models/MODEL_pdd_copd.rds")
disp_pdd_copd <- deviance(model_pdd_copd)/df.residual(model_pdd_copd)
disp_pdd_copd
rm(model_pdd_copd)
# all-resp
model_pdd_allcause_resp <- readRDS("output/Risk estimate models/MODEL_pdd_allcause_resp.rds")
disp_pdd_allcause_resp <- deviance(model_pdd_allcause_resp)/df.residual(model_pdd_allcause_resp)
disp_pdd_allcause_resp
rm(model_pdd_allcause_resp)
# all-CV
model_pdd_allcause_cv <- readRDS("output/Risk estimate models/MODEL_pdd_allcause_cv.rds")
disp_pdd_allcause_cv <- deviance(model_pdd_allcause_cv)/df.residual(model_pdd_allcause_cv)
disp_pdd_allcause_cv
rm(model_pdd_allcause_cv)

