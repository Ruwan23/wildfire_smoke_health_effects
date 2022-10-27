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

# Reading in 2008 and 2009 NOAA and Harvard data
ca_zip3 <- as.character(seq(900,961,1))

#reading in additional NOAA data from 2008 - provided as daily for United States
list <- dir(path="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/NOAA_ZIP_2008/",pattern = "*.rds")
noaa_2008 <- data.frame()
for (i in 1:length(list)){
  tic <- Sys.time()
  data <- readRDS(paste("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/NOAA_ZIP_2008/",list[i],sep=""))
  data <- data %>% filter(substr(ZIP,1,3) %in% ca_zip3)
  date <- as.Date(substr(list[i],10,19))
  data$serv_dt <- rep(date,nrow(data))
  noaa_2008 <- rbind(noaa_2008,data)
  print(i)
  toc <- Sys.time()
  print(toc-tic)
}
noaa_2008 <- noaa_2008 %>% rename(PATZIP = ZIP,
                                  POP.Wght.noaa = POP.Wght,
                                  N.tracts.noaa = N.tracts)
noaa_2008$PATZIP <- as.numeric(noaa_2008$PATZIP)

#reading in additional NOAA data from 2009
list <- dir(path="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/NOAA_ZIP_2009/",pattern = "*.rds")
lista <- list[seq(1,142)]
listb <- list[seq(144,length(list))]
list <- c(lista,listb)
noaa_2009 <- data.frame()
for (i in 1:length(list)){
  tic <- Sys.time()
  data <- readRDS(paste("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/NOAA_ZIP_2009/",list[i],sep=""))
  data <- data %>% filter(substr(ZIP,1,3) %in% ca_zip3)
  date <- as.Date(substr(list[i],10,19))
  data$serv_dt <- rep(date,nrow(data))
  noaa_2009 <- rbind(noaa_2009,data)
  print(i)
  toc <- Sys.time()
  print(toc-tic)
}
noaa_2009 <- noaa_2009 %>% rename(PATZIP = ZIP,
                                  POP.Wght.noaa = POP.Wght,
                                  N.tracts.noaa = N.tracts)
noaa_2009$PATZIP <- as.numeric(noaa_2009$PATZIP)

#saving data
noaa_2008_2009 <- rbind(noaa_2008,noaa_2009)
saveRDS(noaa_2008_2009,file="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/ca_noaa_2008_2009.rds")

#reading in additional Harvard data from 2008
list <- dir(path="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/Harvard_ZIP_2008/",pattern = "*.rds")
harvard_2008 <- data.frame()
for (i in 1:length(list)){
  tic <- Sys.time()
  data <- readRDS(paste("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/Harvard_ZIP_2008/",list[i],sep=""))
  data <- data %>% filter(substr(ZIP,1,3) %in% ca_zip3)
  harvard_2008 <- rbind(harvard_2008,data)
  print(i)
  toc <- Sys.time()
  print(toc-tic)
}
harvard_2008 <- harvard_2008 %>% rename(PATZIP = ZIP,
                                        POP.Wght.harv = POP.Wght,
                                        N.tracts.harv = N.tracts,
                                        serv_dt = Date)
harvard_2008$PATZIP <- as.numeric(harvard_2008$PATZIP)

#reading in additional Harvard data from 2009
list <- dir(path="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/Harvard_ZIP_2009/",pattern = "*.rds")
harvard_2009 <- data.frame()
for (i in 1:length(list)){
  tic <- Sys.time()
  data <- readRDS(paste("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/Harvard_ZIP_2009/",list[i],sep=""))
  data <- data %>% filter(substr(ZIP,1,3) %in% ca_zip3)
  harvard_2009 <- rbind(harvard_2009,data)
  print(i)
  toc <- Sys.time()
  print(toc-tic)
}
harvard_2009 <- harvard_2009 %>% rename(PATZIP = ZIP,
                                        POP.Wght.harv = POP.Wght,
                                        N.tracts.harv = N.tracts,
                                        serv_dt = Date)
harvard_2009$PATZIP <- as.numeric(harvard_2009$PATZIP)

#saving data
harvard_2008_2009 <- rbind(harvard_2008,harvard_2009)
saveRDS(harvard_2008_2009,file="S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data For Analysis/ca_harvard_2008_2009.rds")








