#Ruwan Thilakaratne
#Concatenating ED datasets

#packages
memory.limit(size=100000) # increasing RAM usage for reading in large datasets
library(haven)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(lubridate)
library(tableone)

full_data <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")
full_data <- full_data %>% filter(year(serv_dt) %in% c(2014,2015,2016))
#totals
full_data$PATZIP3 <- substr(full_data$PATZIP,1,3)
full_data_total <- full_data %>%
  group_by(serv_dt) %>%
  summarize(ed_allcause_cv=sum(ed_allcause_cv),
            ed_allcause_resp = sum(ed_allcause_resp),
            ed_copd = sum(ed_copd),
            ed_asthma = sum(ed_asthma))

dates_vline <- as.Date(c("2015-10-01")) 
dates_vline <- which(full_data_total$serv_dt %in% dates_vline)
ggplot(data=full_data_total, aes(x=serv_dt,y=ed_allcause_cv)) +
  geom_point() +
  labs(x="Date",y="All-cause CV ED visits (count)") +
  ggtitle("CV ED time series") + 
  geom_vline(xintercept=as.numeric(full_data_total$serv_dt[dates_vline]),
             linetype=4, colour="black")
ggplot(data=full_data_total, aes(x=serv_dt,y=ed_allcause_resp)) +
  geom_point() +
  labs(x="Date",y="All-cause Resp ED visits (count)") +
  ggtitle("Resp ED time series") + 
  geom_vline(xintercept=as.numeric(full_data_total$serv_dt[dates_vline]),
             linetype=4, colour="black")
ggplot(data=full_data_total, aes(x=serv_dt,y=ed_copd)) +
  geom_point() +
  labs(x="Date",y="All-cause COPD ED visits (count)") +
  ggtitle("COPD ED time series") + 
  geom_vline(xintercept=as.numeric(full_data_total$serv_dt[dates_vline]),
             linetype=4, colour="black")
ggplot(data=full_data_total, aes(x=serv_dt,y=ed_asthma)) +
  geom_point() +
  labs(x="Date",y="All-cause asthma ED visits (count)") +
  ggtitle("Asthma ED time series") + 
  geom_vline(xintercept=as.numeric(full_data_total$serv_dt[dates_vline]),
             linetype=4, colour="black")











