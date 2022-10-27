# "PHIRE Data Analysis - COPD ED model building"
# "Ruwan T"
# "5/11/2021"

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
library(sjPlot)
library(arm)
library(broom.mixed)


## Reading in data and creating covariates
data4 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")

data4$icd10 <- ifelse(data4$serv_dt>=ymd(20151001),1,0)

# Data header
head(data4)

#------------------------
# VISUALIZING TIME SERIES
#------------------------
data4$zip3 <- substr(data4$PATZIP,1,3)
locations <- c("San Francisco","Sacramento","Redding","Los Angeles","Fresno")
locations_zip3 <-  c(940,958,960,900,937)
data4_ts <- data4 %>% filter(zip3 %in% locations_zip3)
data4_ts$serv_dt <- as.Date(data4_ts$serv_dt)
for (i in 1:length(locations_zip3)){
  df <- data4_ts %>% filter(zip3==locations_zip3[i])
  df2 <- df %>% group_by(serv_dt)%>%summarize(copd_sum = sum(ed_copd))
  y <- ggplot(data=df2) + geom_point(aes(x=serv_dt,y=copd_sum)) + 
    xlab("Date") + ylab("Total visits") + ggtitle(paste(locations[i],"All-cause COPD ED (any dx) visits over time",sep=" ")) +
    scale_x_date(date_breaks = "1 year",date_labels = "%Y")
  print(y)
}

data_col <- data4 %>% group_by(serv_dt)%>%summarize(copd_sum = sum(ed_copd))
dates_vline <- as.Date(c("2015-10-01")) 
dates_vline <- which(data_col$serv_dt %in% dates_vline)
ggplot(data=data_col) + geom_point(aes(x=serv_dt,y=copd_sum)) + 
  xlab("Date") + ylab("Total visits") + ggtitle(paste("California","COPD ED (any dx) visits over time",sep=" ")) +
  scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
  geom_vline(xintercept=as.numeric(data_col$serv_dt[dates_vline]),
             linetype=4, colour="black")


#---------------
# MODEL BUILDING
#---------------


unname <- function(x) {
  names(x) <- NULL  
  x
}

# preallocating storage dataframe
results <- data.frame(matrix(NA,nrow=0,ncol=9))
names(results) <- c("Model","Model change","Time (s)","Std error","RR","95%CI LB","95%CI UB","AIC","BIC")


## Model 1

# pm25 residuals
tic <- Sys.time()
model_1 <- speedglm(ed_copd ~ pm25_residual_cz,data = data4,family=poisson(link="log"))
toc <- Sys.time()
model_1_results <- c("copd ~ pm25_residual_cz",
                     "exposure only",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_1)))[2])*10),8),
                     round(unname(exp(coef(model_1)[2]*10)),5),
                     round(unname(exp(confint.default(model_1)[2,1]*10)),5),
                     round(unname(exp(confint.default(model_1)[2,2]*10)),5),
                     prettyNum(round(unname(AIC(model_1)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_1)),2),big.mark=","))
results <- rbind(results,model_1_results)
summary(model_1)



## Model 2

# pm25 + zip code
tic <- Sys.time()
model_2 <- glmer(ed_copd ~ pm25_residual_cz + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_2_results <- c("copd ~ pm25_residual_cz + (1|PATZIP)",
                     "+ zip code RI",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_2)))[2])*10),8),
                     round(exp((unname(fixef(model_2))*10)[2]),5),
                     round(exp((unname(fixef(model_2))*10)[2] - 1.96*(unname(se.fixef(model_2))*10)[2]),5),
                     round(exp((unname(fixef(model_2))*10)[2] + 1.96*(unname(se.fixef(model_2))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_2)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_2)),2),big.mark=","))

results <- rbind(results,model_2_results)
summary(model_2)



## Model 3

# pm25 + spline(4/y) + zip code RI
tic <- Sys.time()
model_3 <- glmer(ed_copd ~ pm25_residual_cz + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_3_results <- c("copd ~ pm25_residual_cz + ns(dayseq,36) + (1|PATZIP)",
                     "+ time function",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_3)))[2])*10),8),
                     round(exp((unname(fixef(model_3))*10)[2]),5),
                     round(exp((unname(fixef(model_3))*10)[2] - 1.96*(unname(se.fixef(model_3))*10)[2]),5),
                     round(exp((unname(fixef(model_3))*10)[2] + 1.96*(unname(se.fixef(model_3))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_3)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_3)),2),big.mark=","))
results <- rbind(results,model_3_results)
summary(model_3)



## Model 4

# pm25 + DOW + spline(4/y) + zip code RI
tic <- Sys.time()
model_4 <- glmer(ed_copd ~ pm25_residual_cz + dow + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_4_results <- c("copd ~ pm25_residual_cz + dow + ns(dayseq,36) + (1|PATZIP)",
                     "+ DOW",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_4)))[2])*10),8),
                     round(exp((unname(fixef(model_4))*10)[2]),5),
                     round(exp((unname(fixef(model_4))*10)[2] - 1.96*(unname(se.fixef(model_4))*10)[2]),5),
                     round(exp((unname(fixef(model_4))*10)[2] + 1.96*(unname(se.fixef(model_4))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_4)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_4)),2),big.mark=","))
results <- rbind(results,model_4_results)
summary(model_4)



## Model 5


# pm25 + DOW + heat index + spline(4/y) + zip code RI
tic <- Sys.time()
model_5 <- glmer(ed_copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_5_results <- c("copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                     "+ heat index",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_5)))[2])*10),8),
                     round(exp((unname(fixef(model_5))*10)[2]),5),
                     round(exp((unname(fixef(model_5))*10)[2] - 1.96*(unname(se.fixef(model_5))*10)[2]),5),
                     round(exp((unname(fixef(model_5))*10)[2] + 1.96*(unname(se.fixef(model_5))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_5)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_5)),2),big.mark=","))
results <- rbind(results,model_5_results)
summary(model_5)

## Model 5's - dfs/year
# pm25 + DOW + heat index + spline(1/y) + zip code RI
# zips <- unique(data4$PATZIP)
# samp <- sample(zips,size=50)
# data4 <- data4 %>% filter(PATZIP %in% samp)
results_dfselection <- c()
dof <- c(18,36,45,54,63,72,81,90,99,108)
for (i in 1:length(dof)){
  tic <- Sys.time()
  model_5 <- glmer(ed_asthma ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,dof[i]) + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
  
  toc <- Sys.time()
  model_5_results <- c("asthma ~ pm25 + ns(hi,3) + dow + ns(dayseq) + (1|PATZIP)",
                       paste("+ heat index; spline ",dof[i]/9," df/yr",sep=""),
                       round(as.numeric(difftime(toc,tic),units="secs"),2),
                       round(unname((sqrt(diag(vcov(model_5)))[2])*10),8),
                       round(exp((unname(fixef(model_5))*10)[2]),5),
                       round(exp((unname(fixef(model_5))*10)[2] - 1.96*(unname(se.fixef(model_5))*10)[2]),5),
                       round(exp((unname(fixef(model_5))*10)[2] + 1.96*(unname(se.fixef(model_5))*10)[2]),5),
                       prettyNum(round(unname(AIC(model_5)),2),big.mark = ","),
                       prettyNum(round(unname(BIC(model_5)),2),big.mark=","))
  results_dfselection <- rbind(results_dfselection,model_5_results)
  print(i)
}
results_dfselection <- data.frame(results_dfselection)
names(results_dfselection) <- c("Model","Change","Time (s)", "Std error x 10","% change in risk","LB","UB","AIC","BIC")
write.csv(results_dfselection,"output/COPD ED/Time spline selection.csv")

## Model 6

# pm25 + DOW + heat index + spline(4/y) + ICD change + zip code RI
tic <- Sys.time()
model_6 <- glmer(ed_copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + icd10 + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_6_results <- c("copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + ICD-10 + (1|PATZIP)",
                     "+ ICD-10 introduction",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_6)))[2])*10),8),
                     round(exp((unname(fixef(model_6))*10)[2]),5),
                     round(exp((unname(fixef(model_6))*10)[2] - 1.96*(unname(se.fixef(model_6))*10)[2]),5),
                     round(exp((unname(fixef(model_6))*10)[2] + 1.96*(unname(se.fixef(model_6))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_6)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_6)),2),big.mark=","))
results <- rbind(results,model_6_results)
summary(model_6)

# pm25 + DOW + heat index + spline(4/y) + ICD change + zip code RI
tic <- Sys.time()
model_6 <- glmer(ed_copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + icd10 + (1|PATZIP),data = data4,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_6_results <- c("copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + ICD-10 + (1|PATZIP)",
                     "+ ICD-10 introduction",
                     round(as.numeric(difftime(toc,tic),units="secs"),2),
                     round(unname((sqrt(diag(vcov(model_6)))[2])*10),8),
                     round(exp((unname(fixef(model_6))*10)[2]),5),
                     round(exp((unname(fixef(model_6))*10)[2] - 1.96*(unname(se.fixef(model_6))*10)[2]),5),
                     round(exp((unname(fixef(model_6))*10)[2] + 1.96*(unname(se.fixef(model_6))*10)[2]),5),
                     prettyNum(round(unname(AIC(model_6)),2),big.mark = ","),
                     prettyNum(round(unname(BIC(model_6)),2),big.mark=","))
results <- rbind(results,model_6_results)
summary(model_6)




## Distributed lag comparisons

# create lags
data5 <- data4 %>%
  group_by(PATZIP) %>%
  mutate(pm25_residual_cz_lag1 = lag(pm25_residual_cz, n = 1),
         pm25_residual_cz_lag2 = lag(pm25_residual_cz, n = 2),
         pm25_residual_cz_lag3 = lag(pm25_residual_cz, n = 3),
         pm25_residual_cz_lag4 = lag(pm25_residual_cz, n = 4),
         pm25_residual_cz_lag5 = lag(pm25_residual_cz, n = 5),
         pm25_residual_cz_lag6 = lag(pm25_residual_cz, n = 6),
         pm25_residual_cz_lag7 = lag(pm25_residual_cz, n = 7),
         pm25_residual_cz_lag8 = lag(pm25_residual_cz, n = 8),
         pm25_residual_cz_lag9 = lag(pm25_residual_cz, n = 9),
         pm25_residual_cz_lag10 = lag(pm25_residual_cz, n = 10),
         pm25_residual_cz_lag11 = lag(pm25_residual_cz, n = 11),
         pm25_residual_cz_lag12 = lag(pm25_residual_cz, n = 12),
         pm25_residual_cz_lag13 = lag(pm25_residual_cz, n = 13),
         pm25_residual_cz_lag14 = lag(pm25_residual_cz, n = 14),
         pm25_residual_cz_lag15 = lag(pm25_residual_cz, n = 15),
         pm25_residual_cz_lag16 = lag(pm25_residual_cz, n = 16),
         pm25_residual_cz_lag17 = lag(pm25_residual_cz, n = 17),
         pm25_residual_cz_lag18 = lag(pm25_residual_cz, n = 18),
         pm25_residual_cz_lag19 = lag(pm25_residual_cz, n = 19),
         pm25_residual_cz_lag20 = lag(pm25_residual_cz, n = 20),
         pm25_residual_cz_lag21 = lag(pm25_residual_cz, n = 21),
         pm25_residual_cz_lag22 = lag(pm25_residual_cz, n = 22),
         pm25_residual_cz_lag23 = lag(pm25_residual_cz, n = 23),
         pm25_residual_cz_lag24 = lag(pm25_residual_cz, n = 24),
         pm25_residual_cz_lag25 = lag(pm25_residual_cz, n = 25),
         pm25_residual_cz_lag26 = lag(pm25_residual_cz, n = 26),
         pm25_residual_cz_lag27 = lag(pm25_residual_cz, n = 27),
         pm25_residual_cz_lag28 = lag(pm25_residual_cz, n = 28),
         pm25_residual_cz_lag29 = lag(pm25_residual_cz, n = 29),
         pm25_residual_cz_lag30 = lag(pm25_residual_cz, n = 30)
         
  )

data6 <- data5 %>% filter(serv_dt >= ymd(20080201))

data6$pm25_residual_cz_lag01 <- (data6$pm25_residual_cz+data6$pm25_residual_cz_lag1)/2
## DLM 0


# pm25 + DOW + heat index + spline(4/y) + ICD change + zip code RI
tic <- Sys.time()
model_lag0 <- glmer(ed_copd ~ pm25_residual_cz + ns(hi,3) + ns(dayseq,36) + dow + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_lag0_results <- c("copd ~ pm25_residual_cz + ns(hi,3) + dow + (1|PATZIP)",
                        "DLM 0",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag0)))[2])*10),8),
                        round(exp((unname(fixef(model_lag0))*10)[2]),5),
                        round(exp((unname(fixef(model_lag0))*10)[2] - 1.96*(unname(se.fixef(model_lag0))*10)[2]),5),
                        round(exp((unname(fixef(model_lag0))*10)[2] + 1.96*(unname(se.fixef(model_lag0))*10)[2]),5),
                        prettyNum(round(unname(AIC(model_lag0)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag0)),2),big.mark=","))
results <- rbind(results,model_lag0_results)
summary(model_lag0)

rr <- as.numeric(round(exp(unname(fixef(model_lag0)[2]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag0))*10)[2] - 1.96*(unname(se.fixef(model_lag0))*10)[2]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag0))*10)[2] + 1.96*(unname(se.fixef(model_lag0))*10)[2]),5))

lag <- c(0)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")


## DLM 0-1

model_lag1 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag1_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-1",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag1)))[2])*10),8),
                        round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag1)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag1)),2),big.mark=","))

results <- rbind(results,model_lag1_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag1)[2:3]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag1))*10)[2:3] - 1.96*(unname(se.fixef(model_lag1))*10)[2:3]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag1))*10)[2:3] + 1.96*(unname(se.fixef(model_lag1))*10)[2:3]),5))

lag <- c(0,1)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-1 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")

#lag0-1 MA
model_lag01 <- glmer(ed_copd ~ pm25_residual_cz_lag01 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag01_results <- c("copd ~pm25_residual_cz_lag01 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "MA 0-1",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag01)))[2])*10),8),
                        round(exp((unname(fixef(model_lag01))*10)[2]),5),
                        round(exp((unname(fixef(model_lag01))*10)[2] - 1.96*(unname(se.fixef(model_lag01))*10)[2]),5),
                        round(exp((unname(fixef(model_lag01))*10)[2] + 1.96*(unname(se.fixef(model_lag01))*10)[2]),5),
                        prettyNum(round(unname(AIC(model_lag01)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag01)),2),big.mark=","))


## DLM 0-2

model_lag2 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag2_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-2",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag2)))[2])*10),8),
                        round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag2)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag2)),2),big.mark=","))

results <- rbind(results,model_lag2_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag2)[2:4]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag2))*10)[2:4] - 1.96*(unname(se.fixef(model_lag2))*10)[2:4]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag2))*10)[2:4] + 1.96*(unname(se.fixef(model_lag2))*10)[2:4]),5))

lag <- c(0,1,2)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-2 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-3

model_lag3 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag3_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-3",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag3)))[2])*10),8),
                        round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag3)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag3)),2),big.mark=","))

results <- rbind(results,model_lag3_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag3)[2:5]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag3))*10)[2:5] - 1.96*(unname(se.fixef(model_lag3))*10)[2:5]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag3))*10)[2:5] + 1.96*(unname(se.fixef(model_lag3))*10)[2:5]),5))

lag <- c(0,1,2,3)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-3 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-4

model_lag4 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag4_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-4",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag4)))[2])*10),8),
                        round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag4)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag4)),2),big.mark=","))

results <- rbind(results,model_lag4_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag4)[2:6]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag4))*10)[2:6] - 1.96*(unname(se.fixef(model_lag4))*10)[2:6]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag4))*10)[2:6] + 1.96*(unname(se.fixef(model_lag4))*10)[2:6]),5))

lag <- c(0,1,2,3,4)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-4 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-5

model_lag5 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag5_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-5",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag5)))[2])*10),8),
                        round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag5)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag5)),2),big.mark=","))

results <- rbind(results,model_lag5_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag5)[2:7]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag5))*10)[2:7] - 1.96*(unname(se.fixef(model_lag5))*10)[2:7]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag5))*10)[2:7] + 1.96*(unname(se.fixef(model_lag5))*10)[2:7]),5))

lag <- c(0,1,2,3,4,5)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-5 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-6
model_lag6 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + pm25_residual_cz_lag6 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag6_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + pm25_residual_cz_lag6 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                        "DLM 0-6",
                        round(as.numeric(difftime(toc,tic),units="secs"),2),
                        round(unname((sqrt(diag(vcov(model_lag6)))[2])*10),8),
                        round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,1],5),
                        round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,3],5),
                        round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,4],5),
                        prettyNum(round(unname(AIC(model_lag6)),2),big.mark = ","),
                        prettyNum(round(unname(BIC(model_lag6)),2),big.mark=","))

results <- rbind(results,model_lag6_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag6)[2:8]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag6))*10)[2:8] - 1.96*(unname(se.fixef(model_lag6))*10)[2:8]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag6))*10)[2:8] + 1.96*(unname(se.fixef(model_lag6))*10)[2:8]),5))

lag <- c(0,1,2,3,4,5,6)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-6 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")

## DLM 0-14
model_lag14 <- glmer(ed_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + pm25_residual_cz_lag6 +
                       pm25_residual_cz_lag7 + pm25_residual_cz_lag8 + pm25_residual_cz_lag9+ pm25_residual_cz_lag10 + pm25_residual_cz_lag11 + pm25_residual_cz_lag12 + pm25_residual_cz_lag13 + pm25_residual_cz_lag14 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag14_results <- c("copd ~ pm25_residual_cz + pm25_residual_cz_lag1-14 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                         "DLM 0-14",
                         round(as.numeric(difftime(toc,tic),units="secs"),2),
                         round(unname((sqrt(diag(vcov(model_lag14)))[2])*10),8),
                         round(exp(deltaMethod(model_lag14, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10 + pm25_residual_cz_lag7*10 +
                                              pm25_residual_cz_lag8*10 + pm25_residual_cz_lag9*10 + pm25_residual_cz_lag10*10 + pm25_residual_cz_lag11*10 + pm25_residual_cz_lag12*10 + pm25_residual_cz_lag13*10 + pm25_residual_cz_lag14*10"))[1,1],5),
                         round(exp(deltaMethod(model_lag14, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10 + pm25_residual_cz_lag7*10 +
                                              pm25_residual_cz_lag8*10 + pm25_residual_cz_lag9*10 + pm25_residual_cz_lag10*10 + pm25_residual_cz_lag11*10 + pm25_residual_cz_lag12*10 + pm25_residual_cz_lag13*10 + pm25_residual_cz_lag14*10"))[1,3],5),
                         round(exp(deltaMethod(model_lag14, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10 + pm25_residual_cz_lag7*10 +
                                              pm25_residual_cz_lag8*10 + pm25_residual_cz_lag9*10 + pm25_residual_cz_lag10*10 + pm25_residual_cz_lag11*10 + pm25_residual_cz_lag12*10 + pm25_residual_cz_lag13*10 + pm25_residual_cz_lag14*10"))[1,4],5),
                         prettyNum(round(unname(AIC(model_lag14)),2),big.mark = ","),
                         prettyNum(round(unname(BIC(model_lag14)),2),big.mark=","))

results <- rbind(results,model_lag14_results)

rr <- as.numeric(round(exp(unname(fixef(model_lag14)[2:16]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag14))*10)[2:16] - 1.96*(unname(se.fixef(model_lag14))*10)[2:16]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag14))*10)[2:16] + 1.96*(unname(se.fixef(model_lag14))*10)[2:16]),5))

lag <- seq(0,14)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-14 for association between PM2.5 (10-unit increase) and COPD ED visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")




## PM2.5 coefficient comparisons

names(results) <- c("Model","Model change","Time (s)","Std error x 10","RR","95%CI LB","95%CI UB","AIC","BIC")
results2 <- cbind(paste("Model ",c("1","2","3","4","5","6",paste("DLM",seq(0,6))),sep=""),results[,c(1,2,4,5,6,7,9)])
                  kable(results2)
                  
                  
## BASELINE INCIDENCE calculation
                  