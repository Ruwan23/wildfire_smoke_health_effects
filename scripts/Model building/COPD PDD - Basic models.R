
# Ruwan Thilakaratne
# COPD PDD DLM models

memory.limit(size=100000)
library(arm)
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
library(weathermetrics)
library(lubridate)
library(effects)
library(sjPlot)
library(arm)
library(broom.mixed)

data4 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")

# restricting to after Feb 1 2008 (i.e. removing Jan 2008) to make BIC's comparable between models with different lag periods
data6 <- data4 %>% filter(serv_dt >= ymd(20080201))

unname <- function(x) {
  names(x) <- NULL  
  x
}


#----------
# COPD
#----------
# preallocating storage dataframe
results_copd <- data.frame(matrix(NA,nrow=0,ncol=9))
names(results_copd) <- c("Model","Model change","Time (s)","Std error","RR","95%CI LB","95%CI UB","AIC","BIC")


## DLM 0

# pm25 + DOW + heat index + spline(4/y) + ICD change + zip code RI
tic <- Sys.time()
model_lag0 <- glmer(pdd_copd ~ pm25_residual_cz + ns(hi,3) + ns(dayseq,36) + dow + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()
model_lag0_results_copd <- c("pdd_copd ~ pm25_residual_cz + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag0)))[2])*10),8),
                             round(exp((unname(fixef(model_lag0))*10)[2]),5),
                             round(exp((unname(fixef(model_lag0))*10)[2] - 1.96*(unname(se.fixef(model_lag0))*10)[2]),5),
                             round(exp((unname(fixef(model_lag0))*10)[2] + 1.96*(unname(se.fixef(model_lag0))*10)[2]),5),
                             prettyNum(round(unname(AIC(model_lag0)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag0)),2),big.mark=","))
results_copd <- rbind(results_copd,model_lag0_results_copd)
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
  ggtitle("DLM 0 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")


## DLM 0-1

model_lag1 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag1_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-1",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag1)))[2])*10),8),
                             round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag1, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag1)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag1)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag1_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag1)[2:3]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag1))*10)[2:3] - 1.96*(unname(se.fixef(model_lag1))*10)[2:3]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag1))*10)[2:3] + 1.96*(unname(se.fixef(model_lag1))*10)[2:3]),5))

lag <- c(0,1)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-1 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-2

model_lag2 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag2_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-2",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag2)))[2])*10),8),
                             round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag2, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag2)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag2)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag2_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag2)[2:4]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag2))*10)[2:4] - 1.96*(unname(se.fixef(model_lag2))*10)[2:4]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag2))*10)[2:4] + 1.96*(unname(se.fixef(model_lag2))*10)[2:4]),5))

lag <- c(0,1,2)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-2 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-3

model_lag3 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag3_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-3",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag3)))[2])*10),8),
                             round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag3, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag3)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag3)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag3_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag3)[2:5]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag3))*10)[2:5] - 1.96*(unname(se.fixef(model_lag3))*10)[2:5]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag3))*10)[2:5] + 1.96*(unname(se.fixef(model_lag3))*10)[2:5]),5))

lag <- c(0,1,2,3)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-3 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-4

model_lag4 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag4_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-4",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag4)))[2])*10),8),
                             round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag4, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag4)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag4)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag4_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag4)[2:6]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag4))*10)[2:6] - 1.96*(unname(se.fixef(model_lag4))*10)[2:6]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag4))*10)[2:6] + 1.96*(unname(se.fixef(model_lag4))*10)[2:6]),5))

lag <- c(0,1,2,3,4)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-4 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-5

model_lag5 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag5_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-5",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag5)))[2])*10),8),
                             round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag5, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag5)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag5)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag5_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag5)[2:7]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag5))*10)[2:7] - 1.96*(unname(se.fixef(model_lag5))*10)[2:7]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag5))*10)[2:7] + 1.96*(unname(se.fixef(model_lag5))*10)[2:7]),5))

lag <- c(0,1,2,3,4,5)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-5 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")



## DLM 0-6
model_lag6 <- glmer(pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + pm25_residual_cz_lag6 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP),data = data6,family=poisson(link="log"),nAGQ = 0)
toc <- Sys.time()

model_lag6_results_copd <- c("pdd_copd ~ pm25_residual_cz + pm25_residual_cz_lag1 + pm25_residual_cz_lag2 + pm25_residual_cz_lag3 + pm25_residual_cz_lag4 + pm25_residual_cz_lag5 + pm25_residual_cz_lag6 + ns(hi,3) + dow + ns(dayseq,36) + (1|PATZIP)",
                             "DLM 0-6",
                             round(as.numeric(difftime(toc,tic),units="secs"),2),
                             round(unname((sqrt(diag(vcov(model_lag6)))[2])*10),8),
                             round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,1],5),
                             round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,3],5),
                             round(exp(deltaMethod(model_lag6, "pm25_residual_cz*10 + pm25_residual_cz_lag1*10 + pm25_residual_cz_lag2*10 + pm25_residual_cz_lag3*10 + pm25_residual_cz_lag4*10 + pm25_residual_cz_lag5*10 + pm25_residual_cz_lag6*10"))[1,4],5),
                             prettyNum(round(unname(AIC(model_lag6)),2),big.mark = ","),
                             prettyNum(round(unname(BIC(model_lag6)),2),big.mark=","))

results_copd <- rbind(results_copd,model_lag6_results_copd)

rr <- as.numeric(round(exp(unname(fixef(model_lag6)[2:8]*10)),5))
lb <- as.numeric(round(exp((unname(fixef(model_lag6))*10)[2:8] - 1.96*(unname(se.fixef(model_lag6))*10)[2:8]),5))
ub <- as.numeric(round(exp((unname(fixef(model_lag6))*10)[2:8] + 1.96*(unname(se.fixef(model_lag6))*10)[2:8]),5))

lag <- c(0,1,2,3,4,5,6)
toplot <- data.frame(rr,lb,ub,lag)
names(toplot) <- c("RR","LB","UB","Lag")

ggplot(data=toplot,aes(x=Lag,y=RR)) +
  geom_point() +
  geom_errorbar(aes(ymin=LB,ymax=UB)) +
  ggtitle("DLM 0-6 for association between PM2.5 (10-unit increase) and copd PDD visits") +
  # ylim(0.99,1.01) +
  geom_hline(yintercept=1.0, linetype="dotted")


names(results_copd) <- c("Model","Model change","Time (s)","Std error","RR","95%CI LB","95%CI UB","AIC","BIC")
write.csv(results_copd,"output/COPD PDD/COPD PDD DLM models.csv")
