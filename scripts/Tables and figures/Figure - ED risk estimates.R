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
library(dlnm)
library(flextable)

results <- readRDS("output/Risk estimate models/ed_all_results.rds")
# results$indicator <- 0
results$indicator <- ifelse(results$Subgroup=="NH Nat Am",1,0)

# line range
asthma_plot <- ggplot(data=results%>%filter(Outcome=="Asthma (Lag 0-2)",indicator==0)) +
  geom_point(aes(x=Subgroup,y=RR)) +
  geom_errorbar(aes(x=Subgroup,ymin=LB,ymax=UB),width=0.3) +
  facet_wrap(~ Outcome, ncol=1,scales="free") +
  ylab("Change in risk (%)") +
  # geom_hline(aes(yintercept=-2),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=0),linetype="solid") +
  geom_hline(aes(yintercept=1),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=2),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=3),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=4),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=5),linetype="dashed",color="grey") +
  scale_y_continuous(name ="Change in risk (%)", 
                     breaks=c(-2,0,2,4,6)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(-1,0,1,2,3,4,5,6))


copd_plot <- ggplot(data=results%>%filter(Outcome=="COPD (Lag 0-1)",indicator==0)) +
  geom_point(aes(x=Subgroup,y=RR)) +
  geom_errorbar(aes(x=Subgroup,ymin=LB,ymax=UB),width=0.3) +
  facet_wrap(~ Outcome, ncol=1,scales="free") +
  ylab("Change in risk (%)") +
  geom_hline(aes(yintercept=-1),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=0),linetype="solid") +
  geom_hline(aes(yintercept=2),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=4),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=6),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=8),linetype="dashed",color="grey") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(-1,0,2,4,6,8))

resp_plot <- ggplot(data=results%>%filter(Outcome=="Respiratory* (Lag 0-4)",indicator==0)) +
  geom_point(aes(x=Subgroup,y=RR)) +
  geom_errorbar(aes(x=Subgroup,ymin=LB,ymax=UB),width=0.3) +
  facet_wrap(~ Outcome, ncol=1,scales="free") +
  ylab("Change in risk (%)") +
  # geom_hline(aes(yintercept=-1),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=0),linetype="solid") +
  geom_hline(aes(yintercept=1),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=2),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=3),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=4),linetype="dashed",color="grey") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(-1,0,1,2,3,4,5,6))

cv_plot <- ggplot(data=results%>%filter(Outcome=="Cardiovascular** (Lag 0)",indicator==0)) +
  geom_point(aes(x=Subgroup,y=RR)) +
  geom_errorbar(aes(x=Subgroup,ymin=LB,ymax=UB),width=0.3) +
  facet_wrap(~ Outcome, ncol=1,scales="free") +
  ylab("Change in risk (%)") +
  # geom_hline(aes(yintercept=-1),linetype="dashed",color="grey") +
  
  geom_hline(aes(yintercept=0),linetype="solid") +
  geom_hline(aes(yintercept=1),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=2),linetype="dashed",color="grey") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(-1,0,1,2,3,4,5,6))


asthma_plot <- arrangeGrob(asthma_plot, top = textGrob("B", x = unit(0, "npc")
                                                       , y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=12, fontfamily="Times Roman")))
copd_plot <- arrangeGrob(copd_plot, top = textGrob("C", x = unit(0, "npc")
                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=12, fontfamily="Times Roman")))
resp_plot <- arrangeGrob(resp_plot, top = textGrob("A", x = unit(0, "npc")
                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=12, fontfamily="Times Roman")))
cv_plot <- arrangeGrob(cv_plot, top = textGrob("D", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12, fontfamily="Times Roman")))

png(filename="output/Paper 1 - Tables and Figures/Figure - ED risk estimates/Figure - ED risk estimates.png",units="in",width=11,height=8.5,res=1200)
grid.arrange(resp_plot, asthma_plot, copd_plot, cv_plot, ncol = 1)
dev.off()


