# Ruwan Thilakaratne
# exploring low PM in 2016
library(dplyr)
library(sp)
library(sf)
library(rgeos)
library(s2)
library(Rcpp)
library(lubridate)

# reading in dataset
data4 <- readRDS("S:/Wildfire-Hoshiko/PHIRE OSHPD DATA/Analytic datasets/FINAL ED and PDD dataset with covariates/Asthma, CV, Resp, COPD ED and PDD data strats with covariates and HMS.rds")
data4$year <- year(data4$serv_dt) # creating year variable to use later

# identifying low vs. high non-smoke PM zips by year
pm_cat <- data4 %>% filter(light==0 & medium==0 & heavy==0 & miss_or_na==0) # subsetting to non-smoke days
pm_cat <- pm_cat %>% group_by(PATZIP,year) %>% summarize(pm_mean = mean(pm25)) # taking zip-year PM2.5 means
pm_cat <- pm_cat %>% mutate(high_low = ifelse(pm_mean>12,1,0)) # indicator for above/below threshold
data5 <- inner_join(pm_cat,data4,by=c("PATZIP","year")) # merging indicator back in with original dataset

# analyzing low non-smoke PM zips
low <- data5 %>% filter(high_low==0,light==0 & medium==0 & heavy==0 & miss_or_na==0) # subsetting to low non-smoke PM zip non-smoke days
low_zip <- low %>% group_by(PATZIP,year) %>% summarize(pm25 = mean(pm25)) # zip-year means
low_zip$year <- as.factor(low_zip$year)
ggplot(data=low_zip) +  # plotting zip-level non-smoke PM means by year
  geom_boxplot(aes(y=pm25,group=year)) + 
  # scale_x_discrete(breaks = c(seq(2008,2016))) +
  labs(x="Year",y="ZIP-level annual average non-smoke PM2.5",title="ZIP-level annual average non-smoke PM2.5 by study year, among ZIPs with LOW non-smoke PM2.5")

# repeating for high non-smoke PM zips
high <- data5 %>% filter(high_low==1,light==0 & medium==0 & heavy==0 & miss_or_na==0)
high_zip <- high %>% group_by(PATZIP,year) %>% summarize(pm25 = mean(pm25))
ggplot(data=high_zip) + 
  geom_boxplot(aes(y=pm25,group=year)) + 
  labs(x="Year",y="ZIP-level annual average non-smoke PM2.5",title="ZIP-level annual average non-smoke PM2.5 by study year, among ZIPs with HIGH non-smoke PM2.5")

low <- data5 %>% filter(high_low==0,(light==1 | medium==1 | heavy==1))
low_zip <- low %>% group_by(PATZIP,year) %>% summarize(pm25 = mean(pm25))
low_zip$year <- as.factor(low_zip$year)
ggplot(data=low_zip) + 
  geom_boxplot(aes(y=pm25,group=year)) + 
  # scale_x_discrete(breaks = c(seq(2008,2016))) +
  labs(x="Year",y="ZIP-level annual average smoke PM2.5",title="ZIP-level annual average smoke PM2.5 by study year, among ZIPs with LOW non-smoke PM2.5")

high <- data5 %>% filter(high_low==1,(light==1 | medium==1 | heavy==1))
high_zip <- high %>% group_by(PATZIP,year) %>% summarize(pm25 = mean(pm25))
ggplot(data=high_zip) + 
  geom_boxplot(aes(y=pm25,group=year)) + 
  labs(x="Year",y="ZIP-level annual average smoke PM2.5",title="ZIP-level annual average smoke PM2.5 by study year, among ZIPs with HIGH non-smoke PM2.5")


# daily time series (averaging across zips)
daily <- data5 %>% filter(high_low==0,light==0 & medium==0 & heavy==0 & miss_or_na==0)

# (excluding 2008 because it's an anomalously polluted year and throws off y axis scale)
ggplot(data=daily %>% filter(high_low==0,year(serv_dt)>2008)) + geom_point(aes(x=serv_dt,y=pm_avg)) + labs(x="date",y="PM2.5 average for CA",title="Daily average PM2.5 in California, among ZIPs with LOW non-smoke PM2.5")
ggplot(data=daily %>% filter(high_low==1,year(serv_dt)>2008)) + geom_point(aes(x=serv_dt,y=pm_avg)) + labs(x="date",y="PM2.5 average for CA",title="Daily average PM2.5 in California, among ZIPs with HIGH non-smoke PM2.5")

daily2 <- data5 %>% group_by(serv_dt, high_low) %>% summarize(pm_avg = mean(pm25))
ggplot(data=daily2 %>% filter(high_low==1,year(serv_dt)>2008)) + 
  geom_point(aes(x=serv_dt,y=pm_avg)) + 
  labs(x="date",y="TOTAL PM2.5 average for CA",title="Daily average TOTAL PM2.5 in California, among ZIPs with HIGH non-smoke PM2.5 (>12)")
ggplot(data=daily2 %>% filter(high_low==0,year(serv_dt)>2008)) + 
  geom_point(aes(x=serv_dt,y=pm_avg)) + 
  labs(x="date",y="TOTAL PM2.5 average for CA",title="Daily average TOTAL PM2.5 in California, among ZIPs with LOW non-smoke PM2.5 (<=12)")

high <- data5 %>% filter(high_low==1)
high_year_zips <- high %>% group_by(PATZIP,year) %>% summarize(n = n())
zips_by_year <- high_year_zips %>% group_by(year) %>% summarize(n = n())

# mapping high vs low PM areas, and PM2.5 means, by year for 2015 and 2016 to compare
# getting zip boundaries and processing for merge
ca_zip3 <- seq(900,961,1)
ca_zip_map <- st_read("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/CA zips/tl_2019_us_zcta510/tl_2019_us_zcta510.shp")
ca_zip_map <- ca_zip_map %>% dplyr::filter(substr(GEOID10,1,3) %in% ca_zip3)
ca_zip_map <- ca_zip_map %>% rename("PATZIP" = "GEOID10")
ca_zip_map$PATZIP <- as.numeric(ca_zip_map$PATZIP)

# zip-level means of 2015 non-smoke days
zips_2015 <- data5 %>% filter(year==2015,(light==0 & medium==0 & heavy==0 & miss_or_na==0))
zips_2015 <- zips_2015 %>% group_by(PATZIP) %>% summarize(high_low = mean(high_low),
                                                          pm_mean = mean(pm25))
map_2015 <- inner_join(zips_2015,ca_zip_map,by=c("PATZIP"))
ggplot() + 
  geom_sf(data = map_2015, aes(geometry=geometry,fill=factor(high_low)),color=NA,lwd=0) + 
  ggtitle("High vs. low non-smoke PM (>12 or <= 12) by zip code, 2015")
ggplot() + 
  geom_sf(data = map_2015, aes(geometry=geometry,fill=pm_mean),color=NA,lwd=0) + 
  ggtitle("Mean NON-SMOKE PM2.5 by zip code, 2015") + 
  scale_fill_gradientn(colours=heat.colors(4),
                       breaks=c(5,10,15,20),labels=c("0-5","5-10","10-15","15-20"),
                       limits=c(0,20))

zips_2016 <- data5 %>% filter(year==2016,(light==0 & medium==0 & heavy==0 & miss_or_na==0))
zips_2016 <- zips_2016 %>% group_by(PATZIP) %>% summarize(high_low = mean(high_low),
                                                          pm_mean = mean(pm25))
map_2016 <- inner_join(zips_2016,ca_zip_map,by=c("PATZIP"))
ggplot() + 
  geom_sf(data = map_2016, aes(geometry=geometry,fill=factor(high_low)),color=NA,lwd=0) + 
  ggtitle("High vs. low non-smoke PM (>12 or <= 12) by zip code, 2016")
ggplot() + 
  geom_sf(data = map_2016, aes(geometry=geometry,fill=pm_mean),color=NA,lwd=0) + 
  ggtitle("Mean NON-SMOKE PM2.5 by zip code, 2016") + 
  scale_fill_gradientn(colours=heat.colors(4),
                       breaks=c(5,10,15,20),labels=c("0-5","5-10","10-15","15-20"),
                       limits=c(0,20))

# histogram of zip-level annual averages in 2015 vs 2016
no_smoke_2015 <- data5 %>% filter(year==2015, (light==0 & medium==0 & heavy==0 & miss_or_na==0))
no_smoke_2015_zip <- no_smoke_2015 %>% group_by(PATZIP) %>% summarize(pm_mean = mean(pm25))
ggplot(data=no_smoke_2015_zip) + 
  geom_histogram(aes(x=pm_mean),color="black",fill="white",binwidth=1) + 
  scale_x_continuous(breaks=c(seq(0,25,5))) +
  labs(x="Non-smoke PM2.5 annual average (zip-level)",title="Histogram of zip-level non-smoke PM2.5 annual average in 2015")

no_smoke_2016 <- data5 %>% filter(year==2016, (light==0 & medium==0 & heavy==0 & miss_or_na==0))
no_smoke_2016_zip <- no_smoke_2016 %>% group_by(PATZIP) %>% summarize(pm_mean = mean(pm25))
ggplot(data=no_smoke_2016_zip) + 
  geom_histogram(aes(x=pm_mean),color="black",fill="white",binwidth=1) + 
  scale_x_continuous(breaks=c(seq(0,25,5))) +
  labs(x="Non-smoke PM2.5 annual average (zip-level)",title="Histogram of zip-level non-smoke PM2.5 annual average in 2016")

# FRM correlations
ca_zip3 <- seq(900,961,1)
frm_16_19 <- read.csv("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data for Analysis/cleaned_daily_FRM_data_2016-2019.csv")
frm_10_15 <- read.csv("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data for Analysis/cleaned_daily_FRM_data_2010-2015.csv") 
frm <- rbind(frm_10_15,frm_16_19) # concatenating the datasets (have same variable names and types)
centroids <- read.csv("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/Data for Analysis/2015_zip_centroids.csv") # 2015 zip centroids from census website
frm_ca <- frm %>% 
          filter(State=="California") %>%
          rename(lon = Lon,
                 lat = Lat) %>%
          mutate(lat_lon = paste(lat,lon,sep="_"))
frm_ca <- data.frame(frm_ca)
frm_ca2 <- frm_ca %>%
          group_by(lat,lon) %>%
          summarize(n = n())
          
centroids_ca <- centroids %>% 
              filter(substr(zcta5,1,3) %in% ca_zip3) %>%
              rename(lat = intptlat,
                     lon = intptlong)
# geoprocessing to link FRMs (n=117) with nearest zip centroid
set1_s2 <- s2_lnglat(frm_ca2$lon, frm_ca2$lat)
set2_s2 <- s2_lnglat(centroids_ca$lon, centroids_ca$lat)
frm_ca2$closest <- s2_closest_feature(set1_s2, set2_s2)
test <- frm_ca2
test$close <- s2_distance_matrix(set1_s2, set2_s2, radius = s2_earth_radius_meters())
centroids_ca$rown <- seq(1,nrow(centroids_ca))
centroids_merge <- centroids_ca %>% dplyr::select(zcta5,rown) %>% rename(closest = rown)
frm_ca3 <- inner_join(frm_ca2, centroids_merge, by=c("closest"))
# merging in temporal frm data
frm_final <- inner_join(frm_ca3, frm_ca, by=c("lat","lon"))
# merging FRM to analytic data on zip and date, subsetting to zips in common (using inner join)
frm_final <- frm_final %>% 
                rename(PATZIP = zcta5,
                       serv_dt = Date) %>%
                mutate(serv_dt = as.Date(serv_dt))
data4_frm <- inner_join(frm_final, data4, by=c("PATZIP","serv_dt"))
data4_frm <- data4_frm %>%
                rename(pm25_frm = PM25,
                       pm25_harv = pm25)
data4_frm_nosmoke <- data4_frm %>% filter(light==0 & medium==0 & heavy==0 & miss_or_na==0)
data4_frm_nosmoke$year <- year(data4_frm_nosmoke$serv_dt)
# all years
correlations_zip <- data4_frm_nosmoke %>%
                  group_by(PATZIP) %>%
                  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
                  mutate(year = "all")
correlations_region <- data4_frm_nosmoke %>%
                  group_by(region) %>%
                  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
                  mutate(year = "all")


# 2011
data4_frm_nosmoke_2011 <- data4_frm_nosmoke %>% filter(year==2011)
correlations_zip_2011 <- data4_frm_nosmoke_2011 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2011")
correlations_region_2011 <- data4_frm_nosmoke_2011 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2011")

# 2012
data4_frm_nosmoke_2012 <- data4_frm_nosmoke %>% filter(year==2012)
correlations_zip_2012 <- data4_frm_nosmoke_2012 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2012")
correlations_region_2012 <- data4_frm_nosmoke_2012 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2012")

# 2013
data4_frm_nosmoke_2013 <- data4_frm_nosmoke %>% filter(year==2013)
correlations_zip_2013 <- data4_frm_nosmoke_2013 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2013")
correlations_region_2013 <- data4_frm_nosmoke_2013 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2013")

# 2014
data4_frm_nosmoke_2014 <- data4_frm_nosmoke %>% filter(year==2014)
correlations_zip_2014 <- data4_frm_nosmoke_2014 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2014")
correlations_region_2014 <- data4_frm_nosmoke_2014 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2014")

# 2015
data4_frm_nosmoke_2015 <- data4_frm_nosmoke %>% filter(year==2015)
correlations_zip_2015 <- data4_frm_nosmoke_2015 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2015")
correlations_region_2015 <- data4_frm_nosmoke_2015 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2015")

# 2016
data4_frm_nosmoke_2016 <- data4_frm_nosmoke %>% filter(year==2016)
correlations_zip_2016 <- data4_frm_nosmoke_2016 %>%
  group_by(PATZIP) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2016")
correlations_region_2016 <- data4_frm_nosmoke_2016 %>%
  group_by(region) %>%
  summarize(cor = cor(pm25_frm,pm25_harv)) %>%
  mutate(year = "2016")

# boxplot
all_cors <- rbind(correlations_zip, correlations_zip_2011, correlations_zip_2012, correlations_zip_2013, correlations_zip_2014, correlations_zip_2015, correlations_zip_2016)
ggplot(data=all_cors) + geom_boxplot(aes(x=year,y=cor)) + labs(x="Time period",y="Correlation coefficient (r)",title="Distributions of ZIP code-level correlation coefficients between FRM and Harvard PM2.5 on non-smoke days by time period")

# map
ca_county_map <- st_read("S:/Wildfire-Hoshiko/PHIRE COVARIATE DATA/CA_Counties/CA_Counties_TIGER2016.shp")
ca_zip_map_cors <- inner_join(ca_zip_map, correlations_zip, by=c("PATZIP"))
ggplot() + 
  geom_sf(data = ca_county_map, aes(geometry=geometry)) + 
  geom_sf(data = ca_zip_map_cors, aes(geometry=geometry,fill=cor),color=NA,lwd=0) + 
  ggtitle("Correlations between Harvard and FRM by zip code, whole study period")

ca_zip_map_cors_2015 <- inner_join(ca_zip_map, correlations_zip_2015, by=c("PATZIP"))
ggplot() + 
  geom_sf(data = ca_county_map, aes(geometry=geometry)) + 
  geom_sf(data = ca_zip_map_cors_2015, aes(geometry=geometry,fill=cor),color=NA,lwd=0) + 
  ggtitle("Correlations between Harvard and FRM by zip code, 2015")

ca_zip_map_cors_2016 <- inner_join(ca_zip_map, correlations_zip_2016, by=c("PATZIP"))
ggplot() + 
  geom_sf(data = ca_county_map, aes(geometry=geometry)) + 
  geom_sf(data = ca_zip_map_cors_2016, aes(geometry=geometry,fill=cor),color=NA,lwd=0) + 
  ggtitle("Correlations between Harvard and FRM by zip code, 2016")



