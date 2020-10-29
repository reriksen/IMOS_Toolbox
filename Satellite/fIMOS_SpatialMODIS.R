# THis software downloads IMOS satellite spatial data given lat/lon limits 
# and a time frame using the OPeNDAP protocol.
# A nice website with a description is here:
# https://opendap.github.io/documentation/QuickStart.html
#
# Inputs: 
# dat is a dataframe or tibble with the following columns:
# Lat_Min: Min Latitude
# Lat_Max: Max Latitude
# Lon_Min: Min Longitude
# Lon_Max: Max Latitude
# Date_Min: Min Date (as POSIX)
# Date_Max: Max Date (as POSIX)
# Prod: are the MODIS products you wish to extract. See below for a list of options
# TempRes: Temporal resolution - What temporal averaging: 1 day (1d), 1 month(1m), 1 year (1y),...
# Monthly climatology (1mNy), Annual climatology (12mNy)

# Possible products to download are: 
# "sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
# "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
# "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490"


# fIMOS_MatchMODIS <- function(dat, pr, ...) {
library(stringr) #install.packages("stringr")  
library(dplyr) #install.packages("dplyr")
library(lubridate) #install.packages("lubridate")
library(raster) #install.packages("raster")
library(tidync)

dat <- tibble(Lat_Min = -50, 
              Lat_Max = 0,
              Lon_Min = 100,
              Lon_Max = 170,
              Date_Min = ymd("2003-01-01"),
              Date_Max = ymd("2020-06-30"),
              Prod = "chl_oc3",
              TempRes = "1m")

n_mths <- interval(dat$Date_Min,dat$Date_Max)%/%months(1)
all_dates <- dat$Date_Min + months(0:n_mths)

pb <- txtProgressBar(min = 1, max = length(all_dates), style = 3)
for (i in 1:length(all_dates)) { # Loop through all rows in the data for each variable you want
  
  # Make sure month and day have a leading zero if less than 10
  mth <- str_pad(month(all_dates[i]),2,"left", pad="0")
  
  url_base <- paste0("http://rs-data1-mel.csiro.au/thredds/dodsC/imos-srs/oc/aqua/",dat$TempRes,"/") # Base URL
  
  if (str_detect(dat$TempRes,"1m")){
    imos_url <- paste0(url_base, year(all_dates[i]),"/",year(all_dates[i]),mth,".",dat$Prod,".nc4") # nc4 starts in October 2017
    vr <- paste0(dat$Prod,"_mean")
  }
  
  out <- tidync(imos_url) %>%
    hyper_filter(longitude = longitude >= dat$Lon_Min & longitude <= dat$Lon_Max,
                 latitude = latitude >= dat$Lat_Min & latitude <= dat$Lat_Max) %>%
    hyper_array(select_var = "chl_oc3_mean")
  
  if(i == 1){
    r <- raster(t(out$chl_oc3_mean), xmn=dat$Lon_Min, xmx=dat$Lon_Max, ymn=dat$Lat_Min, ymx=dat$Lat_Max, crs = "+proj=longlat +datum=WGS84")
  } else{
    r <- stack(r, raster(t(out$chl_oc3_mean), xmn=dat$Lon_Min, xmx=dat$Lon_Max, ymn=dat$Lat_Min, ymx=dat$Lat_Max, crs = "+proj=longlat +datum=WGS84"))
  }
  
  setTxtProgressBar(pb, i)
}
names(r) <- all_dates[1:192]


myRaster <- writeRaster(r,"myStack.grd", format="raster")

lm_fun = function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }} # slope

slope <- calc(r, lm_fun)


