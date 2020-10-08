# Downloads GSLA, GSL, u and v from CSIRO Blulink DODS server.
# http)//www.marine.csiro.au/dods/nph-dods/dods-data/bluelink/GSLA/'
#
# Useage: Alt = sat_plot_Altimetry(date_num,Lat,Lon,count,bkgrnd,coast)
# where:
#       date_num: matlab datenum
#       Lat and Lon: 2 value vectors
#       Alt and count: not used yet.
#
# U - East-West Velocity
# V - Noth-South Velocity
#
# Written by Jason Everett (UQ) July 2019


fIMOS_MatchAltimetry <- function(dat, ...) {
  library(stringr) #install.packages("stringr")  
  library(dplyr) #install.packages("dplyr")
  library(tibble) #install.packages("tibble")
  library(ncdf4) #install.packages("ncdf4")
  library(yaImpute) #install.packages("yaImpute")
  library(lubridate) #install.packages("lubridate")
  
  dat <- dat %>% 
    add_column(GSLA = NA, GSL = NA, UCUR = NA, VCUR = NA, Longitude_Actual = NA, Latitude_Actual = NA)
  
  if (!exists("res_spat")){
    print("Defaulting to 1 pixel x 1 pixel. Provide res_spat if you want to increase")
    res_spat <-  1
  }
  
  if (sum(c("Day","Month","Year") %in% colnames(dat)) != 3) { # Check that Day, Month, Year exists
    if (sum(str_detect(colnames(dat),"Date")) == 1) { # Otherwise check that Date exists
      dat <- dat %>% 
        mutate(Day = day(Date),
               Month = month(Date),
               Year = year(Date))
    } else {
      print("Missing Date or Day/Month/Year columns")
    }
  }
  
  
  pb <- txtProgressBar(min = 0, max = length(dat$Latitude), style = 3)
  for (i in 1:length(dat$Latitude)) { # Loop through all rows in the data for each variable you want
    
    # Make sure month and day have a leading zero if less than 10
    mth <- str_pad(dat$Month[i],2,"left",pad="0")
    dy <- str_pad(dat$Day[i],2,"left",pad="0")
    yr <- dat$Year[i]
    
    url_base = 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/OceanCurrent/GSLA/DM01/yearfiles/';
    # 
    # if (yr==2019) {file = 'IMOS_OceanCurrent_HV_2019_C-20190605T233144Z.nc.gz'}
    # if (yr==2018) {file = 'IMOS_OceanCurrent_HV_2018_C-20190129T223442Z.nc.gz'}
    # if (yr==2017) {file = 'IMOS_OceanCurrent_HV_2017_C-20180129T224442Z.nc.gz'}	
    # if (yr==2016) {file = 'IMOS_OceanCurrent_HV_2016_C-20170129T222646Z.nc.gz'}
    # if (yr==2015) {file = 'IMOS_OceanCurrent_HV_2015_C-20160903T062732Z.nc.gz'}
    # if (yr==2014) {file = 'IMOS_OceanCurrent_HV_2014_C-20150413T032427Z.nc.gz'}
    # if (yr==2013) {file = 'IMOS_OceanCurrent_HV_2013_C-20141028T045015Z.nc.gz'}
    # if (yr==2012) {file = 'IMOS_OceanCurrent_HV_2012_C-20140106T022022Z.nc.gz'}
    # if (yr==2011) {file = 'IMOS_OceanCurrent_HV_2011_C-20141028T044752Z.nc.gz'}
    # if (yr==2010) {file = 'IMOS_OceanCurrent_HV_2010_C-20140120T231433Z.nc.gz'}
    # if (yr==2009) {file = 'IMOS_OceanCurrent_HV_2009_C-20140120T232133Z.nc.gz'}
    # if (yr==2008) {file = 'IMOS_OceanCurrent_HV_2008_C-20150521T043304Z.nc.gz'}
    # if (yr==2007) {file = 'IMOS_OceanCurrent_HV_2007_C-20150521T042748Z.nc.gz'}
    # if (yr==2006) {file = 'IMOS_OceanCurrent_HV_2006_C-20150521T042256Z.nc.gz'}
    # if (yr==2005) {file = 'IMOS_OceanCurrent_HV_2005_C-20150521T041806Z.nc.gz'}
    # if (yr==2004) {file = 'IMOS_OceanCurrent_HV_2004_C-20150521T041319Z.nc.gz'}
    # if (yr==2003) {file = 'IMOS_OceanCurrent_HV_2003_C-20150521T040831Z.nc.gz'}
    # if (yr==2002) {file = 'IMOS_OceanCurrent_HV_2002_C-20150521T040330Z.nc.gz'}
    # if (yr==2001) {file = 'IMOS_OceanCurrent_HV_2001_C-20150521T035914Z.nc.gz'}
    # if (yr==2000) {file = 'IMOS_OceanCurrent_HV_2000_C-20150521T035421Z.nc.gz'}
    # if (yr==1999) {file = 'IMOS_OceanCurrent_HV_1999_C-20150521T034853Z.nc.gz'}
    # if (yr==1998) {file = 'IMOS_OceanCurrent_HV_1998_C-20150521T034357Z.nc.gz'}
    # if (yr==1997) {file = 'IMOS_OceanCurrent_HV_1997_C-20150521T033755Z.nc.gz'}
    # if (yr==1996) {file = 'IMOS_OceanCurrent_HV_1996_C-20150521T033049Z.nc.gz'}
    # if (yr==1995) {file = 'IMOS_OceanCurrent_HV_1995_C-20150521T032432Z.nc.gz'}
    # if (yr==1994) {file = 'IMOS_OceanCurrent_HV_1994_C-20150521T031623Z.nc.gz'}
    # if (yr==1993) {file = 'IMOS_OceanCurrent_HV_1993_C-20150521T030649Z.nc.gz'}
    # 
    # 
    
    if (yr==1993) {file = 'IMOS_OceanCurrent_HV_DM01_1993_C-20200615T093812Z.nc.gz'}
    if (yr==1994) {file = 'IMOS_OceanCurrent_HV_DM01_1994_C-20200615T110033Z.nc.gz'}
    if (yr==1995) {file = 'IMOS_OceanCurrent_HV_DM01_1995_C-20200615T123704Z.nc.gz'}
    if (yr==1996) {file = 'IMOS_OceanCurrent_HV_DM01_1996_C-20200615T141744Z.nc.gz'}
    if (yr==1997) {file = 'IMOS_OceanCurrent_HV_DM01_1997_C-20200615T154026Z.nc.gz'}
    if (yr==1998) {file = 'IMOS_OceanCurrent_HV_DM01_1998_C-20200615T164237Z.nc.gz'}
    if (yr==1999) {file = 'IMOS_OceanCurrent_HV_DM01_1999_C-20200615T174520Z.nc.gz'}
    if (yr==2000) {file = 'IMOS_OceanCurrent_HV_DM01_2000_C-20200615T184648Z.nc.gz'}
    if (yr==2001) {file = 'IMOS_OceanCurrent_HV_DM01_2001_C-20200615T194911Z.nc.gz'}
    if (yr==2002) {file = 'IMOS_OceanCurrent_HV_DM01_2002_C-20200615T205247Z.nc.gz'}
    if (yr==2003) {file = 'IMOS_OceanCurrent_HV_DM01_2003_C-20200615T215257Z.nc.gz'}
    if (yr==2004) {file = 'IMOS_OceanCurrent_HV_DM01_2004_C-20200615T225422Z.nc.gz'}
    if (yr==2005) {file = 'IMOS_OceanCurrent_HV_DM01_2005_C-20200615T235235Z.nc.gz'}
    if (yr==2006) {file = 'IMOS_OceanCurrent_HV_DM01_2006_C-20200616T005101Z.nc.gz'}
    if (yr==2007) {file = 'IMOS_OceanCurrent_HV_DM01_2007_C-20200616T015143Z.nc.gz'}
    if (yr==2008) {file = 'IMOS_OceanCurrent_HV_DM01_2008_C-20200616T025018Z.nc.gz'}
    if (yr==2009) {file = 'IMOS_OceanCurrent_HV_DM01_2009_C-20200616T034913Z.nc.gz'}
    if (yr==2010) {file = 'IMOS_OceanCurrent_HV_DM01_2010_C-20200616T044857Z.nc.gz'}
    if (yr==2011) {file = 'IMOS_OceanCurrent_HV_DM01_2011_C-20200616T054753Z.nc.gz'}
    if (yr==2012) {file = 'IMOS_OceanCurrent_HV_DM01_2012_C-20200616T064615Z.nc.gz'}
    if (yr==2013) {file = 'IMOS_OceanCurrent_HV_DM01_2013_C-20200616T074441Z.nc.gz'}
    if (yr==2014) {file = 'IMOS_OceanCurrent_HV_DM01_2014_C-20200616T085821Z.nc.gz'}
    if (yr==2015) {file = 'IMOS_OceanCurrent_HV_DM01_2015_C-20200616T101243Z.nc.gz'}
    if (yr==2016) {file = 'IMOS_OceanCurrent_HV_DM01_2016_C-20200616T114334Z.nc.gz'}
    if (yr==2017) {file = 'IMOS_OceanCurrent_HV_DM01_2017_C-20200616T131008Z.nc.gz'}
    if (yr==2018) {file = 'IMOS_OceanCurrent_HV_DM01_2018_C-20200616T144027Z.nc.gz'}
    if (yr==2019) {file = 'IMOS_OceanCurrent_HV_DM01_2019_C-20200616T160523Z.nc.gz'}
    if (yr==2020) {file = 'IMOS_OceanCurrent_HV_DM01_2020_C-20200616T163221Z.nc.gz'}
    
    alt_url <- paste0(url_base,file)
    nc <- nc_open(alt_url, write=FALSE, readunlim=TRUE, verbose=FALSE)
    
    # Approximate nearest neighbour
    idx_lon <- ann(as.matrix(nc$dim$LONGITUDE$vals), as.matrix(dat$Longitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
    idx_lat <- ann(as.matrix(nc$dim$LATITUDE$vals), as.matrix(dat$Latitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
    
    idx_time <- ann(as.matrix(yday(days(nc$dim$TIME$vals) + ymd("1985-01-01"))), as.matrix(yday(dat$Date[i])), k = 1, verbose = FALSE)$knnIndexDist[,1]
    
    cnt <- c(1,1,1)
    if (res_spat > 1) { # If more than 1x1 pixel is requested we adjust the idx by res_spat/2 and count by res_spa
      idx_lon <- idx_lon - floor(res_spat/2)
      idx_lon <- idx_lon - floor(res_spat/2)
      cnt <- c(res_spat, res_spat, 1)
    }
    
    # Extract, average and add Altimetry data to dataframe
    GSLA <- ncvar_get(nc, "GSLA", start=c(idx_lon, idx_lat, idx_time), count = cnt)
    dat$GSLA[i] <- mean(GSLA, na.rm = TRUE)
    GSL <- ncvar_get(nc, "GSL", start=c(idx_lon, idx_lat, idx_time), count = cnt)
    dat$GSL[i] <- mean(GSL, na.rm = TRUE)
    
    UCUR <- ncvar_get(nc, "UCUR", start=c(idx_lon, idx_lat, idx_time), count = cnt)
    dat$UCUR[i] <- mean(UCUR, na.rm = TRUE)
    VCUR <- ncvar_get(nc, "VCUR", start=c(idx_lon, idx_lat, idx_time), count = cnt)
    dat$VCUR[i] <- mean(VCUR, na.rm = TRUE) 
    
    dat$Longitude_Actual[i] <- nc$dim$LONGITUDE$vals[idx_lon]
    dat$Latitude_Actual[i] <- nc$dim$LATITUDE$vals[idx_lat]
    
    setTxtProgressBar(pb, i)
    cat("\n")
    print(i)
    nc_close(nc)
  }
  
  return(dat)
}
