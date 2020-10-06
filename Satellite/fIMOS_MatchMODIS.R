# This software downloads IMOS satellite data at a vector of 
# given locations and times using the OPeNDAP protocol.
# A nice website with a description is here:
# https://opendap.github.io/documentation/QuickStart.html
#
# Inputs: 
# dat is a dataframe or tibble with the following columns:
# Latitude
# Longitude
# and either
# Date (as POSIX)
# or
# Day &
# Month &
# Year
# pr are the MODIS products you wish to extract. See below for a list of options
#
# Optional Inputs:
# res_spat - Spatial resolution. How many pixels (n x n) to download in each direction
# res_temp - What temporal averaging: 1 day (1d), 1 month(1m), 1 year (1y),...
# Monthly climatology (1mNy), Annual climatology (12mNy)

# Possible products to download are: 
# "sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
# "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
# "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oci", "chl_oc3", "chl_gsm", "K_490"

# Options for res_temp
# "1d", "1m", "1y", "1mNy", "12mNy"

fIMOS_MatchMODIS <- function(dat, pr, ...) {
  library(stringr) #install.packages("stringr")  
  library(dplyr) #install.packages("dplyr")
  library(ncdf4) #install.packages("ncdf4")
  library(yaImpute) #install.packages("yaImpute")
  library(lubridate) #install.packages("lubridate")
  
  # Set resolution
  if (!exists("res_temp")){
    print("Defaulting to daily satellite data")
    res_temp <-  "1d"
  }
  
  if (!exists("res_spat")){
    print("Defaulting to 1 pixel x 1 pixel. Provide res_spat if you want to increase")
    res_spat <-  1
  }
  
  if (("Date" %in% colnames(dat))==FALSE) {
    stop("No Date column found in data. Please include a Date column in POSIXct format")
  }
  
  # If Day, Month, Year doesn't exist we create them
  if (sum(c("Day","Month","Year") %in% colnames(dat)) != 3) {   
    dat <- dat %>% 
      mutate(Day = day(Date),
             Month = month(Date),
             Year = year(Date))
  }
  
  # Create a NA matrix of size length(Latitude) x no_products and fill it incrementally
  mat <- matrix(data=NA,nrow=length(dat$Latitude),ncol=length(pr))
  
  pb <- txtProgressBar(min = 0, max = length(dat$Latitude), style = 3)
  for (i in 1:length(dat$Latitude)) { # Loop through all rows in the data for each variable you want
    
    # Make sure month and day have a leading zero if less than 10
    mth <- str_pad(dat$Month[i],2,"left",pad="0")
    dy <- str_pad(dat$Day[i],2,"left",pad="0")
    
    for (j in 1:length(pr)) { # Loop through all rows in the data for each variable you want
      url_base <- paste0("http://rs-data1-mel.csiro.au/thredds/dodsC/imos-srs/oc/aqua/",res_temp,"/") # Base URL
      
      if (str_detect(res_temp,"1d")){
        imos_url <- paste0(url_base, dat$Year[i],"/",mth,"/A",dat$Year[i],mth,dy,".L2OC_BASE.aust.",pr[j],".nc")
        vr <- pr[j]
      }
      if (str_detect(res_temp,"1m")){
        imos_url <- paste0(url_base, dat$Year[i],"/",dat$Year[i],mth,".",pr[j],".nc")
        vr <- paste0(pr[j],"_mean")
      }
      if (str_detect(res_temp,"1y")){
        imos_url <- paste0(url_base,dat$Year[i],".",pr[j],"_mean.nc4")
        vr <- paste0(pr[j],"_mean_mean")
      }
      if (str_detect(res_temp, "1mNy")){
        imos_url <- paste0(url_base,"2003-2014.",mth,".",pr[j],"_mean.nc4")
        vr <- paste0(pr[j],"_mean_mean")
      }
      if (str_detect(res_temp, "12mNy")){
        imos_url <- paste0(url_base,"2003-2014x01-12.",pr[j],"_mean_mean.nc4")
        vr <- paste0(pr[j],"_mean_mean")
      }
      
      
      tryCatch({ # Not all dates will exist
        nc <- nc_open(imos_url, write=FALSE, readunlim=TRUE, verbose=FALSE)
        # Approximate nearest neighbour
        idx_lon <- ann(as.matrix(nc$dim$lon$vals), as.matrix(dat$Longitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
        idx_lat <- ann(as.matrix(nc$dim$lat$vals), as.matrix(dat$Latitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
        cnt <- c(1,1,1)
        if (res_spat > 1) { # If more than 1x1 pixel is requested we adjust the idx by res_spat/2 and count by res_spa
          idx_lon <- idx_lon - floor(res_spat/2)
          idx_lat <- idx_lat - floor(res_spat/2)
          cnt <- c(res_spat, res_spat, 1)
        }
        out <- ncvar_get(nc, vr, start=c(idx_lon, idx_lat, 1), count = cnt)
        
        mat[i,j] <- mean(out, na.rm = TRUE)
      }, 
      error = function(cond) {
        print(paste0('Date (',dat$Day[i],'/',dat$Month[i],'/',dat$Year[i],') not available.'))
        }
      )
    }
    setTxtProgressBar(pb, i)
  }
  
  # Now lookp through and assign columns to the variables
  for (j in 1:length(pr)) {
    nm <- paste0(pr[j],"_",res_temp)
    dat <- dat %>% 
      mutate(!!nm := mat[,j])
  }
  return(dat)
}
