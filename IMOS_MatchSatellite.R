# THis software downloads IMOS satellite data at a vector of 
# given locations and times using the OPeNDAP protocol.
# A nice website with a description is here:
# https://opendap.github.io/documentation/QuickStart.html
# 
library(readr)
library(stringr)
library(dplyr)
library(ncdf4)
library(yaImpute)

filename <- "TestData_IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance_HTL.csv"

# Get the latitude/longitude and date from the file
dat <- read_csv(filename)
dat <- dat[600:610,]

# Decide on the products I want.
# Possible values are: 
# "sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
# "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
# "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490"

pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
        "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
        "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")



# Set resolution
res_temp <-  1 #day  Decide on the temporal resolution. Lets start with raw daily data served by IMOS

res_spat <-  10 #pixels Decide on the spatial resolution. Lets start with raw ~1 km data served by IMOS

# Create a NA matrix of size length(Latitude) x no_products and fill it incrementally
mat <- matrix(data=NA,nrow=length(dat$Latitude),ncol=length(pr))

pb <- txtProgressBar(min = 0, max = length(dat$Latitude), style = 3)
for (i in 1:length(dat$Latitude)) { # Loop through all rows in the data for each variable you want
  
  # Make sure month and day have a leading zero if less than 10
  mth <- str_pad(dat$Month[i],2,"left",pad="0")
  dy <- str_pad(dat$Day[i],2,"left",pad="0")
  
  for (j in 1:length(pr)) { # Loop through all rows in the data for each variable you want
    
    url_base <- "http://rs-data1-mel.csiro.au/thredds/dodsC/imos-srs/oc/aqua/" # Base URL
    url <- paste0(url_base,res_temp,"d/",dat$Year[i],"/",mth,"/A",dat$Year[i],mth,dy,".L2OC_BASE.aust.",pr[j],".nc")
    
    if (dat$SampleDateLocal[i] > as.POSIXct("2002-07-03 00:00:00 AEST")){ # MODIS is only available from 4/7/2002
      nc <- nc_open(url, write=FALSE, readunlim=TRUE, verbose=FALSE)
      
      # Approximate nearest neighbour
      idx_lon <- ann(as.matrix(nc$dim$lon$vals), as.matrix(dat$Longitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
      idx_lat <- ann(as.matrix(nc$dim$lat$vals), as.matrix(dat$Latitude[i]), k = 1, verbose = FALSE)$knnIndexDist[,1]
      cnt <- c(1,1,1)
      if (res_spat > 1) { # If more than 1x1 pixel is requested we adjust the idx by res_spat/2 and count by res_spa
        idx_lon <- idx_lon - floor(res_spat/2)
        idx_lon <- idx_lon - floor(res_spat/2)
        cnt <- c(res_spat, res_spat, 1)
      }
      
      out <- ncvar_get(nc, pr[j], start=c(idx_lon, idx_lat, 1), count = cnt)
      
      mat[i,j] <- mean(out,na.rm = TRUE)
      
    }
  }
  setTxtProgressBar(pb, i)
}

# Now lookp through and assign columns to the variables
for (j in 1:length(pr)) {
  dat <- dat %>% 
    mutate(!!pr[j] := mat[,j])
}

