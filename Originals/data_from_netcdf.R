
########################################
#claire getting ghrsst data, Product 110
########################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)
library(lubridate)

dat <- read.csv("H://Data_exports/cpr_segments.csv")
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)

dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d",tz="GMT"))
dat$file <-   paste("d", format(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d"),format="%Y%m%d"),".nc" , sep="")

dat <- subset(dat, date > as.POSIXct(strptime("2009-02-23",format="%Y-%m-%d",tz="GMT")))

nz <- length(dat$LATITUDE)

for (i in 1:nz) {
  file <- dat$file[i]
  name <- sprintf("//argos-hba.it.csiro.au/sdode-data/Product110/%s", file)
  nc <- nc_open(name)	
  x2 <- (ceiling(dat$LONGITUDE/0.5)*0.5)[i] #max long
  x1 <- (floor(dat$LONGITUDE/0.5)*0.5)[i]   #min long
  y1 <- (floor(dat$LATITUDE/0.5)*0.5)[i]    #min lat
  y2 <- (ceiling(dat$LATITUDE/0.5)*0.5)[i]  #max lat  
  x <- dat$LON[i]
  y <- dat$LAT[i]
  LonIdxc <- which( nc$dim$lon$vals == x2)
  LonIdxf <- which( nc$dim$lon$vals == x1)
  LatIdxc <- which( nc$dim$lat$vals == y1)
  LatIdxf <- which( nc$dim$lat$vals == y2)

  meanicc <- ncvar_get(nc, "analysed_sst", start=c(LonIdxc, LatIdxc, 1), count = c(1,1,1))
  meaniff <- ncvar_get(nc, "analysed_sst", start=c(LonIdxf, LatIdxf, 1), count = c(1,1,1))
  meanifc <- ncvar_get(nc, "analysed_sst", start=c(LonIdxf, LatIdxc, 1), count = c(1,1,1))
  meanicf <- ncvar_get(nc, "analysed_sst", start=c(LonIdxc, LatIdxf, 1), count = c(1,1,1))

  v1 <- ((x2-x)/(x2-x1))*meanifc + ((x-x1)/(x2-x1)*meanicc)
  v2 <- ((x2-x)/(x2-x1))*meaniff + ((x-x1)/(x2-x1)*meanicf)
  vxy <- ((y2-y)/(y2-y1))*v1 + ((y-y1)/(y2-y1)*v2)
  
  dat[i,8] <- as.numeric(vxy-273.15)
  
}

write.csv(dat, "H://Data_exports/cpr_segments.csv")


########################################
#claire getting noaa reynolds weekly data, Product 120
########################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)
library(lubridate)

dat <- read.csv("~/Data_exports/VJD_records_extract_20180802_qld.csv")
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)

dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d %H:%M:%S",tz="GMT"))
dat$file <- paste("sst.wkmean.1981-1989.nc")
dat$file[format(dat$date, "%Y") > 1989] <- paste("sst.wkmean.1990-present.nc")

#nc <-nc_open("//argos-hba.it.csiro.au/sdode-data/Product120/sst.wkmean.1990-present.nc")

ref_date <- dmy("01-01-1800")
dat$dd <- round(as.numeric(difftime(dat$date,ref_date)),0)
summary(dat$dd)
dat <- subset(dat, dat$dd > 66410)

nz <- length(dat$LATITUDE)
name <- sprintf("//argos-hba.it.csiro.au/sdode-data/Product120/%s", file)
nc <- nc_open(name)

for (i in 1:nz) {
  file <- dat$file[i]
  x2 <- min(which(nc$dim$lon$vals >= dat$LONGITUDE[i])) #max long
  x1 <- max(which(nc$dim$lon$vals <= dat$LONGITUDE[i]))  #min long
  y1 <- min(which(nc$dim$lat$vals <= dat$LATITUDE[i]))   #min lat
  y2 <- max(which(nc$dim$lat$vals >= dat$LATITUDE[i]))  #max lat  
    ttt <- dat$dd[i]
  tt <- (round((dat$dd[i]-66410)/7)*7)+66410
  t <- which(nc$dim$time$vals > tt)
  tIdx <- min(t)
  
  vxy <- ncvar_get(nc, "sst", start=c(x1, y2, tIdx), count = c(1,1,1))
  
  dat[i,9] <- as.numeric(vxy)
  
}

write.csv(dat, "H://Data_exports/cpr_segments.csv")



########################################
#claire getting MODIS chla 4km 8 day, Product 416
########################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)
library(lubridate)

dat <- read.csv("H://data_exports/cpr_segments.csv")
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)
dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d",tz="GMT"))
#dat <- dat[,c(1,8:12,22:23)]
dat$day <- format(dat$date, "%j")
dat$YEAR <- format(dat$date, "%Y")
dat$day[dat$day > 361] <- 361
dat$day <- as.numeric(dat$day)
dat$daye <- round(dat$day/8,0)*8
dat$daye[dat$daye == 0] <- 8

dat$day2 <- dat$daye - 7
dat$day2[dat$day2 > 362] <- 365
dat$day2[dat$YEAR == 2004 & dat$day > 360] <- 366
dat$day2[dat$YEAR == 2008 & dat$day > 360] <- 366
dat$day2[dat$YEAR == 2012 & dat$day > 360] <- 366
dat$day2[dat$YEAR == 2016 & dat$day > 360] <- 366

dat$daye <- sprintf("%03d", dat$daye)
dat$days <- sprintf("%03d", dat$day2)

dat$file <- paste0("A",dat$YEAR,dat$days,dat$YEAR,dat$daye,".L3m_8D_CHL_chlor_a_4km.nc")
#nc <-nc_open("//argos-hba.it.csiro.au/sdode-data/Product416/A20021772002184.L3m_8D_CHL_chlor_a_4km.nc")

nz <- length(dat$LATITUDE)

for (i in 1:nz) {
  file <- dat$file[i]
  name <- sprintf("//argos-hba.it.csiro.au/sdode-data/Product416/%s", file)
  nc <- nc_open(name)	
  x2 <- min(which(nc$dim$lon$vals >= dat$LONGITUDE[i])) #max long
  x1 <- max(which(nc$dim$lon$vals <= dat$LONGITUDE[i]))  #min long
  y1 <- min(which(nc$dim$lat$vals <= dat$LATITUDE[i]))   #min lat
  y2 <- max(which(nc$dim$lat$vals >= dat$LATITUDE[i]))  #max lat  
  vxy <- ncvar_get(nc, "chlor_a", start=c(x1, y2), count = c(1,1))
  
  dat[i,14] <- as.numeric(vxy)
  
}

write.csv(dat, "cpr_segments.csv")

########################################
#claire getting MODIS chla 4km monthly, Product 417 
########################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)
library(lubridate)

dat <- read.csv("~/Data_exports/VJD_records_extract_20180802_qld.csv")
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)
dat <- subset(dat, (dat$YEAR > 2002))
dat <- dat[,c(1,8:13,22:23)]
dat$date <-   as.POSIXct(strptime(dat$EVENT_DATE,format="%Y-%m-%d %H:%M:%S",tz="GMT"))

dat$datm <- floor_date(dat$date, "month")
dat$doys <- format(dat$datm, "%j")
dat$datme <- ceiling_date(dat$date, "month")-ddays(1)
dat$doye <- format(dat$datme, "%j")
dat$doye[dat$CSIRO_ID==21243] <- "121"
dat$doye[dat$CSIRO_ID==20206] <- "060"
dat$doye[dat$doys=="001"] <- "031"

dat$file <- paste0("A",dat$YEAR,dat$doys,dat$YEAR,dat$doye,".L3m_MO_CHL_chlor_a_4km.nc")
#nc <-nc_open("//argos-hba.it.csiro.au/sdode-data/Product416/A20021772002184.L3m_8D_CHL_chlor_a_4km.nc")

nz <- length(dat$LATITUDE)

for (i in 1:nz) {
  file <- dat$file[i]
  name <- sprintf("//argos-hba.it.csiro.au/sdode-data/Product417/%s", file)
  nc <- nc_open(name)	
  x2 <- min(which(nc$dim$lon$vals >= dat$LONGITUDE[i])) #max long
  x1 <- max(which(nc$dim$lon$vals <= dat$LONGITUDE[i]))  #min long
  y1 <- min(which(nc$dim$lat$vals <= dat$LATITUDE[i]))   #min lat
  y2 <- max(which(nc$dim$lat$vals >= dat$LATITUDE[i]))  #max lat  
  vxy <- ncvar_get(nc, "chlor_a", start=c(x1, y2), count = c(1,1))
  
  dat[i,19] <- as.numeric(vxy)
  
}

write.csv(dat, "H://Data_exports/cpr_segments.csv")


########################################
#claire getting bat ga data, Product 40
########################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)

dat <- read.csv("H://Data_exports/cpr_segments.csv")
dat <- trim
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)

dat$date <-   as.POSIXct(strptime(dat$dates,format="%Y-%m-%d",tz="GMT"))
dat$file <-   paste("d", format(strptime(dat$date,format="%Y-%m-%d"),format="%Y%m%d"),".nc" , sep="")

nc <- nc_open("//argos-hba.it.csiro.au/sdode-data/Product40/ga_9sec_bath.nc")
nc 

nz <- length(dat$Lat)

for (i in 1:nz) {
  x2 <- (ceiling(dat$Lon/0.25)*0.25)[i] #max long
  x1 <- (floor(dat$Lon/0.25)*0.25)[i]   #min long
  y1 <- (floor(dat$Lat/0.25)*0.25)[i]    #min lat
  y2 <- (ceiling(dat$Lat/0.25)*0.25)[i]  #max lat  
  x <- dat$Lon[i]
  y <- dat$Lat[i]
  LonIdxc <- which( nc$dim$lon$vals == x2)
  LonIdxf <- which( nc$dim$lon$vals == x1)
  LatIdxc <- which( nc$dim$lat$vals == y1)
  LatIdxf <- which( nc$dim$lat$vals == y2)
  
  meanicc <- ncvar_get(nc, "height", start=c(LonIdxc, LatIdxc), count = c(1,1))
  meaniff <- ncvar_get(nc, "height", start=c(LonIdxf, LatIdxf), count = c(1,1))
  meanifc <- ncvar_get(nc, "height", start=c(LonIdxf, LatIdxc), count = c(1,1))
  meanicf <- ncvar_get(nc, "height", start=c(LonIdxc, LatIdxf), count = c(1,1))
  
  v1 <- ((x2-x)/(x2-x1))*meanifc + ((x-x1)/(x2-x1)*meanicc)
  v2 <- ((x2-x)/(x2-x1))*meaniff + ((x-x1)/(x2-x1)*meanicf)
  vxy <- ((y2-y)/(y2-y1))*v1 + ((y-y1)/(y2-y1)*v2)
  
  trim[i,30] <- as.numeric(vxy)
  
}

write.csv(dat, "H://Data_exports/cpr_segments.csv")



##################################
#TURBIDITY LAYER FROM NETCDF Product 407
##################################
rm(list = ls())
setwd("H://Data_exports")

library(ncdf4)
library(chron)
library(lubridate)

dat <- read.csv("cpr_segments.csv")
#opts = curlOptions(maxconnects = 3, netrc = TRUE) # Set options (fixes error 421 - too many connections)
nc <- nc_open("//argos-hba.it.csiro.au/sdode-data/Product407/A20170012017031.L3m_MO_KD490_Kd_490_4km.nc")
nc # Calls the netCDF, which lists all of the variables

dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d %H:%M:%S",tz="GMT"))
dat$YEAR <- format(dat$date, "%Y")
dat$datm <- floor_date(dat$date, "month")
dat$doys <- format(dat$datm, "%j")
dat$datme <- ceiling_date(dat$date, "month")-ddays(1)
dat$doye <- format(dat$datme, "%j")
dat$doye[dat$CSIRO_ID==21243] <- "121"
dat$doye[dat$CSIRO_ID==20206] <- "060"
dat$doye[dat$doys=="001"] <- "031"

dat$file <- paste0("A",dat$YEAR,dat$doys,dat$YEAR,dat$doye,".L3m_MO_KD490_Kd_490_4km.nc")

dat <- subset(dat, date > as.POSIXct(strptime("2018-01-01",format="%Y-%m-%d",tz="GMT")))
nz <- length(dat$LATITUDE)

for (i in 1:nz) {
  file <- dat$file[i]
  name <- sprintf("//argos-hba.it.csiro.au/sdode-data/Product407/%s", file)
  nc <- nc_open(name)	
  x2 <- (dat$LONGITUDE+0.025)[i] #max long
  x1 <- (dat$LONGITUDE-0.025)[i]   #min long
  y1 <- (dat$LATITUDE-0.025)[i]    #min lat
  y2 <- (dat$LATITUDE+0.025)[i]  #max lat  
  x <- dat$LONGITUDE[i]
  y <- dat$LATITUDE[i]
  LonIdx <- which( nc$dim$lon$vals > x1 & nc$dim$lon$vals < x2)
  LatIdx <- which( nc$dim$lat$vals > y1 & nc$dim$lat$vals < y2)
  LatIdx <- min(LatIdx)
  LonIdx <- min(LonIdx)
  mean <- ncvar_get(nc, "Kd_490", start=c(LonIdx, LatIdx), count = c(1,1))
  
  dat[i,20] <- as.numeric(mean)
  
}

write.csv(dat, "H://Data_exports/cpr_segments.csv")


##################################
#DUST INDEX FROM NETCDF giovanni nasa
##################################
rm(list = ls())
setwd("~/Data_exports")

library(ncdf4)
library(chron)
library(RCurl)
library(lubridate)

dat <- read.csv("nrs_segments.csv")

dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d %H:%M:%S",tz="GMT"))
dat$YEAR <- format(dat$date, "%Y")
dat$datm <- floor_date(dat$date, "month")
dat$doys <- format(dat$datm, "%j")
dat$datme <- ceiling_date(dat$date, "month")-ddays(1)
dat$doye <- format(dat$datme, "%j")
dat$doye[dat$CSIRO_ID==21243] <- "121"
dat$doye[dat$CSIRO_ID==20206] <- "060"
dat$doye[dat$doys=="001"] <- "031"

dat$file <- paste0("A",dat$YEAR,dat$doys,dat$YEAR,dat$doye,".L3m_MO_RRS_angstrom_9km.nc")

#get the correct files downloaded from ocean colour, only download those that we haven't already got
lev <- length(dat$doys)

for (i in 1:lev){
  year <- dat$YEAR[i]
  doy1 <- dat$doys[i]
  doy2 <- dat$doye[i]
  fname <- paste0("A", year, doy1, year, doy2,".L3m_MO_RRS_angstrom_9km.nc" , sep="")
  file <-  paste0("https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/", year,"/", doy1,"/A" ,year, doy1,year,doy2,".L3m_MO_RRS_angstrom_9km.nc.nc4" , sep="")
  download.file(file, fname, method = "auto",
                quiet = FALSE, mode="wb", cacheOK = TRUE)
}

nz <- length(dat$LATITUDE)

for (i in 1:nz) {
  file <- dat$file[i]
  name <- sprintf("H://Data_exports/%s", file)
  nc <- nc_open(name)	
  x2 <- (dat$LONGITUDE+0.05)[i] #max long
  x1 <- (dat$LONGITUDE-0.05)[i]   #min long
  y1 <- (dat$LATITUDE-0.05)[i]    #min lat
  y2 <- (dat$LATITUDE+0.05)[i]  #max lat  
  x <- dat$LONGITUDE[i]
  y <- dat$LATITUDE[i]
  LonIdx <- which( nc$dim$lon$vals > x1 & nc$dim$lon$vals < x2)
  LatIdx <- which( nc$dim$lat$vals > y1 & nc$dim$lat$vals < y2)
  LatIdx <- min(LatIdx)
  LonIdx <- min(LonIdx)
  v1 <- nc$var[[1]]
  mean <- ncvar_get(nc, v1 , start=c(LonIdx, LatIdx), count = c(1,1))
  
  dat[i,12] <- as.numeric(mean)

}

write.csv(dat, "H://Data_exports/nrs_segments.csv")

##################################
#PAR INDEX FROM NETCDF 
##################################
rm(list = ls())
setwd("~/R/data")

library(ncdf4)
library(chron)
library(RCurl)
library(lubridate)

#get the correct files downloaded from ocean colour
setwd("~/R/data")
dates <- read.csv("cpr_dates.csv")

lev <- length(dat$doys)

for (i in 1:lev){
  year <- dat$YEAR[i]
  doy1 <- dat$doys[i]
  doy2 <- dat$doye[i]
  fname <- paste0("A", year, doy1, year, doy2,".L3m_MO_RRS_angstrom_9km.nc" , sep="")
  file <-  paste0("https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/", year,"/", doy1,"/A" ,year, doy1,year,doy2,".L3m_MO_RRS_angstrom_9km.nc.nc4" , sep="")
  download.file(file, fname, method = "auto",
                quiet = FALSE, mode="wb", cacheOK = TRUE)
}

download.file("https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2007/335/A20073352007365.L3m_MO_RRS_angstrom_9km.nc.nc4", "test.nc4", method = "auto",
                    quiet = FALSE, mode="wb", cacheOK = TRUE)
nc <- nc_open("test.nc4")

## get data from downloaded file
setwd("~/R/data")
dat <- read.csv("~/Data_exports/cpr_segments.csv")
dat$date <-   as.POSIXct(strptime(dat$SAMPLE_DATE,format="%Y-%m-%d %H:%M:%S",tz="GMT"))
dat$doy1 <- format(as.POSIXct(strptime(dat$DOY1,format="%Y-%m-%d")),"%j")
dat$doy2 <- format(as.POSIXct(strptime(dat$DOY2,format="%Y-%m-%d")),"%j")

nz <- length(dat$LATITUDE)
val_id <- matrix(0, nrow=nz, ncol=3)
colnames(val_id)=c("silk_id", "segment_id", "Ang")


for (i in 1:nz) {
  year <- paste(format(strptime(dat$date[i],format="%Y-%m-%d %H:%M:%S"),format="%Y"))
  doy1 <- dat$doy1[i]
  doy2 <- dat$doy2[i]
  
  fname <- paste0("~/R/data/A", year, doy1, year, doy2,".L3m_MO_RRS_angstrom_9km.nc" , sep="")
  
  nc <- nc_open(fname)
  
  x2 <- (dat$LONGITUDE+0.05)[i] #max long
  x1 <- (dat$LONGITUDE-0.05)[i]   #min long
  y1 <- (dat$LATITUDE-0.05)[i]    #min lat
  y2 <- (dat$LATITUDE+0.05)[i]  #max lat  
  x <- dat$LONGITUDE[i]
  y <- dat$LATITUDE[i]
  LonIdx <- which( nc$dim$lon$vals > x1 & nc$dim$lon$vals < x2)
  LatIdx <- which( nc$dim$lat$vals > y1 & nc$dim$lat$vals < y2)
  LatIdx <- min(LatIdx)
  LonIdx <- min(LonIdx)
  v1 <- nc$var[[1]]
  mean <- ncvar_get(nc, v1 , start=c(LonIdx, LatIdx), count = c(1,1))
  
  val_id[i] <- as.numeric(dat$SILK_ID[i])
  val_id[i,2] <- as.numeric(dat$SEGMENT_NO[i])
  val_id[i,3] <- as.numeric(mean)
  nc_close(nc)
  
}

write.csv(val_id, "~/Data_exports/cpr_ang.csv")

#######################################
##  files from ocean colour - never got this bit working correctly
#######################################
library(sp)
library(raster)
library(ncdf4)
#install.packages(c("curl", "httr"))
library(devtools)
#install_github("btupper/threddscrawler")
#install_github("BigelowLab/obpgcrawler")
install_github("btupper/spnc")
withr::with_libpaths(.Library, install_github("btupper/spnc"))
library(obpgcrawler)
library(curl)
library(spnc)
pkg_path <- '/path/to/spnc'
install(pkg_path)

query <- obpg_query(top = 'https://oceandata.sci.gsfc.nasa.gov/opendap/catalog.xml',
                    platform = 'MODISA', 
                    product = 'L3SMI',
                    what = 'most_recent',
                    greplargs = list(pattern='L3m_8D_RRS_angstrom_9km', fixed = TRUE))
query

query <- obpg_query(year = c('2014', '2015'),day = 1:30, greplargs = list(pattern='8D_CHL_chlor_a_4km', fixed = TRUE))
chl <- nc_open(query[['A20140012014008.L3m_8D_CHL_chlor_a_4km.nc']]$url)
chl

list(dates$doy1)
