## IMOS LTS plotting from gridded timeseries products 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Oct 2020
## Updated: 
## 8 Oct 2020 (Written to Git)


# if on a windows machind you need this version of netcdf
#install.packages("devtools")
#devtools::install_github("mdsumner/ncdf4")

suppressPackageStartupMessages({
  library(stats)
  library(ncdf4)
  library(ggplot2)
  library(colorRamps) # for Matlab like colour scheme
  library(data.table)
  library(MBA) # Does bilinear interpolation
  library(reshape2)
  library(lubridate)
  library(tidyverse) ## should go near last to put on tip of search path
})
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_bw() + theme(legend.position = "bottom"))
options(stringsAsFactors = FALSE)   ### the default from R 4.0.0 on

### Helper functions
### 

source("../Plankton/IMOS_Plankton_functions.R")

# access gridded product files from IMOS thredds server. If you want to only use that file you can save it to a location and change the address below
# if these break then the file has probably been updated and the name changed - annoying .....

PHB <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NSW/PH100/gridded_timeseries/IMOS_ANMN-NSW_TZ_20091029_PH100_FV02_TEMP-gridded-timeseries_END-20190828_C-20200108.nc")
MAI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSMAI/gridded_timeseries/IMOS_ANMN-NRS_TZ_20080411_NRSMAI_FV02_TEMP-gridded-timeseries_END-20200522_C-20201007.nc")
NSI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSNSI/gridded_timeseries/IMOS_ANMN-NRS_TZ_20101213_NRSNSI_FV02_TEMP-gridded-timeseries_END-20191214_C-20201007.nc")
ROT <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSROT/gridded_timeseries/IMOS_ANMN-NRS_TZ_20081120_NRSROT_FV02_TEMP-gridded-timeseries_END-20200703_C-20201007.nc")
KAI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSKAI/gridded_timeseries/IMOS_ANMN-NRS_TZ_20080212_NRSKAI_FV02_TEMP-gridded-timeseries_END-20191127_C-20200526.nc")

# create a list with the filenames and station codes
stations <- list("PHB" = PHB, "MAI" = MAI, "NSI" = NSI, "KAI" = KAI, "ROT" = ROT)
names <- c("PHB", "MAI", "NSI", "KAI", "ROT")
stations <- mapply(c, stations, names, SIMPLIFY = FALSE) 

#############################################################################################################################################
## to plot mooring data as contour plot

mooringplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f)) # open each element of list in turn
  nc <- eval(nc)

  temp <- ncvar_get(nc, "TEMP") # access variables from netcdf to plot
  depth <- ncvar_get(nc, "DEPTH")
  time <- ncvar_get(nc, "TIME")
  
  rownames(temp) <- depth
  colnames(temp) <- time
  
  nc_close(nc)
  
  longdata <- melt(temp) %>% arrange(Var2, Var1) %>% drop_na()
  colnames(longdata) <- c("Depth", "Time", "Temp")
  
  refDate <- ymd_hms("1950-01-01 00:00:00") # netcdf time in days since 1950
  
  longdata$Date <- as.POSIXct(longdata$Time * 3600 * 24, origin = refDate)
  longdata$DOY <- as.numeric(strftime(longdata$Date, format = "%j"))
  
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  longdata_sum <- longdata %>% # organise data for plotting
    mutate(Date = decimal_date(Date)) %>%
    select(Date, Depth, Temp) %>% arrange(Date, Depth) %>% 
    drop_na(Temp)
  
  mba <- mba.surf(longdata_sum, no.X = 200, no.Y = 200, n = 1, m = 5)
  # This is just to organise dataframe for plotting
  dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
  df <- melt(mba$xyz.est$z, varnames = c('Time', 'Depth'), value.name = 'Temperature')
  
  tpp <- ggplot(data = df, aes(Time, Depth)) +
    geom_raster(aes(fill = Temperature), interpolate = F) +
    scale_fill_gradientn(colours = matlab.like(7), na.value = "white") +
    scale_x_continuous(breaks = seq(2010,2020,2), labels = seq(2010,2020,2), expand = c(0, 0)) +
    theme_bw(base_size = 16) +
    xlab("Time") +
    ylab("Depth (m)") +
    labs(fill = expression("Temperature ("* degree *"C)")) +
    theme(legend.position = "bottom") +
    scale_y_reverse(expand = c(0, 0))  # tight axis scale
  
  ggsave(paste0("LTS_plots/gridded_mooringsT", NRS, ".png", sep=''), tpp, width = 12, height = 6, dpi = 600)
})

##################################################################################################################################################
# Plots mooring data as a time series 

timeseriesplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f))
  nc <- eval(nc)

  temp <- ncvar_get(nc, "TEMP")
  depth <- ncvar_get(nc, "DEPTH")
  time <- ncvar_get(nc, "TIME")
  
  rownames(temp) <- depth
  colnames(temp) <- time
  
  nc_close(nc)
  
  refDate <- ymd_hms("1950-01-01 00:00:00")
  timeNow <- ymd_hms("2019-08-28 00:00:00")
  
  longdata <- reshape2::melt(temp) %>% arrange(Var2, Var1) %>% drop_na()
  colnames(longdata) <- c("Depth", "Time", "Temp")
  
  longdata <- longdata %>% mutate(DepthBin = cut(Depth, breaks = c(0, 20, 40, 60, 80, 100)),
                                  DepthBin = recode(DepthBin, "(0,20]" = "0 - 20m",
                                                    "(20,40]" = "20 - 40m",
                                                    "(40,60]" = "40 - 60m",
                                                    "(60,80]" = "60 - 80m",
                                                    "(80,100]" = "80 - 1000m"),
                                  Date = as.POSIXct(Time * 3600 * 24, origin = refDate),
                                  Dated = as.Date(Date, format = "%y%m%d"),
                                  DOY = as.numeric(strftime(Date, format = "%j"))) %>% droplevels()   %>%
    drop_na() %>% 
    group_by(DepthBin, Time, Date, Dated, DOY) %>% summarise(Temp = mean(Temp)) %>% untibble()

  P <- ggplot() + geom_line(data = longdata, aes(Date, Temp)) + 
    facet_grid(DepthBin ~ .) + labs(y = expression("Temperature ("* degree *"C)")) +
    scale_x_datetime(breaks = seq(as.POSIXct("2010-01-01 00:00:00"),
                                             as.POSIXct("2020-01-01 00:00:00"), "2 years"), 
                                  labels = c(2010,2012,2014,2016,2018,2020), expand = c(0,0))
    
  
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  ggsave(paste("LTS_plots/gridded_timeseries_",NRS,".png", sep = ""), P, width = 6, height = 8,  dpi=600)
  
})

##############################################################################################################################################################
# plot a climatology at each station

climatologyplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f))
  nc <- eval(nc)
  
  temp <- ncvar_get(nc, "TEMP")
  depth <- ncvar_get(nc, "DEPTH")
  time <- ncvar_get(nc, "TIME")
  
  rownames(temp) <- depth
  colnames(temp) <- time
  
  nc_close(nc)
  
  refDate <- ymd_hms("1950-01-01 00:00:00")
  timeNow <- ymd_hms("2019-08-28 00:00:00")
  
  longdata <- reshape2::melt(temp) %>% arrange(Var2, Var1) %>% drop_na()
  colnames(longdata) <- c("Depth", "Time", "Temp")
  
  longdata <- longdata %>% mutate(DepthBin = cut(Depth, breaks = c(0, 20, 40, 60, 80, 100)),
                                  DepthBin = recode(DepthBin, "(0,20]" = "0 - 20m",
                                                    "(20,40]" = "20 - 40m",
                                                    "(40,60]" = "40 - 60m",
                                                    "(60,80]" = "60 - 80m",
                                                    "(80,100]" = "80 - 1000m"),
                                  Date = as.POSIXct(Time * 3600 * 24, origin = refDate),
                                  Dated = as.Date(Date, format = "%y%m%d"),
                                  DOY = as.numeric(strftime(Date, format = "%j"))) %>% droplevels()   %>%
    drop_na() %>% 
    group_by(DepthBin, DOY) %>% summarise(Temp = mean(Temp)) %>% untibble()
  
  P <- ggplot() + geom_smooth(data = longdata, aes(DOY, Temp)) + 
    facet_grid(DepthBin ~ ., scales = "free") + labs(y = expression("Temperature ("* degree *"C)"),
                                    x = "Day of Year") +
    scale_x_continuous(expand = c(0,0))
  
  
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  ggsave(paste("LTS_plots/gridded_climatology_",NRS,".png", sep = ""), P, width = 6, height = 8,  dpi=600)
  
})
