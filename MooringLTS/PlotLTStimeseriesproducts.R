## IMOS LTS plotting from hourly timeseries products 
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

# access hourly time series product files from IMOS thredds server. 
# If the web links don't work (windows quirks) you may want to download and save it, change the address below to that location
# if these break then the file has probably been updated and the name changed - annoying .....
# no similar file available for PH100 / PHB checking this???

#PHB <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NSW/PH100/hourly_timeseries/IMOS_ANMN-NSW_TZ_20091029_PH100_FV01_TEMP-aggregated-timeseries_END-20190612_C-20190819.nc")
MAI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSMAI/hourly_timeseries/IMOS_ANMN-NRS_BOSTUZ_20080411_NRSMAI_FV02_hourly-timeseries_END-20200522_C-20201007.nc")
NSI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSMAI/hourly_timeseries/IMOS_ANMN-NRS_BOSTZ_20101213_NRSNSI_FV02_hourly-timeseries_END-20191214_C-20201007.nc")
ROT <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSMAI/hourly_timeseries/IMOS_ANMN-NRS_STZ_20081120_NRSROT_FV02_hourly-timeseries_END-20200703_C-20201007.nc")
KAI <- paste0("http://thredds.aodn.org.au/thredds/dodsC/IMOS/ANMN/NRS/NRSKAI/hourly_timeseries/IMOS_ANMN-NRS_BOSTZ_20080212_NRSKAI_FV02_hourly-timeseries_END-20200527_C-20201007.nc")

# create a list with the filenames and station codes
stations <- list(#"PHB" = PHB, #add back in if PHB or PH100 product
                 "MAI" = MAI, "NSI" = NSI, "KAI" = KAI, "ROT" = ROT)
names <- c(#"PHB", #add back in if PHB or PH100 product
           "MAI", "NSI", "KAI", "ROT")
stations <- mapply(c, stations, names, SIMPLIFY = FALSE) 

#############################################################################################################################################
## to plot mooring data as contour plot

mooringplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f)) # open each element of list in turn
  nc <- eval(nc)

  refDate <- ymd_hms("1950-01-01 00:00:00") # netcdf time in days since 1950
 
  df <- data.frame("Time" = ncvar_get(nc, "TIME"),
                   "Depth" = ncvar_get(nc, "DEPTH"),
                   "Temp" = ncvar_get(nc, "TEMP")) %>%
    mutate(Date = as.POSIXct(Time * 3600 * 24, origin = refDate),
           DOY = as.numeric(strftime(Date, format = "%j")))
  
  nc_close(nc)
   
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  df_sum <- df %>% # organise data for plotting
    mutate(Date = decimal_date(Date)) %>%
    select(Date, Depth, Temp) %>% arrange(Date, Depth) %>% 
    drop_na()
  
  mba <- mba.surf(df_sum, no.X = 200, no.Y = 200, n = 1, m = 5)
  # This is just to organise dataframe for plotting
  dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
  df_sum <- melt(mba$xyz.est$z, varnames = c('Time', 'Depth'), value.name = 'Temperature')
  
  tpp <- ggplot(data = df_sum, aes(Time, Depth)) +
    geom_raster(aes(fill = Temperature), interpolate = F) +
    scale_fill_gradientn(colours = matlab.like(7), na.value = "white") +
    scale_x_continuous(breaks = seq(2010,2020,2), labels = seq(2010,2020,2), expand = c(0, 0)) +
    theme_bw(base_size = 12) +
    xlab("Time") +
    ylab("Depth (m)") +
    labs(fill = expression("Temperature ("* degree *"C)")) +
    theme(legend.position = "bottom") +
    scale_y_reverse(expand = c(0, 0))  # tight axis scale
  x11(width = 12, height = 6)
  ggsave(paste0("LTS_plots/hts_mooringsT", NRS, ".png", sep=''), tpp, width = 12, height = 6, dpi = 600)
})

##################################################################################################################################################
# Plots mooring data as a time series 

timeseriesplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f))
  nc <- eval(nc)

  refDate <- ymd_hms("1950-01-01 00:00:00") # netcdf time in days since 1950
  
  df <- data.frame("Time" = ncvar_get(nc, "TIME"),
                   "Depth" = ncvar_get(nc, "DEPTH"),
                   "Temp" = ncvar_get(nc, "TEMP")) %>%
    mutate(Date = as.POSIXct(Time * 3600 * 24, origin = refDate),
           DOY = as.numeric(strftime(Date, format = "%j")))
  
  nc_close(nc)
  
  df_sum <- df %>% mutate(DepthBin = cut(Depth, breaks = c(0, 20, 40, 60, 80, 100)), # separate into depth bins
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

  P <- ggplot() + geom_line(data = df_sum, aes(Date, Temp)) + 
    facet_grid(DepthBin ~ .) + labs(y = expression("Temperature ("* degree *"C)")) +
    scale_x_datetime(breaks = seq(as.POSIXct("2010-01-01 00:00:00"),
                                             as.POSIXct("2020-01-01 00:00:00"), "2 years"), 
                                  labels = c(2010,2012,2014,2016,2018,2020), expand = c(0,0))
    
  
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  ggsave(paste("LTS_plots/hts_timeseries_",NRS,".png", sep = ""), P, width = 6, height = 8,  dpi=600)
  
})

##############################################################################################################################################################
# plot a climatology at each station

climatologyplots <- lapply(X = stations, FUN = function(f) {
  nc <- substitute(nc_open(FORM[[1]]), list(FORM = f))
  nc <- eval(nc)
  
  refDate <- ymd_hms("1950-01-01 00:00:00") # netcdf time in days since 1950
  
  df <- data.frame("Time" = ncvar_get(nc, "TIME"),
                   "Depth" = ncvar_get(nc, "DEPTH"),
                   "Temp" = ncvar_get(nc, "TEMP")) %>%
    mutate(Date = as.POSIXct(Time * 3600 * 24, origin = refDate),
           DOY = as.numeric(strftime(Date, format = "%j")))
  
  nc_close(nc)
  
  df_sum <- df %>% mutate(DepthBin = cut(Depth, breaks = c(0, 20, 40, 60, 80, 100)),
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
  
  P <- ggplot() + geom_smooth(data = df_sum, aes(DOY, Temp)) + 
    facet_grid(DepthBin ~ ., scales = "free") + labs(y = expression("Temperature ("* degree *"C)"),
                                    x = "Day of Year") +
    scale_x_continuous(expand = c(0,0))
 
  NRS <- substitute(as.character(FORM)[[2]], list(FORM = f))
  NRS <- eval(NRS)
  
  ggsave(paste("LTS_plots/hts_climatology_",NRS,".png", sep = ""), P, width = 6, height = 8,  dpi=600)
  
})
