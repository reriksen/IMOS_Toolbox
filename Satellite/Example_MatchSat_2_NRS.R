library(readr)
library(dplyr)

source("fIMOS_MatchGHRSST.R")
source("fIMOS_MatchMODIS.R")
source("fIMOS_MatchAltimetry.R")

# If on Windows you will need to install a development 
# version of ncdf4 which allows the use of OpenDAP
if(.Platform$OS.type == "windows") {
  warning("It looks like you are on a Windows PC - You will need to install a 
  development version of ncdf4 which allows the use of OpenDAP. Please 
  run devtools::install_github('mdsumner/ncdf4') to install or 
  see 'https://github.com/mdsumner/ncdf4' for more information.")
}

# install.packages("devtools")
# devtools::install_github("mdsumner/ncdf4")

# Get the latitude/longitude and date from the file
filename <- "TestData_IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance_HTL.csv"
dat <- read_csv(filename)
dat <- dat[floor(runif(20,1, length(dat$Latitude))),] # Subset to 20 random values for testing

dat <- dat %>% 
  rename(Date = SampleDateLocal)


# Get GHRSST SST Data
# Possible products to download are: 
# dt_analysis, l2p_flags, quality_level, satellite_zenith_angle, sea_ice_fraction, sea_ice_fraction_dtime_from_sst, 
# sea_surface_temperature, sea_surface_temperature_day_night, sses_bias, sses_count,sses_standard_deviation,
# sst_count, sst_dtime, sst_mean, sst_standard_deviation, wind_speed, wind_speed_dtime_from_sst,
res_temp <- "1d"
res_spat <- 10 # Return the average of res_spat x res_spat pixels
pr <- ("sea_surface_temperature")
dat <- fIMOS_MatchGHRSST(dat, pr, res_temp, res_spat)


# Get MODIS Data
# Possible products
# pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
#         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
#         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")

pr <- c("sst", "chl_oci")
res_temp <- "1d"
res_spat <- 10 # Return the average of res_spat x res_spat pixels

dat <- fIMOS_MatchMODIS(dat, pr, res_temp, res_spat)


# Get Altimetry
dat <- fIMOS_MatchAltimetry(dat, res_spat)