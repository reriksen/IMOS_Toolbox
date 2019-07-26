library(readr)
library(dplyr)

source("fIMOS_MatchMODIS.R")

# Get the latitude/longitude and date from the file
filename <- "TestData_IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance_HTL.csv"
dat <- read_csv(filename)
dat <- dat[floor(runif(20,1, length(dat$Latitude))),] # Subset to 20 random values for subsetting

# # Test Date
# dat <- dat %>% 
#   select(-c("Day", "Month", "Year")) %>% 
#   rename(Date = SampleDateLocal)

# pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
#         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
#         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")

pr <- c("sst", "chl_oc3")
res_temp <- "12mNy"

dat <- fIMOS_MatchMODIS(dat, pr, res_temp)