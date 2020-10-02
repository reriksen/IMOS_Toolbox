suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(reshape)
  library(vegan)
  library(data.table)  
  library(tidyverse)
})

rawD <- "RawData"
outD <- "Output"

untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense

# uses mostly the same raw data from IMOS_PlanktonProducts_Create.R

# ensure we have all trips accounted for 
# note there are circumstances where a trip won't have a phyto and a zoo samples due to loss of sample etc.

TripsIMOS <- read_csv("RawData/nrs_trips_imosind.csv", na = "(null)") %>% 
  dplyr::rename("Station" = "STATION_NAME", "Latitude" = "Y_COORD", "Longitude" = "X_COORD", "SampleDateLocal" = "TRIP_START_DATETIME_LOCAL", 
                "NRScode" = "NRS_CODE", "SST_C" = "SST", "ChlorophyllMonthlyClimatology_mg_m3" = "CHL_AVG", "ChlorophyllSatellite_mg_m3" = "CHL") %>%
  mutate(Year = year(SampleDateLocal),
         Month = month(SampleDateLocal),
         Day = day(SampleDateLocal),
         Time_24hr = str_sub(SampleDateLocal, -8, -1), # hms doesn"t seem to work on 00:00:00 times
         SampleDateLocal = as.character(SampleDateLocal)) %>% distinct()

IMSONRS_ind <- TripsIMOS[,c(1:5,9:12,6:8)] %>%
  left_join(TZoo, by = ("NRScode")) %>%
  left_join(TCope, by = ("NRScode")) %>%
  left_join(ACopeSize, by = ("NRScode")) %>%
  left_join(HCrat %>% select(-c('CO', 'CC')), by = ("NRScode")) %>% 
  left_join(CopepodEvenness,  by = ("NRScode")) %>%
  left_join(biomass %>% select(-SampleDepth_m), by = ("NRScode")) %>%
  left_join(nuts[,c(1,10,2,5,4,3,8,9,7)], by = ("NRScode")) %>%
  left_join(PhytoC, by = ("NRScode")) %>%
  left_join(TPhyto, by = ("NRScode")) %>%
  left_join(DDrat %>% select(-c('Diatom', 'Dinoflagellate')), by = ("NRScode")) %>%
  left_join(AvgCellVol, by = ("NRScode")) %>%
  left_join(PhytoEven, by = ("NRScode")) %>%
  left_join(DiaEven, by = ("NRScode")) %>%
  left_join(DinoEven, by = ("NRScode"))   


