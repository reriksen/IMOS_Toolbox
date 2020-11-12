## IMOS BGC Combined Water Quality Parameters
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Aug 2020
## Updated: 
## 24 Sept 2020 (Written to Git)
## 6th October 2020

suppressPackageStartupMessages({
  library(lubridate)
  library(data.table)
  library(tidyverse)
})

source("IMOS_Plankton_functions.R")

rawD <- "RawData"
outD <- "Output"

################################
## Bring in data for combined water quality
################################

# Each trip and depth combination for water quality parameters
# the number of rows in this table should equal that in comb, if not look out for duplicates and replicates
NRSTrips <- get_NRSTrips()

# Hydrochemistry data 
Chemistry <- read_csv(paste0(rawD,.Platform$file.sep,"chemistry.csv"), na = c("(null)", NaN)) %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M, Silicate_umol_L = SILICATE_UMOL_PER_L, Nitrate_umol_L =  NITRATE_UMOL_PER_L,
         Phosphate_umol_L =  PHOSPHATE_UMOL_PER_L, Salinity = SALINITY, 
         Ammonium_umol_L =  AMMONIUM_UMOL_PER_L,
         Nitrite_umol_L =  NITRITE_UMOL_PER_L,
         TCO2_umol_kg =  TCO2_UMOL_PER_KG,
         TAlkalinity_umol_kg =  TALKALINITY_UMOL_PER_KG,
         Oxygen_umol_L =  OXYGEN_UMOL_PER_L) %>% 
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         Silicate_umol_L = ifelse(SILICATE_FLAG %in% c(3,4,9), NA, Silicate_umol_L), # remove all data flagged as bad or probably bad
         Phosphate_umol_L = ifelse(PHOSPHATE_FLAG %in% c(3,4,9), NA, Phosphate_umol_L),
         Ammonium_umol_L = ifelse(AMMONIUM_FLAG %in% c(3,4,9), NA, Ammonium_umol_L),
         Nitrate_umol_L = ifelse(NITRATE_FLAG %in% c(3,4,9), NA, Nitrate_umol_L),
         Nitrite_umol_L = ifelse(NITRITE_FLAG %in% c(3,4,9), NA, Nitrite_umol_L),
         Oxygen_umol_L = ifelse(OXYGEN_FLAG %in% c(3,4,9), NA, Oxygen_umol_L),
         TCO2_umol_kg = ifelse(CARBON_FLAG %in% c(3,4,9), NA, TCO2_umol_kg),
         TAlkalinity_umol_kg = ifelse(ALKALINITY_FLAG %in% c(3,4,9), NA, TAlkalinity_umol_kg),
         Salinity = ifelse(SALINITY_FLAG %in% c(3,4,9), NA, Salinity)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(Silicate_umol_L = mean(Silicate_umol_L, na.rm = TRUE), # some replicated samples from error picking up PHB data, needs addressing in database
            Phosphate_umol_L = mean(Phosphate_umol_L, na.rm = TRUE),
            Ammonium_umol_L = mean(Ammonium_umol_L, na.rm = TRUE),
            Nitrate_umol_L = mean(Nitrate_umol_L, na.rm = TRUE),
            Nitrite_umol_L = mean(Nitrite_umol_L, na.rm = TRUE),
            Oxygen_umol_L = mean(Oxygen_umol_L, na.rm = TRUE),
            TCO2_umol_kg = mean(TCO2_umol_kg, na.rm = TRUE),
            TAlkalinity_umol_kg = mean(TAlkalinity_umol_kg, na.rm = TRUE),
            Salinity = mean(Salinity, na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  mutate_all(~ replace(., is.na(.), NA)) %>% 
  untibble()

# Zooplankton biomass
ZBiomass <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_biomass.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_CODE, Biomass_mgm3 = BIOMASS_MGM3) %>% 
  mutate(SampleDepth_m = "WC",
         NRScode = gsub('^.{3}|.{9}$', '', NRScode)) %>% 
  untibble()

# Pigments data
Pigments <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_pigments.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4)) %>% 
  filter(QC_FLAG %in% c(0,1,2,5,8)) %>% # keep data flagged as good
  untibble()

# Flow cytometry picoplankton data
Pico <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_picoplankton.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE,
         SampleDepth_m = SAMPLE_DEPTH_M, Prochlorococcus_cells_ml = PROCHLOROCOCCUS_CELLSPERML, Synecochoccus_cells_ml = SYNECOCHOCCUS_CELLSPERML, 
         Picoeukaryotes_cells_ml = PICOEUKARYOTES_CELLSPERML) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         Prochlorococcus_cells_ml = ifelse(PROCHLOROCOCCUS_FLAG %in% c(3,4,9), NA, Prochlorococcus_cells_ml), # remove bad data
         Synecochoccus_cells_ml = ifelse(SYNECOCHOCCUS_FLAG %in% c(3,4,9), NA, Synecochoccus_cells_ml),
         Picoeukaryotes_cells_ml = ifelse(PICOEUKARYOTES_FLAG %in% c(3,4,9), NA, Picoeukaryotes_cells_ml)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(Prochlorococcus_cells_ml = mean(Prochlorococcus_cells_ml, na.rm = TRUE), # mean of replicates
            Synecochoccus_cells_ml = mean(Synecochoccus_cells_ml, na.rm = TRUE),
            Picoeukaryotes_cells_ml = mean(Picoeukaryotes_cells_ml, na.rm = TRUE),
            .groups = "drop") %>% 
  untibble()

# Total suspended solid data
TSS <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_TSS.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = SAMPLE_DEPTH_M, TSS_mg_L = TSS_MG_PER_L, 
         InorganicFraction_mg_L = INORGANIC_FRACTION_MG_PER_L, 
         OrganicFraction_mg_L = ORGANIC_FRACTION_MG_PER_L, Secchi_m = SECCHI_DEPTH_M) %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m),
         NRScode = substring(NRScode,4),
         TSS_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, TSS_mg_L), # remove bad data
         InorganicFraction_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, InorganicFraction_mg_L),
         OrganicFraction_mg_L = ifelse(TSS_FLAG %in% c(3,4,9), NA, OrganicFraction_mg_L)) %>%
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(TSS_mg_L = mean(TSS_mg_L, na.rm = TRUE), # mean of replicates
            InorganicFraction_mg_L = mean(InorganicFraction_mg_L, na.rm = TRUE),
            OrganicFraction_mg_L = mean(OrganicFraction_mg_L, na.rm = TRUE),
            .groups = "drop") %>% 
  untibble()

# Secchi Disc        
Secchi <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_TSS.csv"), na = "(null)") %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = SAMPLE_DEPTH_M, Secchi_m = SECCHI_DEPTH_M) %>%
  select(NRScode, Secchi_m, SampleDepth_m) %>% 
  distinct() %>%
  mutate(SampleDepth_m = "WC",
         NRScode = substring(NRScode,4))

# CTD Cast Data
CTD <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                col_types = cols(PRES = col_double(), # columns start with nulls so tidyverse annoyingly assigns col_logical()
                                 PAR = col_double(),
                                 SPEC_CNDC = col_double())) %>% 
  rename(NRScode = NRS_TRIP_CODE, SampleDepth_m = PRES_REL, CTDDensity_kgm3 = DENS, 
         CTDTemperature = TEMP, CTDPAR_umolm2s = PAR,
         CTDConductivity_sm = CNDC, CTDSpecificConductivity_Sm = SPEC_CNDC, 
         CTDSalinity = PSAL, CTDTurbidity_ntu = TURB, CTDChlF_mgm3 = CHLF) %>%
    mutate(SampleDepth_m = as.character(round(SampleDepth_m, 0))) %>% 
  untibble()

notrips <-  read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                     col_types = cols(PRES = col_double(), # columns start with nulls so tidyverse annoyingly assigns col_logical()
                                      PAR = col_double(),
                                      SPEC_CNDC = col_double())) %>% select(NRS_TRIP_CODE) %>% distinct()


# Combined BGC data for each station at the sample depth
BGC <- NRSTrips %>% 
  left_join(ZBiomass %>% 
              select(NRScode, SampleDepth_m, Biomass_mgm3), by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Secchi,  by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Chemistry, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Pico, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(Pigments, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(TSS, by = c("NRScode", "SampleDepth_m")) %>%
  left_join(CTD, by = c("NRScode", "SampleDepth_m")) 

# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- BGC %>% 
  group_by(NRScode, SampleDepth_m) %>% 
  summarise(n = n(),
            .groups = "drop")

# Check
max(test$n)

# save to github
fwrite(BGC, file = paste0(outD,.Platform$file.sep,"NRS_CombinedWaterQuality.csv"), row.names = FALSE)
