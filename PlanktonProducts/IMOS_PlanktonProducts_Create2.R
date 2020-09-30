## IMOS plankton data products Indices 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Sept 2020
## Updated: 
## 1 Oct 2020 (Written to Git)

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

Trips <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_trips.csv"), na = "(null)") %>% 
  dplyr::rename("Station" = "STATION_NAME", "Latitude" = "Y_COORD", "Longitude" = "X_COORD", "SampleDateLocal" = "TRIP_START_DATETIME_LOCAL", 
                "NRScode" = "NRS_CODE", "StationDepth_m" = "STATION_DEPTH") %>%
  select(-SAMPLEDEPTH_M) %>%
  mutate(Year = year(SampleDateLocal),
         Month = month(SampleDateLocal),
         Day = day(SampleDateLocal),
         Time_24hr = str_sub(SampleDateLocal, -8, -1), # hms doesn"t seem to work on 00:00:00 times
        SampleDateLocal = as.character(SampleDateLocal)) %>% distinct()

# SST and Chlorophyll from CTD
ctd <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_CTD.csv"), na = "(null)",
                col_types = cols(PRES = col_double(), # columns start with  nulls so tidyverse annoyingly assigns col_logical()
                                 PAR = col_double(),
                                 SPEC_CNDC = col_double())) %>% 
  dplyr::rename("NRScode" = "NRS_TRIP_CODE", "SampleDepth_m" = "PRES_REL", "CTDDensity_kgm3" = "DENS", "CTDTemperature" = "TEMP", "CTDPAR_umolm2s" = "PAR",
                "CTDConductivity_sm" = "CNDC", "CTDSpecificConductivity_Sm" = "SPEC_CNDC", "CTDSalinity" = "PSAL", "CTDTurbidity_ntu" = "TURB", "CTDChlF_mgm3" = "CHLF") %>%
  mutate(SampleDepth_m = as.character(SampleDepth_m, 0)) %>% 
  filter(SampleDepth_m <11) %>% group_by(NRScode) %>% summarise(CTD_SST_C = mean(CTDTemperature, na.rm = TRUE),
                                                                CTDChlF_mgm3 = mean(CTDChlF_mgm3, na.rm = TRUE)) %>% untibble()


# Access satellite data for the sample dates using the IMOS_Toolbox

source("https://raw.githubusercontent.com/jaseeverett/IMOS_Toolbox/master/fIMOS_MatchAltimetry.R") 
source("https://raw.githubusercontent.com/jaseeverett/IMOS_Toolbox/master/fIMOS_MatchMODIS.R")

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

library(ncdf4)

dat <- Trips %>% 
  rename(Date = SampleDateLocal)

# Possible products
# pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
#         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
#         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")

pr <- c("sst", "chl_oc3")
res_temp <- "1d"
res_spat <- 10 # Return the average of res_spat x res_spat pixels

# Get MODIS Data
dat <- fIMOS_MatchMODIS(dat, pr, res_temp, res_spat)

# Get Altimetry (Gridded sea level anomaly, Gridded sea level, Surface geostrophic velocity)
dat <- fIMOS_MatchAltimetry(dat, res_spat)

# nutrient data
nuts <- chemistry %>% group_by(NRScode) %>% summarise(Silicate_umol_L = mean(Silicate_umol_L, na.rm = TRUE),
                                                      Phosphate_umol_L = mean(Phosphate_umol_L, na.rm = TRUE),
                                                      Ammonium_umol_L = mean(Ammonium_umol_L, na.rm = TRUE),
                                                      Nitrate_umol_L = mean(Nitrate_umol_L, na.rm = TRUE),
                                                      Nitrite_umol_L = mean(Nitrite_umol_L, na.rm = TRUE),
                                                      Oxygen_umol_L = mean(Oxygen_umol_L, na.rm = TRUE),
                                                      TCO2_umol_kg = mean(TCO2_umol_kg, na.rm = TRUE),
                                                      TAlkalinity_umol_kg = mean(TAlkalinity_umol_kg, na.rm = TRUE),
                                                      Salinity_umol_L = mean(Salinity, na.rm = TRUE))

# Total zoop abundance
zoodata <-  NRSZsamp %>% left_join(NRSZdat, by = "Sample")

TZoo <-  zoodata %>% group_by(NRScode) %>% summarise(ZoopAbundance_m3 = sum(ZAbund_m3, na.rm = TRUE))
TCope <- zoodata %>% filter(Copepod == 'COPEPOD') %>% group_by(NRScode) %>% summarise(CopeAbundance_m3 = sum(ZAbund_m3, na.rm = TRUE))

# Bring in copepod information table with sizes etc.
Zinfo <- read_csv(paste0(rawD,.Platform$file.sep,"taxon_info.csv"), na = "(null)") %>% 
  dplyr::rename( "TaxonName" = "TAXON_NAME") %>% untibble()

ACopeSize <- zoodata %>% filter(Copepod == 'COPEPOD') %>%  inner_join(Zinfo %>% select(SIZE_AVE_MM, TaxonName, DIET), by = "TaxonName") %>%
  mutate(abunSize = SIZE_AVE_MM * ZAbund_m3, 
         DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>%
  group_by(NRScode) %>% summarise(AvgTotalLengthCopepod_mm = sum(abunSize, na.rm = TRUE)/sum(ZAbund_m3, na.rm = TRUE))

HCrat <- zoodata %>% filter(Copepod == 'COPEPOD') %>%  inner_join(Zinfo %>% select(TaxonName, DIET), by = "TaxonName") %>%
  mutate(DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>% drop_na() %>%
  select(NRScode, DIET, ZAbund_m3) %>% 
  group_by(NRScode, DIET) %>% summarise(sumdiet = sum(ZAbund_m3 , na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumdiet, names_from = DIET) %>%
  mutate(HerbivoreCarnivoreCopepodRatio = CO / (CO + CC)) %>% untibble

# Diversity, evenness etc.     

# Bring in plankton data
NRSZcount <- read_csv(paste0(rawD,.Platform$file.sep,"NRS_zoop_count_raw.csv"), na = "(null)") %>%
  dplyr::rename("TaxonName" = "TAXON_NAME", "Copepod" = "TAXON_GROUP", "TaxonGroup" = "TAXON_GRP01", "NRScode" = "NRS_CODE",
                "Genus" = "GENUS", "Species" = "SPECIES", "TaxonCount" = "TAXON_COUNT")

zooCount <-  NRSZsamp %>% left_join(NRSZcount, by = "NRScode")

n <-  zooCount %>% filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% summarise(NoCopepodSpecies_Sample = n())
ShannonCopepodDiversity <- zooCount %>% filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% summarise(ZCount = sum(TaxonCount, na.rm = TRUE)) %>%
  pivot_wider(values_from = ZCount, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')
CopepodEvenness <- n %>% cbind(ShannonCopepodDiversity) %>% mutate(CopepodEvenness = ShannonCopepodDiversity / log(NoCopepodSpecies_Sample))

# Total Phyto abundance
phytodata <-  NRSPsamp %>% left_join(NRSPdat, by = "Sample") %>% filter(TaxonGroup != 'Other')

PhytoC <- phytodata %>% select(NRScode, TaxonGroup, Cells_L, Biovolume_uM3_L) %>% 
  mutate(BV_Cell = Biovolume_uM3_L / Cells_L, # biovolume of one cell
         Carbon = ifelse(TaxonGroup == 'Dinoflagellate', 0.76*(BV_Cell)^0.819, # conversion to Carbon based on taxongroup and biovolume of cell
                         ifelse(TaxonGroup == 'Ciliate', 0.22*(BV_Cell)^0.939,
                                ifelse(TaxonGroup == 'Cyanobacteria', 0.2, 0.288*(BV_Cell)^0.811 ))),
         Carbon_L = Cells_L * Carbon) %>% # Carbon per litre
  group_by(NRScode) %>% summarise(PhytoBiomassCarbon_pg_L = sum(Carbon_L))

TPhyto <-  phytodata %>% group_by(NRScode) %>% summarise(AbundancePhyto_cells_L = sum(Cells_L, na.rm = TRUE))

DDrat <- phytodata %>%  filter(TaxonGroup %in% c('Centric diatom', "Pennate diatom", 'Dinoflagellate')) %>% 
  mutate(TaxonGroup = recode(TaxonGroup, 'Centric diatom' = 'Diatom', 'Pennate diatom' = 'Diatom')) %>%
  select(NRScode, TaxonGroup, Cells_L) %>% 
  group_by(NRScode, TaxonGroup) %>% summarise(sumTG = sum(Cells_L, na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumTG, names_from = TaxonGroup) %>%
  mutate(DiatomDinoflagellateRatio = Diatom / (Diatom + Dinoflagellate)) %>% untibble

AvgCellVol <- phytodata %>% filter(!is.na(Biovolume_uM3_L)) %>% 
  group_by(NRScode) %>% summarise(AvgCellVol_um3 = mean(sum(Biovolume_uM3_L)/sum(Cells_L)))

# Diversity (phyto, diatoms, dinos)
# stick to abundance data here or we lose all the data that Pru counted which we don't have counts for.

np <-  phytodata %>% filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% summarise(NoPhytoSpecies_Sample = n())
ShannonPhytoDiversity <- phytodata %>% filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% summarise(Pdata = sum(Cells_L, na.rm = TRUE)) %>%
  pivot_wider(values_from = Pdata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')
PhytoEven <- np %>% cbind(ShannonPhytoDiversity) %>% mutate(PhytoEvenness = ShannonPhytoDiversity / log(NoPhytoSpecies_Sample))

ndia <-  phytodata %>% filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% summarise(NoDiatomSpecies_Sample = n())
ShannonDiatomDiversity <- phytodata %>% filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom')  & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% summarise(Diadata = sum(Cells_L, na.rm = TRUE)) %>%
  pivot_wider(values_from = Diadata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')
DiaEven <- ndia %>% cbind(ShannonDiatomDiversity) %>% mutate(DiatomEvenness = ShannonDiatomDiversity / log(NoDiatomSpecies_Sample))

ndino <-  phytodata %>% filter(TaxonGroup == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode) %>% summarise(NoDinoSpecies_Sample = n())
ShannonDinoDiversity <- phytodata %>% filter(TaxonGroup  == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(NRScode, TaxonName) %>% summarise(Dinodata = sum(Cells_L, na.rm = TRUE)) %>%
  pivot_wider(values_from = Dinodata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-NRScode) %>%
  diversity('shannon')
DinoEven <- ndino %>% cbind(ShannonDinoDiversity) %>% mutate(DinoflagellateEvenness = ShannonDinoDiversity / log(NoDinoSpecies_Sample))


# make indices table (nrows must always equal nrows of Trips)
Indices <-  Trips  %>%
  left_join(TZoo, by = ("NRScode")) %>%
  left_join(TCope, by = ("NRScode")) %>%
  left_join(biomass %>% select(-SampleDepth_m), by = ("NRScode")) %>%
  left_join(ACopeSize, by = ("NRScode")) %>%
  left_join(HCrat %>% select(-c('CO', 'CC')), by = ("NRScode")) %>% 
  left_join(CopepodEvenness,  by = ("NRScode")) %>%
  left_join(PhytoC, by = ("NRScode")) %>%
  left_join(TPhyto, by = ("NRScode")) %>%
  left_join(DDrat %>% select(-c('Diatom', 'Dinoflagellate')), by = ("NRScode")) %>%
  left_join(AvgCellVol, by = ("NRScode")) %>%
  left_join(PhytoEven, by = ("NRScode")) %>%
  left_join(DiaEven, by = ("NRScode")) %>%
  left_join(DinoEven, by = ("NRScode")) %>%   
  left_join(ctd, by = ("NRScode")) %>%
  left_join(dat %>% select(NRScode, sst_1d, chl_oc3_1d, GSLA, GSL, UCUR, VCUR), by = ("NRScode")) %>%
  left_join(nuts, by = ("NRScode")) 
  


# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- Indices %>% group_by(NRScode) %>% summarise(n = n())


