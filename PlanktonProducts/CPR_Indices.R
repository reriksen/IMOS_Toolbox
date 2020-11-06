## IMOS plankton data products Indices 
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: Sept 2020
## Updated: 
## 1 Oct 2020 (Written to Git)

suppressPackageStartupMessages({
  library(readr)
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

cprTrips <- read_csv(paste0(rawD,.Platform$file.sep,"PSampCPR.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, Route = ROUTE, Region = REGION, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC) %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1), # hms doesn"t seem to work on 00:00:00 times
         SampleDateUTC = as.character(SampleDateUTC)) %>% 
  select(c(Sample, Latitude:Time_24hr, Region, Route))

cprProps <- read_csv(paste0(rawD,.Platform$file.sep,"AllSampCPR.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, ChlorophyllMonthlyClimatology_mg_m3 = CHL_AVG, ChlorophyllSatellite_mg_m3 = CHL, WaterDepth_m = WATERDEPTH_M)

# accessing the satelitte data from MODIS
datcpr <- cprTrips %>% 
  rename(Date = SampleDateUTC) 

# Possible products
# pr <- c("sst_quality", "sst", "picop_brewin2012in", "picop_brewin2010at", "par", 
#         "owtd", "npp_vgpm_eppley_oc3", "npp_vgpm_eppley_gsm", "nanop_brewin2012in",
#         "nanop_brewin2010at", "l2_flags", "ipar", "dt", "chl_oc3", "chl_gsm", "K_490")

pr <- c("sst", "chl_oc3")
res_temp <- "1d"
res_spat <- 10 # Return the average of res_spat x res_spat pixels

# Get MODIS Data
datcpr <- fIMOS_MatchMODIS(datcpr, pr, res_temp, res_spat)

# Get Altimetry (Gridded sea level anomaly, Gridded sea level, Surface geostrophic velocity)
datcpr <- fIMOS_MatchAltimetry(dat, res_spat)

# Total zoop abundance
zoodatacpr <-  cprZsamp %>% left_join(cprZdat, by = "Sample")

TZoocpr <-  zoodatacpr %>% group_by(Sample) %>% summarise(ZoopAbundance_m3 = sum(ZAbun_m3, na.rm = TRUE))
TCopecpr <- zoodatacpr %>% filter(Copepod == 'COPEPOD') %>% group_by(Sample) %>% summarise(CopeAbundance_m3 = sum(ZAbun_m3, na.rm = TRUE))

# Bring in copepod information table with sizes etc.
ACopeSizeCpr <- zoodatacpr %>% filter(Copepod == 'COPEPOD') %>%  inner_join(Zinfo %>% select(SIZE_AVE_MM, TaxonName, DIET), by = "TaxonName") %>%
  mutate(abunSize = SIZE_AVE_MM * ZAbun_m3, 
         DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>%
  group_by(Sample) %>% summarise(AvgTotalLengthCopepod_mm = sum(abunSize, na.rm = TRUE)/sum(ZAbun_m3, na.rm = TRUE))

HCratCpr <- zoodatacpr %>% filter(Copepod == 'COPEPOD') %>%  inner_join(Zinfo %>% select(TaxonName, DIET), by = "TaxonName") %>%
  mutate(DIET = ifelse(DIET == 'CC', 'CC', 'CO')) %>% drop_na() %>%
  select(Sample, DIET, ZAbun_m3) %>% 
  group_by(Sample, DIET) %>% summarise(sumdiet = sum(ZAbun_m3 , na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumdiet, names_from = DIET) %>%
  mutate(HerbivoreCarnivoreCopepodRatio = CO / (CO + CC)) %>% untibble

CPRbiomass <- read_csv(paste0(rawD,.Platform$file.sep,"CPR_biomass.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, Biomass_mgCarbon_m3 = BIOMASS_MG_M3) %>% untibble()

# Diversity, evenness etc.     

# Bring in plankton data
CPRZcount <- read_csv(paste0(rawD,.Platform$file.sep,"CPR_zoo_count_raw.csv"), na = "(null)") %>%
  rename(TaxonName = TAXON_NAME, Copepod = TAXON_GROUP, TaxonGroup = TAXON_GRP01, Sample = SAMPLE,
                Genus= GENUS, Species = SPECIES, TaxonCount = COUNTS)

zooCountCpr <-  cprZsamp %>% left_join(CPRZcount, by = "Sample")

nCPR <-  zooCountCpr %>% filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoCopepodSpecies_Sample = n())
ShannonCopepodDiversityCPR <- zooCountCpr %>% filter(Copepod == 'COPEPOD' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% summarise(ZCount = sum(TaxonCount, na.rm = TRUE)) %>%
  pivot_wider(values_from = ZCount, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')
CopepodEvennessCPR <- nCPR %>% cbind(ShannonCopepodDiversityCPR) %>% mutate(CopepodEvenness = ShannonCopepodDiversityCPR / log(NoCopepodSpecies_Sample))

# Total Phyto abundance
phytodatacpr <-  cprPsamp %>% left_join(cprPdat, by = "Sample") %>% filter(TaxonGroup != 'Other')

PhytoCcpr <- phytodatacpr %>% select(Sample, TaxonGroup, PAbun_m3, BioVolume_um3_m3) %>% 
  mutate(BV_Cell = BioVolume_um3_m3 / PAbun_m3, # biovolume of one cell
         Carbon = ifelse(TaxonGroup == 'Dinoflagellate', 0.76*(BV_Cell)^0.819, # conversion to Carbon based on taxongroup and biovolume of cell
                         ifelse(TaxonGroup == 'Ciliate', 0.22*(BV_Cell)^0.939,
                                ifelse(TaxonGroup == 'Cyanobacteria', 0.2, 0.288*(BV_Cell)^0.811 ))),
         Carbon_m3 = PAbun_m3 * Carbon) %>% # Carbon per m3
  group_by(Sample) %>% summarise(PhytoBiomassCarbon_pg_m3 = sum(Carbon_m3))

TPhytoCpr <-  phytodatacpr %>% group_by(Sample) %>% summarise(AbundancePhyto_cells_m3 = sum(PAbun_m3, na.rm = TRUE))

DDratcpr <- phytodatacpr %>%  filter(TaxonGroup %in% c('Centric diatom', "Pennate diatom", 'Dinoflagellate')) %>% 
  mutate(TaxonGroup = recode(TaxonGroup, 'Centric diatom' = 'Diatom', 'Pennate diatom' = 'Diatom')) %>%
  select(Sample, TaxonGroup, PAbun_m3) %>% 
  group_by(Sample, TaxonGroup) %>% summarise(sumTG = sum(PAbun_m3, na.rm = TRUE)) %>% 
  pivot_wider(values_from = sumTG, names_from = TaxonGroup) %>%
  mutate(DiatomDinoflagellateRatio = Diatom / (Diatom + Dinoflagellate)) %>% untibble

AvgCellVolcpr <- phytodatacpr %>% filter(!is.na(BioVolume_um3_m3)) %>% 
  group_by(Sample) %>% summarise(AvgCellVol_um3 = mean(sum(BioVolume_um3_m3)/sum(PAbun_m3)))

# Diversity (phyto, diatoms, dinos)
# stick to abundance data here as otherwise we have FOV counts

npcpr <-  phytodatacpr %>% filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoPhytoSpecies_Sample = n())
ShannonPhytoDiversitycpr <- phytodatacpr %>% filter(TaxonGroup != 'Other' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% summarise(Pdata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Pdata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')
PhytoEvencpr <- npcpr %>% cbind(ShannonPhytoDiversitycpr) %>% mutate(PhytoEvenness = ShannonPhytoDiversitycpr / log(NoPhytoSpecies_Sample))

ndiacpr <-  phytodatacpr %>% filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom') & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoDiatomSpecies_Sample = n())
ShannonDiatomDiversitycpr <- phytodatacpr %>% filter(TaxonGroup %in% c('Centric diatom', 'Pennate diatom')  & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% summarise(Diadata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Diadata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')
DiaEvencpr <- ndiacpr %>% cbind(ShannonDiatomDiversitycpr) %>% mutate(DiatomEvenness = ShannonDiatomDiversitycpr / log(NoDiatomSpecies_Sample))

ndinocpr <-  phytodatacpr %>% filter(TaxonGroup == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample) %>% summarise(NoDinoSpecies_Sample = n())
ShannonDinoDiversitycpr <- phytodatacpr %>% filter(TaxonGroup  == 'Dinoflagellate' & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(TaxonName = paste0(Genus," ", word(Species,1))) %>% # bin complexes 
  group_by(Sample, TaxonName) %>% summarise(Dinodata = sum(PAbun_m3, na.rm = TRUE)) %>%
  pivot_wider(values_from = Dinodata, names_from = TaxonName, values_fill = 0) %>% ungroup() %>%
  select(-Sample) %>%
  diversity('shannon')
DinoEvencpr <- ndinocpr %>% cbind(ShannonDinoDiversitycpr) %>% mutate(DinoflagellateEvenness = ShannonDinoDiversitycpr / log(NoDinoSpecies_Sample))

# make indices table (nrows must always equal nrows of Trips)
IndicesCPR <-  cprTrips  %>%
  left_join(TZoocpr, by = ("Sample")) %>%
  left_join(TCopecpr, by = ("Sample")) %>%
  left_join(ACopeSizeCpr, by = ("Sample")) %>%
  left_join(HCratCpr %>% select(-c('CO', 'CC')), by = ("Sample")) %>% 
  left_join(CPRbiomass, by = ("Sample")) %>%
  left_join(CopepodEvennessCPR,  by = ("Sample")) %>%
  left_join(PhytoCcpr, by = ("Sample")) %>%
  left_join(TPhytoCpr, by = ("Sample")) %>%
  left_join(DDratcpr %>% select(-c('Diatom', 'Dinoflagellate')), by = ("Sample")) %>%
  left_join(AvgCellVolcpr, by = ("Sample")) %>%
  left_join(PhytoEvencpr, by = ("Sample")) %>%
  left_join(DiaEvencpr, by = ("Sample")) %>%
  left_join(DinoEvencpr, by = ("Sample")) %>%   
  left_join(datcpr %>% select(Sample, sst_1d, chl_oc3_1d), by = ("Sample")) %>%  #add once run , GSLA, GSL, UCUR, VCUR
  select(-Sample)

# make indices table (nrows must always equal nrows of Trips) - old one for IMOS
IMOSCPR_ind <-  cprTrips  %>%
  left_join(cprProps, by = ("Sample")) %>%
  left_join(TZoocpr, by = ("Sample")) %>%
  left_join(TCopecpr, by = ("Sample")) %>%
  left_join(ACopeSizeCpr, by = ("Sample")) %>%
  left_join(HCratCpr %>% select(-c('CO', 'CC')), by = ("Sample")) %>% 
  left_join(CPRbiomass, by = ("Sample")) %>%
  left_join(CopepodEvennessCPR,  by = ("Sample")) %>%
  left_join(PhytoCcpr, by = ("Sample")) %>%
  left_join(TPhytoCpr, by = ("Sample")) %>%
  left_join(DDratcpr %>% select(-c('Diatom', 'Dinoflagellate')), by = ("Sample")) %>%
  left_join(AvgCellVolcpr, by = ("Sample")) %>%
  left_join(PhytoEvencpr, by = ("Sample")) %>%
  left_join(DiaEvencpr, by = ("Sample")) %>%
  left_join(DinoEvencpr, by = ("Sample")) %>%   
  select(-Sample)



# test table
# n should be 1, replicates or duplicate samples will have values > 1
test <- Indices %>% group_by(Latitude, Longitude, SampleDateUTC) %>% summarise(n = n())


