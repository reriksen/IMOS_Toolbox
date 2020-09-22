## IMOS plankton data products
## claire
## May 2020

suppressPackageStartupMessages({
library(tidyverse)
  library(lubridate)
library(reshape)
library(data.table)
})

#### CPR Phytoplankton #######################################################################################################################################################
# Bring in all CPR phytoplankton samples
cprPsamp <- read_csv("PSampCPR.csv", na = "(null)") %>% 
  dplyr::rename("Sample" = "SAMPLE", "Route" = "ROUTE", "Latitude" = "LATITUDE", "Longitude" = "LONGITUDE", "SampleDateUTC" = "SAMPLEDATEUTC") %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn't seem to work on 00:00:00 times

# Bring in plankton data
cprPdat <- read_csv("CPR_phyto_raw.csv", na = "(null)") %>%
  dplyr::rename("Sample" = "SAMPLE", "TaxonName" = "TAXON_NAME", "TaxonGroup" = "TAXON_GROUP", "Genus" = "GENUS", "Species" = "SPECIES", "PAbun_m3" = "PHYTO_ABUNDANCE_M3")

# Bring in Change Log
cprPcl <- read_csv("ChangeLogCPRP.csv", na = "(null)") %>%
  dplyr::rename("TaxonName" = "TAXON_NAME", "StartDate" = "START_DATE", "ParentName" = "PARENT_NAME")

#### CPR PHYTO RAW ####

cprRawP1 <- left_join(cprPsamp, cprPdat, by = "Sample") %>% select(c(2:10,14)) %>% arrange(-desc(TaxonName)) %>%
  mutate(TaxonName = ifelse(is.na(TaxonName), 'No taxa found', TaxonName)) # for segments where no phyto was found
cprRawP <- cprRawP1 %>% pivot_wider(names_from = TaxonName, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) %>%
  select(-'No taxa found')

fwrite(cprRawP, file = "CPR_phyto_raw_mat.csv", row.names = FALSE)

#### CPR PHYTO HTG ####

cprHTGP1 <- cprPdat %>% group_by(Sample, TaxonGroup) %>% summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE)) %>%
  filter(!TaxonGroup %in% c('Other','Coccolithophore', 'Diatom','Protozoa')) 
cprHTGP1 <-  cprPsamp %>% left_join(cprHTGP1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), 'Ciliate', TaxonGroup),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3)) %>% arrange(-desc(TaxonGroup))
cprHTGP <-  cprHTGP1 %>% pivot_wider(names_from = TaxonGroup, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprHTGP[,-1], file = "CPR_phyto_HTG_mat.csv", row.names = FALSE)

#### CPR PHYTO GENUS ####

# Check genus are effected by change log
clg <- cprPcl %>% mutate(genus1 = word(TaxonName, 1),
                         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species

cprGenP1 <- cprPdat %>% filter(!TaxonName %in% levels(as.factor(clg$TaxonName))) %>% group_by(Sample, Genus) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE)) %>% drop_na(Genus) 
cprGenP1 <- cprPsamp %>% left_join(cprGenP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), 'Acanthoica', Genus),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3))  %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(PAbun_m3 = sum(PAbun_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
cprGenP2 <- cprPdat %>% filter(TaxonName %in% levels(as.factor(clg$TaxonName))) %>% 
  left_join(cprPcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE)) 

for (i in 1:nlevels(cprGenP2$Genus)) {
  Dates <- as.data.frame(cprGenP2) %>% filter(Genus == Genus[i]) %>% slice(1)  %>% droplevels()
  gen <- as.data.frame(cprGenP2) %>% filter(Genus == Genus[i]) %>% droplevels() 
  gen <- cprPsamp %>% left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           PAbun_m3 = replace(PAbun_m3, StartDate>SampleDateUTC, -999), 
           PAbun_m3 = replace(PAbun_m3, StartDate<SampleDateUTC & is.na(PAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(PAbun_m3 = sum(PAbun_m3)) %>% as.data.frame()     
  cprGenP1 <- rbind(cprGenP1, gen)
  }

cprGenP1 <- cprGenP1 %>% group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(PAbun_m3 = max(PAbun_m3)) %>% arrange(-desc(Genus)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprGenP <-  cprGenP1 %>% pivot_wider(names_from = Genus, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprGenP, file = "CPR_phyto_genus_mat.csv", row.names = FALSE)

#### CPR PHYTO SPECIES ####

# Check at what level we need change log
cls <- cprPcl %>% mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

cprSpecP1 <- cprPdat %>% filter(!TaxonName %in% levels(as.factor(clg$TaxonName))
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species)) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE))
cprSpecP1 <- cprPsamp %>% left_join(cprSpecP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         TaxonName = ifelse(is.na(TaxonName), 'Paralia sulcata', TaxonName),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3))  %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(PAbun_m3 = sum(PAbun_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
cprSpecP2 <- cprPdat %>% filter(TaxonName %in% levels(as.factor(clg$TaxonName))
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species)) %>% 
  left_join(cprPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE)) 

for (i in 1:nlevels(cprSpecP2$TaxonName)) {
  Dates <- as.data.frame(cprSpecP2) %>% filter(TaxonName == TaxonName[i]) %>% slice(1)  %>% droplevels()
  spec <- as.data.frame(cprSpecP2) %>% filter(TaxonName == TaxonName[i]) %>% droplevels() 
  spec <- cprPsamp %>% left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           PAbun_m3 = replace(PAbun_m3, StartDate>SampleDateUTC, -999), 
           PAbun_m3 = replace(PAbun_m3, StartDate<SampleDateUTC & is.na(PAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
    summarise(PAbun_m3 = sum(PAbun_m3)) %>% as.data.frame()     
  cprSpecP1 <- rbind(cprSpecP1, spec)
}

cprSpecP1 <- cprSpecP1 %>% group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(PAbun_m3 = max(PAbun_m3)) %>% arrange(-desc(TaxonName)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprSpecP <-  cprSpecP1 %>% pivot_wider(names_from = TaxonName, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprSpecP, file = "CPR_phyto_species_mat.csv", row.names = FALSE)

#### CPR Zoopplankton #### ################################################################################################################################
# Bring in all CPR phytoplankton samples
cprZsamp <- read_csv("ZSampCPR.csv", na = "(null)") %>% 
  dplyr::rename("Sample" = "SAMPLE", "Route" = "ROUTE", "Latitude" = "LATITUDE", "Longitude" = "LONGITUDE", "SampleDateUTC" = "SAMPLEDATEUTC") %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn't seem to work on 00:00:00 times

# Bring in plankton data
cprZdat <- read_csv("CPR_zoo_raw.csv", na = "(null)") %>%
  dplyr::rename("Sample" = "SAMPLE", "TaxonName" = "TAXON_NAME", "Copepod" = "TAXON_GROUP", "TaxonGroup" = "TAXON_GRP01",
                "Genus" = "GENUS", "Species" = "SPECIES", "ZAbun_m3" = "ZOOP_ABUNDANCE_M3")

# Bring in Change Log
cprZcl <- read_csv("ChangeLogCPRZ.csv", na = "(null)") %>%
  dplyr::rename("TaxonName" = "TAXON_NAME", "StartDate" = "START_DATE", "ParentName" = "PARENT_NAME")

#### CPR ZOOP RAW ####

cprRawZ1 <- left_join(cprZsamp, cprZdat, by = "Sample") %>% select(c(2:10,15)) %>% arrange(-desc(TaxonName)) %>%
  mutate(TaxonName = ifelse(is.na(TaxonName), 'No taxa found', TaxonName))
cprRawZ <- cprRawZ1 %>% pivot_wider(names_from = TaxonName, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC))  %>%
  select(-'No taxa found')

fwrite(cprRawZ, file = "CPR_zoop_raw_mat.csv", row.names = FALSE)

#### CPR ZOOP HTG ####

cprHTGZ1 <- cprZdat %>% group_by(Sample, TaxonGroup) %>% summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE)) %>%
  filter(!TaxonGroup %in% c('Other')) 
cprHTGZ1 <-  cprZsamp %>% left_join(cprHTGZ1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), 'Copepod', TaxonGroup),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3)) %>% arrange(-desc(TaxonGroup))
cprHTGZ <-  cprHTGZ1 %>% pivot_wider(names_from = TaxonGroup, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprHTGZ[,-1], file = "CPR_zoop_HTG_mat.csv", row.names = FALSE)

#### CPR ZOOP GENUS ####

# Check genus are effected by change log
clgz <- cprZcl %>% mutate(genus1 = word(TaxonName, 1),
                         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species

cprGenZ1 <- cprZdat %>% filter(!TaxonName %in% levels(as.factor(clgz$TaxonName))) %>% group_by(Sample, Genus) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE)) %>% drop_na(Genus) 
cprGenZ1 <- cprZsamp %>% left_join(cprGenZ1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), 'Calanus', Genus),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3))  %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
cprGenZ2 <- cprZdat %>% filter(TaxonName %in% levels(as.factor(clgz$TaxonName))) %>% 
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE)) 

for (i in 1:nlevels(cprGenZ2$Genus)) {
  Datesz <- as.data.frame(cprGenZ2) %>% filter(Genus == Genus[i]) %>% slice(1)  %>% droplevels()
  genz <- as.data.frame(cprGenZ2) %>% filter(Genus == Genus[i]) %>% droplevels() 
  genz <- cprZsamp %>% left_join(genz, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Datesz$StartDate),
           Genus = replace(Genus, is.na(Genus), Datesz$Genus),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()     
  cprGenZ1 <- rbind(cprGenZ1, genz)
}

cprGenZ1 <- cprGenZ1 %>% group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3)) %>% arrange(-desc(Genus)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprGenZ <-  cprGenZ1 %>% pivot_wider(names_from = Genus, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprGenZ, file = "CPR_zoop_genus_mat.csv", row.names = FALSE)

#### CPR ZOOP COPEPODS ####

# Check at what level we need change log
clc <- cprZcl %>% mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

cprCop1 <- cprZdat %>% filter(!TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod =='COPEPOD'
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE))
cprCop1 <- cprZsamp %>% left_join(cprCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), 'Calanus Australis', Species), # avoids nulls in pivot
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3))  %>%  # avoids nulls in pivot
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
cprCop2 <- cprZdat %>% filter(TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod =='COPEPOD'
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE)) 

for (i in 1:nlevels(cprCop2$Species)) {
  Dates <- as.data.frame(cprCop2) %>% filter(Species == Species[i]) %>% slice(1)  %>% droplevels()
  copes <- as.data.frame(cprCop2) %>% filter(Species == Species[i]) %>% droplevels() 
  copes <- cprZsamp %>% left_join(copes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()     
  cprCop1 <- rbind(cprCop1, copes)
}

cprCop1 <- cprCop1 %>% group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3)) %>% arrange(-desc(Species)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprCop <-  cprCop1 %>% pivot_wider(names_from = Species, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprCop, file = "CPR_zoop_copes_mat.csv", row.names = FALSE)

#### CPR ZOOP NON-COPEPODS ####

# for non change logspecies

cprnCop1 <- cprZdat %>% filter(!TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod !='COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE))
cprnCop1 <- cprZsamp %>% left_join(cprnCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), 'Evadne spinifera', Species),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3))  %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
cprnCop2 <- cprZdat %>% filter(TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod !='COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>% 
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE)) 

for (i in 1:nlevels(cprnCop2$Species)) {
  Dates <- as.data.frame(cprnCop2) %>% filter(Species == Species[i]) %>% slice(1)  %>% droplevels()
  ncopes <- as.data.frame(cprnCop2) %>% filter(Species == Species[i]) %>% droplevels() 
  ncopes <- cprZsamp %>% left_join(ncopes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3)) %>% as.data.frame()     
  cprnCop1 <- rbind(cprnCop1, ncopes)
}

cprnCop1 <- cprnCop1 %>% group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3)) %>% arrange(-desc(Species)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprnCop <-  cprnCop1 %>% pivot_wider(names_from = Species, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprnCop, file = "CPR_zoop_noncopes_mat.csv", row.names = FALSE)

#### NRS Phytoplankton #### #################################################################################################################################

# Bring in all NRS phytoplankton samples
NRSPsamp <- read_csv("PSampNRS.csv", na = "(null)") %>% 
  dplyr::rename("Sample" = "SAMPLE", "Station" = "STATION", "Latitude" = "LATITUDE", "Longitude" = "LONGITUDE", "SampleDateLocal" = "SAMPLEDATE", "NRScode" = "NRS_CODE") %>%
  mutate(Year = year(SampleDateLocal),
         Month = month(SampleDateLocal),
         Day = day(SampleDateLocal),
         Time_24hr = str_sub(SampleDateLocal, -8, -1)) # hms doesn't seem to work on 00:00:00 times

# Bring in plankton data
NRSPdat <- read_csv("NRS_phyto_raw.csv", na = "(null)") %>%
  dplyr::rename("Sample" = "SAMPLE", "TaxonName" = "TAXON_NAME", "TaxonGroup" = "TAXON_GROUP", "Genus" = "GENUS", "Species" = "SPECIES", 
                "Cells_L" = "CELL_PER_LITRE", "Biovolume_uM3_L" = "BIOVOLUME_UM3_PER_L")

# Bring in Change Log
NRSPcl <- read_csv("ChangeLogNRSP.csv", na = "(null)") %>%
  dplyr::rename("TaxonName" = "TAXON_NAME", "StartDate" = "START_DATE", "ParentName" = "PARENT_NAME")

#### NRS PHYTO RAW ####

NRSRawP1 <- left_join(NRSPsamp, NRSPdat, by = "Sample") %>% select(c(1:11,15)) %>% arrange(-desc(TaxonName)) 
NRSRawP <- NRSRawP1 %>% pivot_wider(names_from = TaxonName, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSRawP, file = "NRS_phyto_raw_mat.csv", row.names = FALSE)

#### NRS PHYTO HTG ####

NRSHTGP1 <- NRSPdat %>% group_by(Sample, TaxonGroup) %>% summarise(Cells_L = sum(Cells_L, na.rm = TRUE)) %>%
  filter(!TaxonGroup %in% c('Other','Coccolithophore', 'Diatom','Protozoa')) 
NRSHTGP1 <-  NRSPsamp %>% left_join(NRSHTGP1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), 'Ciliate', TaxonGroup),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L)) %>% arrange(-desc(TaxonGroup))
NRSHTGP <-  NRSHTGP1 %>% pivot_wider(names_from = TaxonGroup, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSHTGP[,-1], file = "NRS_phyto_HTG_mat.csv", row.names = FALSE)

#### NRS PHYTO GENUS ####

# Check genus are effected by change log
nrslg <- NRSPcl %>% mutate(genus1 = word(TaxonName, 1),
                         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species

NRSGenP1 <- NRSPdat %>% filter(!TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% group_by(Sample, Genus) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE)) %>% drop_na(Genus) 
NRSGenP1 <- NRSPsamp %>% left_join(NRSGenP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), 'Acanthoica', Genus),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L))  %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(Cells_L = sum(Cells_L)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
NRSGenP2 <- NRSPdat %>% filter(TaxonName %in% levels(as.factor(nrslg$TaxonName))) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE)) 

for (i in 1:nlevels(NRSGenP2$Genus)) {
  Dates <- as.data.frame(NRSGenP2) %>% filter(Genus == Genus[i]) %>% slice(1)  %>% droplevels()
  gen <- as.data.frame(NRSGenP2) %>% filter(Genus == Genus[i]) %>% droplevels() 
  gen <- NRSPsamp %>% left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           Cells_L = replace(Cells_L, StartDate>SampleDateLocal, -999), 
           Cells_L = replace(Cells_L, StartDate<SampleDateLocal & is.na(Cells_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(Cells_L = sum(Cells_L)) %>% as.data.frame()     
  NRSGenP1 <- rbind(NRSGenP1, gen)
}

NRSGenP1 <- NRSGenP1 %>% group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(Cells_L = max(Cells_L)) %>% arrange(-desc(Genus)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSGenP <-  NRSGenP1 %>% pivot_wider(names_from = Genus, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSGenP, file = "NRS_phyto_genus_mat.csv", row.names = FALSE)

#### NRS PHYTO SPECIES ####

# Check at what level we need change log
nrsls <- NRSPcl %>% mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

NRSSpecP1 <- NRSPdat %>% filter(!TaxonName %in% levels(as.factor(nrsls$TaxonName))
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species)) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE))
NRSSpecP1 <- NRSPsamp %>% left_join(NRSSpecP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         TaxonName = ifelse(is.na(TaxonName), 'Paralia sulcata', TaxonName),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L))  %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(Cells_L = sum(Cells_L)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
NRSSpecP2 <- NRSPdat %>% filter(TaxonName %in% levels(as.factor(nrsls$TaxonName))
                                & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species)) %>% 
  left_join(NRSPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE)) 

for (i in 1:nlevels(NRSSpecP2$TaxonName)) {
  Dates <- as.data.frame(NRSSpecP2) %>% filter(TaxonName == TaxonName[i]) %>% slice(1)  %>% droplevels()
  spec <- as.data.frame(NRSSpecP2) %>% filter(TaxonName == TaxonName[i]) %>% droplevels() 
  spec <- NRSPsamp %>% left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           Cells_L = replace(Cells_L, StartDate>SampleDateLocal, -999), 
           Cells_L = replace(Cells_L, StartDate<SampleDateLocal & is.na(Cells_L), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
    summarise(Cells_L = sum(Cells_L)) %>% as.data.frame()     
  NRSSpecP1 <- rbind(NRSSpecP1, spec)
}

NRSSpecP1 <- NRSSpecP1 %>% group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(Cells_L = max(Cells_L)) %>% arrange(-desc(TaxonName)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSSpecP <-  NRSSpecP1 %>% pivot_wider(names_from = TaxonName, values_from = Cells_L, values_fill = list(Cells_L = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSSpecP, file = "NRS_phyto_species_mat.csv", row.names = FALSE)

#### NRS Zooplankton #### #################################################################################################################################
# Bring in all NRS phytoplankton samples
NRSZsamp <- read_csv("ZSampNRS.csv", na = "(null)") %>% 
  dplyr::rename("Sample" = "SAMPLE", "Station" = "STATION", "Latitude" = "LATITUDE", "Longitude" = "LONGITUDE", "SampleDateLocal" = "SAMPLEDATE", "NRScode" = "NRS_CODE") %>%
  mutate(Year = year(SampleDateLocal),
         Month = month(SampleDateLocal),
         Day = day(SampleDateLocal),
         Time_24hr = str_sub(SampleDateLocal, -8, -1)) # hms doesn't seem to work on 00:00:00 times

# Bring in plankton data
NRSZdat <- read_csv("NRS_zoop_raw.csv", na = "(null)") %>%
  dplyr::rename("Sample" = "SAMPLE", "TaxonName" = "TAXON_NAME", "Copepod" = "TAXON_GROUP", "TaxonGroup" = "TAXON_GRP01", 
                "Genus" = "GENUS", "Species" = "SPECIES", "ZAbund_m3" = "TAXON_PER_M3")

# Bring in Change Log
NRSZcl <- read_csv("ChangeLogNRSZ.csv", na = "(null)") %>%
  dplyr::rename("TaxonName" = "TAXON_NAME", "StartDate" = "START_DATE", "ParentName" = "PARENT_NAME")

#### NRS ZOOP RAW ####

NRSRawZ1 <- left_join(NRSZsamp, NRSZdat, by = "Sample") %>% select(c(1,3:11,16)) %>% arrange(-desc(TaxonName)) 
NRSRawZ <- NRSRawZ1 %>% pivot_wider(names_from = TaxonName, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSRawP, file = "NRS_zoop_raw_mat.csv", row.names = FALSE)

#### NRS ZOOP HTG ####

nrsHTGZ1 <- NRSZdat %>% group_by(Sample, TaxonGroup) %>% summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE)) %>%
  filter(!TaxonGroup %in% c('Other')) 
nrsHTGZ1 <-  NRSZsamp %>% left_join(nrsHTGZ1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), 'Copepod', TaxonGroup),
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3)) %>% arrange(-desc(TaxonGroup))
nrsHTGZ <-  nrsHTGZ1 %>% pivot_wider(names_from = TaxonGroup, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(nrsHTGZ[,-1], file = "NRS_zoop_HTG_mat.csv", row.names = FALSE)

#### NRS ZOOP GENUS ####

# Check genus are effected by change log
nrszlg <- NRSZcl %>% mutate(genus1 = word(TaxonName, 1),
                           genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species

NRSGenZ1 <- NRSZdat %>% filter(!TaxonName %in% levels(as.factor(nrszlg$TaxonName))) %>% group_by(Sample, Genus) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE)) %>% drop_na(Genus) 
NRSGenZ1 <- NRSZsamp %>% left_join(NRSGenZ1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), 'Acanthoica', word(Genus,1)), # bin subgenera together
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3))  %>% 
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
NRSGenP2 <- NRSZdat %>% filter(TaxonName %in% levels(as.factor(nrszlg$TaxonName))) %>% 
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE)) 

for (i in 1:nlevels(NRSGenP2$Genus)) {
  Dates <- as.data.frame(NRSGenP2) %>% filter(Genus == Genus[i]) %>% slice(1)  %>% droplevels()
  gen <- as.data.frame(NRSGenP2) %>% filter(Genus == Genus[i]) %>% droplevels() 
  gen <- NRSZsamp %>% left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()     
  NRSGenZ1 <- rbind(NRSGenZ1, gen)
}

NRSGenZ1 <- NRSGenZ1 %>% group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3)) %>% arrange(-desc(Genus)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSGenZ <-  NRSGenZ1 %>% pivot_wider(names_from = Genus, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSGenZ, file = "NRS_zoop_genus_mat.csv", row.names = FALSE)

#### NRS ZOOP COPEPODS ####

# Check at what level we need change log
nrsclc <- NRSZcl %>% mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

NRSCop1 <- NRSZdat %>% filter(!TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =='COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE))
NRSCop1 <- NRSZsamp %>% left_join(NRSCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), 'Calanus Australis', Species), # avoids nulls in pivot
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3))  %>%  # avoids nulls in pivot
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
NRSCop2 <- NRSZdat %>% filter(TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =='COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE)) 

for (i in 1:nlevels(NRSCop2$Species)) {
  Dates <- as.data.frame(NRSCop2) %>% filter(Species == Species[i]) %>% slice(1)  %>% droplevels()
  copes <- as.data.frame(NRSCop2) %>% filter(Species == Species[i]) %>% droplevels() 
  copes <- NRSZsamp %>% left_join(copes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()     
  NRSCop1 <- rbind(NRSCop1, copes)
}

NRSCop1 <- NRSCop1 %>% group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3)) %>% arrange(-desc(Species)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSCop <-  NRSCop1 %>% pivot_wider(names_from = Species, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSCop, file = "NRS_zoop_copes_mat.csv", row.names = FALSE)

#### NRS ZOOP NON-COPEPODS ####

# for non change log species

NRSnCop1 <- NRSZdat %>% filter(!TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =='NON-COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE))
NRSnCop1 <- NRSZsamp %>% left_join(NRSnCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), 'Calanus Australis', Species), # avoids nulls in pivot
         ZAbund_m3 = ifelse(is.na(ZAbund_m3), 0, ZAbund_m3))  %>%  # avoids nulls in pivot
  group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()

# add change log species with -999 for NA's and real absences as 0's
NRSnCop2 <- NRSZdat %>% filter(TaxonName %in% levels(as.factor(nrsclc$TaxonName)) & Copepod =='NON-COPEPOD'
                              & Species != 'spp.' & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(NRSZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE)) 

for (i in 1:nlevels(NRSnCop2$Species)) {
  Dates <- as.data.frame(NRSnCop2) %>% filter(Species == Species[i]) %>% slice(1)  %>% droplevels()
  ncopes <- as.data.frame(NRSnCop2) %>% filter(Species == Species[i]) %>% droplevels() 
  ncopes <- NRSZsamp %>% left_join(ncopes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbund_m3 = replace(ZAbund_m3, StartDate>SampleDateLocal, -999), 
           ZAbund_m3 = replace(ZAbund_m3, StartDate<SampleDateLocal & is.na(ZAbund_m3), 0)) %>% 
    group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3)) %>% as.data.frame()     
  NRSnCop1 <- rbind(NRSnCop1, ncopes)
}

NRSnCop1 <- NRSnCop1 %>% group_by(NRScode, Station, Latitude, Longitude, SampleDateLocal, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3)) %>% arrange(-desc(Species)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

NRSnCop <-  NRSnCop1 %>% pivot_wider(names_from = Species, values_from = ZAbund_m3, values_fill = list(ZAbund_m3 = 0)) %>% 
  arrange(desc(SampleDateLocal)) 

fwrite(NRSnCop, file = "NRS_zoop_noncopes_mat.csv", row.names = FALSE)

