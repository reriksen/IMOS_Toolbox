## IMOS plankton data products
## Claire Davies (CSIRO) and Jason D Everett (UQ/CSIRO)

## Created: May 2020
## Updated: 
## 21 July 2020 (Written to Git)
## 22 September 2020 (Updated data file structure)
## 13th November 2020 (Split CPR and NRS)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(data.table)
})

rawD <- "RawData"
outD <- "Output"

#### CPR Phytoplankton #######################################################################################################################################################
# Bring in all CPR phytoplankton samples
cprPsamp <- read_csv(paste0(rawD,.Platform$file.sep,"PSampCPR.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, Route = ROUTE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC) %>%
  select(-REGION) %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn"t seem to work on 00:00:00 times

# Bring in plankton data
cprPdat <- read_csv(paste0(rawD,.Platform$file.sep,"CPR_phyto_raw.csv"), na = "(null)") %>%
  rename(Sample = SAMPLE, TaxonName = TAXON_NAME, TaxonGroup = TAXON_GROUP, Genus = GENUS, Species = SPECIES, PAbun_m3 = PHYTO_ABUNDANCE_M3, BioVolume_um3_m3 = BIOVOL_UM3_M3)

# Bring in Change Log
cprPcl <- read_csv(paste0(rawD,.Platform$file.sep,"ChangeLogCPRP.csv"), na = "(null)") %>%
  rename(TaxonName = TAXON_NAME, StartDate = START_DATE, ParentName = PARENT_NAME)

#### CPR PHYTO RAW ####

cprRawP1 <- left_join(cprPsamp, cprPdat, by = "Sample") %>% 
  select(c(2:10,14)) %>% 
  arrange(-desc(TaxonName)) %>%
  mutate(TaxonName = ifelse(is.na(TaxonName), "No taxa found", TaxonName)) # for segments where no phyto was found

cprRawP <- cprRawP1 %>% 
  pivot_wider(names_from = TaxonName, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) %>%
  select(-"No taxa found")

fwrite(cprRawP, file = paste0(outD,.Platform$file.sep,"CPR_phyto_raw_mat.csv"), row.names = FALSE)

#### CPR PHYTO HTG ####

cprHTGP1 <- cprPdat %>% group_by(Sample, TaxonGroup) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop") %>%
  filter(!TaxonGroup %in% c("Other","Coccolithophore", "Diatom","Protozoa"))

cprHTGP1 <-  cprPsamp %>% left_join(cprHTGP1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), "Ciliate", TaxonGroup),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3)) %>% 
  arrange(-desc(TaxonGroup))

cprHTGP <-  cprHTGP1 %>% 
  pivot_wider(names_from = TaxonGroup, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) %>% 
  select(-Sample)

fwrite(cprHTGP, file = paste0(outD,.Platform$file.sep,"CPR_phyto_HTG_mat.csv"), row.names = FALSE)

#### CPR PHYTO GENUS ####

# Check genus are effected by change log
clg <- cprPcl %>% mutate(genus1 = word(TaxonName, 1),
                         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species

cprGenP1 <- cprPdat %>% 
  filter(!TaxonName %in% levels(as.factor(clg$TaxonName))) %>% 
  group_by(Sample, Genus) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop") %>% 
  drop_na(Genus)

cprGenP1 <- cprPsamp %>% 
  left_join(cprGenP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), "Acanthoica", Genus),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3)) %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprGenP2 <- cprPdat %>% 
  filter(TaxonName %in% levels(as.factor(clg$TaxonName))) %>% 
  left_join(cprPcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprGenP2$Genus)) {
  Dates <- as.data.frame(cprGenP2) %>% 
    filter(Genus == Genus[i]) %>% 
    slice(1)  %>% 
    droplevels()
  
  gen <- as.data.frame(cprGenP2) %>% 
    filter(Genus == Genus[i]) %>% 
    droplevels() 
  
  gen <- cprPsamp %>% 
    left_join(gen, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Genus = replace(Genus, is.na(Genus), Dates$Genus),
           PAbun_m3 = replace(PAbun_m3, StartDate>SampleDateUTC, -999), 
           PAbun_m3 = replace(PAbun_m3, StartDate<SampleDateUTC & is.na(PAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  
  cprGenP1 <- rbind(cprGenP1, gen)
}

cprGenP1 <- cprGenP1 %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(PAbun_m3 = max(PAbun_m3), .groups = "drop") %>% 
  arrange(-desc(Genus)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprGenP <-  cprGenP1 %>% 
  pivot_wider(names_from = Genus, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprGenP, file = paste0(outD,.Platform$file.sep,"CPR_phyto_genus_mat.csv"), row.names = FALSE)

#### CPR PHYTO SPECIES ####

# Check at what level we need change log
cls <- cprPcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species
cprSpecP1 <- cprPdat %>% 
  filter(!TaxonName %in% levels(as.factor(clg$TaxonName))
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species)) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop")

cprSpecP1 <- cprPsamp %>% 
  left_join(cprSpecP1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         TaxonName = ifelse(is.na(TaxonName), "Paralia sulcata", TaxonName),
         PAbun_m3 = ifelse(is.na(PAbun_m3), 0, PAbun_m3))  %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprSpecP2 <- cprPdat %>% 
  filter(TaxonName %in% levels(as.factor(clg$TaxonName))
         & Species != "spp." & !is.na(Species) 
         & !grepl("cf.", Species)) %>% 
  left_join(cprPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% 
  drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprSpecP2$TaxonName)) {
  Dates <- as.data.frame(cprSpecP2) %>% 
    filter(TaxonName == TaxonName[i]) %>% 
    slice(1) %>% 
    droplevels()
  
  spec <- as.data.frame(cprSpecP2) %>% 
    filter(TaxonName == TaxonName[i]) %>% 
    droplevels() 
  
  spec <- cprPsamp %>% 
    left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           PAbun_m3 = replace(PAbun_m3, StartDate>SampleDateUTC, -999), 
           PAbun_m3 = replace(PAbun_m3, StartDate<SampleDateUTC & is.na(PAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
    summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprSpecP1 <- rbind(cprSpecP1, spec)
}

cprSpecP1 <- cprSpecP1 %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(PAbun_m3 = max(PAbun_m3), .groups = "drop") %>% 
  arrange(-desc(TaxonName)) %>% 
  as.data.frame() 

# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified
cprSpecP <-  cprSpecP1 %>% 
  pivot_wider(names_from = TaxonName, values_from = PAbun_m3, values_fill = list(PAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprSpecP, file = paste0(outD,.Platform$file.sep,"CPR_phyto_species_mat.csv"), row.names = FALSE)

#### CPR Zoopplankton #### ################################################################################################################################
# Bring in all CPR zooplankton samples
cprZsamp <- read_csv(paste0(rawD,.Platform$file.sep,"ZSampCPR.csv"), na = "(null)") %>% 
  rename(Sample = SAMPLE, Route = ROUTE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC) %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn"t seem to work on 00:00:00 times

# Bring in plankton summary data
cprZdat <- read_csv(paste0(rawD,.Platform$file.sep,"CPR_zoo_raw.csv"), na = "(null)") %>%
  rename(Sample = SAMPLE, TaxonName = TAXON_NAME, Copepod = TAXON_GROUP, TaxonGroup = TAXON_GRP01,
         Genus = GENUS, Species = SPECIES, ZAbun_m3 = ZOOP_ABUNDANCE_M3)

# Bring in Change Log
cprZcl <- read_csv(paste0(rawD,.Platform$file.sep,"ChangeLogCPRZ.csv"), na = "(null)") %>%
  rename(TaxonName = TAXON_NAME, StartDate = START_DATE, ParentName = PARENT_NAME)

#### CPR ZOOP RAW ####
cprRawZ1 <- left_join(cprZsamp, cprZdat, by = "Sample") %>% 
  select(-c("Copepod", "TaxonGroup", "Genus", "Species")) %>% 
  arrange(-desc(TaxonName)) %>%
  mutate(TaxonName = ifelse(is.na(TaxonName), "No taxa found", TaxonName))

cprRawZ <- cprRawZ1 %>% 
  pivot_wider(names_from = TaxonName, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) %>%
  select(-"No taxa found") %>% 
  select(-Sample)

fwrite(cprRawZ, file = paste0(outD,.Platform$file.sep,"CPR_zoop_raw_mat.csv"), row.names = FALSE)

#### CPR ZOOP HTG ####
cprHTGZ1 <- cprZdat %>% 
  group_by(Sample, TaxonGroup) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop") %>%
  filter(!TaxonGroup %in% c("Other")) 

cprHTGZ1 <- cprZsamp %>% 
  left_join(cprHTGZ1, by = "Sample") %>% 
  mutate(TaxonGroup = ifelse(is.na(TaxonGroup), "Copepod", TaxonGroup),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3)) %>% 
  arrange(-desc(TaxonGroup))

cprHTGZ <-  cprHTGZ1 %>% 
  pivot_wider(names_from = TaxonGroup, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) %>% 
  select(-Sample)

fwrite(cprHTGZ, file = paste0(outD,.Platform$file.sep,"CPR_zoop_HTG_mat.csv"), row.names = FALSE)

#### CPR ZOOP GENUS ####

# Check genus that are affected by change log
clgz <- cprZcl %>% 
  mutate(genus1 = word(TaxonName, 1),
         genus2 = word(ParentName, 1)) %>%
  mutate(same = ifelse(genus1==genus2, "yes", "no")) %>%
  filter(same == "no")# no changes at genera level

# for non change log species
cprGenZ1 <- cprZdat %>% 
  filter(!TaxonName %in% levels(as.factor(clgz$TaxonName))) %>% 
  group_by(Sample, Genus) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop") %>% 
  drop_na(Genus) 

cprGenZ1 <- cprZsamp %>% 
  left_join(cprGenZ1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Genus = ifelse(is.na(Genus), "Calanus", Genus),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3)) %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprGenZ2 <- cprZdat %>% 
  filter(TaxonName %in% levels(as.factor(clgz$TaxonName))) %>% 
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Genus = as_factor(Genus)) %>% 
  drop_na(Genus) %>%
  group_by(Sample, StartDate, Genus) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprGenZ2$Genus)) {
  Datesz <- as.data.frame(cprGenZ2) %>% 
    filter(Genus == Genus[i]) %>% 
    slice(1) %>% 
    droplevels()
  
  genz <- as.data.frame(cprGenZ2) %>% 
    filter(Genus == Genus[i]) %>% 
    droplevels() 
  
  genz <- cprZsamp %>% 
    left_join(genz, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Datesz$StartDate),
           Genus = replace(Genus, is.na(Genus), Datesz$Genus),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprGenZ1 <- rbind(cprGenZ1, genz)
}

cprGenZ1 <- cprGenZ1 %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Genus) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3), .groups = "drop") %>% 
  arrange(-desc(Genus)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurrences as not regularly identified

cprGenZ <-  cprGenZ1 %>% 
  pivot_wider(names_from = Genus, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprGenZ, file = paste0(outD,.Platform$file.sep,"CPR_zoop_genus_mat.csv"), row.names = FALSE)

#### CPR ZOOP COPEPODS ####

# Check at what level we need change log
clc <- cprZcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

cprCop1 <- cprZdat %>% 
  filter(!TaxonName %in% levels(as.factor(clc$TaxonName)) & 
           Copepod =="COPEPOD" & 
           Species != "spp." & 
           !is.na(Species) & 
           !grepl("cf.", Species) & 
           !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop")

cprCop1 <- cprZsamp %>% 
  left_join(cprCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), "Calanus Australis", Species), # avoids nulls in pivot
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3)) %>% # avoids nulls in pivot
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprCop2 <- cprZdat %>% 
  filter(TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod =="COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>%
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprCop2$Species)) {
  Dates <- as.data.frame(cprCop2) %>% 
    filter(Species == Species[i]) %>% 
    slice(1) %>% 
    droplevels()
  
  copes <- as.data.frame(cprCop2) %>% 
    filter(Species == Species[i]) %>%
    droplevels() 
  
  copes <- cprZsamp %>% 
    left_join(copes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprCop1 <- rbind(cprCop1, copes)
}

cprCop1 <- cprCop1 %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3), .groups = "drop") %>% 
  arrange(-desc(Species)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurrences as not regularly identified

cprCop <- cprCop1 %>% 
  pivot_wider(names_from = Species, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

fwrite(cprCop, file = paste0(outD,.Platform$file.sep,"CPR_zoop_copes_mat.csv"), row.names = FALSE)

#### CPR ZOOP NON-COPEPODS ####

# for non change logspecies

cprnCop1 <- cprZdat %>% 
  filter(!TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod !="COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop")

cprnCop1 <- cprZsamp %>% 
  left_join(cprnCop1, by = "Sample") %>% 
  mutate(StartDate = ymd("2007-12-19"),
         Species = ifelse(is.na(Species), "Evadne spinifera", Species),
         ZAbun_m3 = ifelse(is.na(ZAbun_m3), 0, ZAbun_m3)) %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprnCop2 <- cprZdat %>% 
  filter(TaxonName %in% levels(as.factor(clc$TaxonName)) & Copepod !="COPEPOD"
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(Species = paste0(Genus," ", word(Species,1))) %>% # bin complexes
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% 
  drop_na(Species) %>% 
  group_by(Sample, StartDate, Species) %>% 
  summarise(ZAbun_m3 = sum(ZAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprnCop2$Species)) {
  Dates <- as.data.frame(cprnCop2) %>% 
    filter(Species == Species[i]) %>% 
    slice(1) %>% 
    droplevels()
  
  ncopes <- as.data.frame(cprnCop2) %>% 
    filter(Species == Species[i]) %>% 
    droplevels() 
  
  ncopes <- cprZsamp %>% 
    left_join(ncopes, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           ZAbun_m3 = replace(ZAbun_m3, StartDate>SampleDateUTC, -999), 
           ZAbun_m3 = replace(ZAbun_m3, StartDate<SampleDateUTC & is.na(ZAbun_m3), 0)) %>% 
    group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
    summarise(ZAbun_m3 = sum(ZAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprnCop1 <- rbind(cprnCop1, ncopes)
}

cprnCop1 <- cprnCop1 %>% 
  group_by(Route, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, Species) %>%
  summarise(ZAbun_m3 = max(ZAbun_m3), .groups = "drop") %>% 
  arrange(-desc(Species)) %>% as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

cprnCop <-  cprnCop1 %>% 
  pivot_wider(names_from = Species, values_from = ZAbun_m3, values_fill = list(ZAbun_m3 = 0)) %>% 
  arrange(desc(SampleDateUTC)) 

# fwrite(cprnCop, file = paste0(outD,.Platform$file.sep,"CPR_zoop_noncopes_mat.csv"), row.names = FALSE)

