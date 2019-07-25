# Netcdf import of CSIRO CARS climatology data
# Claire
# Making CARS layers and selecting data

library(raster)
library(ncdf4)
library(raster)
library(ncdf4)
library(Matrix)
library(plyr)
library(nlme)
library(mgcv)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(grid)
library(effects)
library(splines)
library(pracma)
library(fields)

###############################
# GET MONTHLY SST DATA INTO R #
###############################
rm(list = ls())
setwd("H://R/data")

# Open netCDF in R
nc <- nc_open("~/R/Data/oxygen_cars2009.nc") # Change filename to any CARS2009 product
#nc <- nc_open("//argos-hba.it.csiro.au/sdode-data/Product407/A20170012017031.L3m_MO_KD490_Kd_490_4km.nc")
nc # Calls the netCDF, which lists all of the variables
dname <- "mean"
# Get lon, lat, and sea_surface_temperature
Lons <- ncvar_get(nc,"lon")
Lats <- ncvar_get(nc,"lat")
Depths <- ncvar_get(nc,"depth")

tmp <- ncvar_get(nc,dname)
dlname <- ncatt_get(nc,dname,"long_name")
dunits <- ncatt_get(nc,dname,"units")
fillvalue <- ncatt_get(nc,dname,"_FillValue")
dim(tmp) # Lon, Lat, Depth

# Which depth? Let's just take the surface! 0 m
tmp <- apply(t(tmp[ , ,1]), 2, rev) # transpose and then flip data for plotting

# Plot mean field
r <- raster(nrow = length(Lats), ncol = length(Lons)) # Make a global raster
r[] <- tmp # Write to raster
cols <- rev(rainbow(20))
plot(r, col = cols) # Check it works!
# Mean field

# To get value for a particular Day of Yr
# Evaluate at day-of-year 45 (mid February)
# CARS stores the mean value and the annual and semi-annual sin and cos coefficients for each lat, lon, depth
# tmp_values <- mean + an_cos*cos(t) + an_sin*sin(t) + sa_cos*cos(2*t) + sa_sin*sin(2*t)
an_cos <- ncvar_get(nc,"an_cos")
an_sin <- ncvar_get(nc,"an_sin")
sa_cos <- ncvar_get(nc,"sa_cos")
sa_sin <- ncvar_get(nc,"sa_sin")

# transpose and then flip data for plotting
an_cos <- apply(t(an_cos[ , , 1]), 2, rev)
an_sin <- apply(t(an_sin[ , , 1]), 2, rev)
sa_cos <- apply(t(sa_cos[ , , 1]), 2, rev)
sa_sin <- apply(t(sa_sin[ , , 1]), 2, rev)

DOY <- 45
t <- 2 * pi * DOY/366
tmp_values <- tmp + an_cos*cos(t) + an_sin*sin(t) + sa_cos*cos(2*t) + sa_sin*sin(2*t)

# Plot variable for particular DOY
#r <- raster(nrow = length(Lats), ncol = length(Lons)) # Make a global raster
r <- raster(nrow = length(Lats), ncol = length(Lons), xmn=0, xmx=360) # Make a global raster
r[] <- tmp_values # Write to raster
plot(r, col = cols) # Check it works!

##################################
#make layer for CARS ncdf
##################################

#get data from satelites
lat <- read.csv("~/R/Satellite_CSV/lat.csv")
long <- read.csv("~/R/Satellite_CSV/lon1.csv")
Lat <- data.matrix(rbind(-8.0208, lat))
long <- data.matrix(long)
nLon <- nrow(long)
nLat <- nrow(Lat)

dim(tmp_values)
fliptv <- tmp_values[c(331:1),] #as lats need to be increasing in interp2

mesh <- meshgrid(long, Lat)
var2 <- interp2(Lons, Lats, fliptv, mesh$X, mesh$Y, method = "linear")

var2m <- matrix(var2, nrow = length(Lat), ncol = length(long))

write.matrix(var2m,"Ox.prn")

MapCoords <- c(105, 178, -55, -8) # LonMin, LonMax, LatMin, LatMax
e <- extent(MapCoords)
ext <- as.vector(e) 
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

var2Layer <- raster(var2m, xmn = e[1], xmx = e[2], ymn = e[3], ymx = e[4])

#to plot var2layer
x11(width = 10, height = 5)
boundaries <- map('worldHires',
                  xlim=ext[1:2], ylim=ext[3:4],
                  plot=FALSE)
boundaries <- map2SpatialLines(boundaries,
                               proj4string=CRS(projection(var2Layer)))
spplot(var2Layer, col.regions=jet.colors,
       sp.layout=list('sp.lines', boundaries, lwd=0.5),
       colorkey=list(height=0.5),
       scales = list(draw = TRUE), xlab = "Longitude", ylab = "Latitude")



##################################
#select value from CARS netcdf for each location/time to the nearest 0.5 degrees
#################################
nc <- nc_open("H://R/Data/phosphate_cars2009.nc")
dat <- read.csv("H://R/data/eReefs_data.csv")
nz <- length(dat$LATITUDE)
val_id <- matrix(0, nrow=nz, ncol=5)
colnames(val_id)=c("lat", "long", "DOY", "mean","Value")

for (i in 1:nz) {
LonIdx <- which( nc$dim$lon$vals == round((dat$LONGITUDE/0.5)*0.5)[i])
LatIdx <- which( nc$dim$lat$vals == round((dat$LATITUDE/0.5)*0.5)[i])
DepIdx <- which( nc$dim$depth$vals == 0)
meani <- ncvar_get(nc, "mean", start=c(LonIdx, LatIdx, DepIdx), count = c(1,1,1))
an_cosi <- ncvar_get(nc,"an_cos", start=c(LonIdx, LatIdx, DepIdx), count = c(1,1,1))
an_sini <- ncvar_get(nc,"an_sin", start=c(LonIdx, LatIdx, DepIdx), count = c(1,1,1))
sa_cosi <- ncvar_get(nc,"sa_cos", start=c(LonIdx, LatIdx, DepIdx), count = c(1,1,1))
sa_sini <- ncvar_get(nc,"sa_sin", start=c(LonIdx, LatIdx, DepIdx), count = c(1,1,1))
DOY <- dat$DOY[i]
t <- 2 * pi * DOY/366
valuei <- meani + an_cosi*cos(t) + an_sini*sin(t) + sa_cosi*cos(2*t) + sa_sini*sin(2*t)
val_id[i] <- as.numeric(dat$LONGITUDE[i])
val_id[i,2] <- as.numeric(dat$LATITUDE[i])
val_id[i,3] <- as.numeric(dat$DOY[i])
val_id[i,4] <- as.numeric(mean)
val_id[i,5] <- as.numeric(valuei)
}

###################################
#as above to the interpolated value 
###################################
dat <- read.csv("H://R/data/eReefs_data.csv", header= TRUE, na.strings="(null)")
dat$date <- as.POSIXct(strptime(dat$date, format="%Y-%m-%d", tz="Etc/GMT-10"))
dat$DOY <- as.numeric(format(dat$date, "%j"))
nz <- length(dat$latitude)

nc <- nc_open("H://R/data/phosphate_cars2009.nc")

for (i in 1:nz) {
  x2 <- (ceiling(dat$longitude/0.5)*0.5)[i] #max long
  x1 <- (floor(dat$longitude/0.5)*0.5)[i]   #min long
  y1 <- (floor(dat$latitude/0.5)*0.5)[i]    #min lat
  y2 <- (ceiling(dat$latitude/0.5)*0.5)[i]  #max lat  
  x <- dat$longitude[i]
  y <- dat$latitude[i]
  LonIdxc <- which( nc$dim$lon$vals == x2)
  LonIdxf <- which( nc$dim$lon$vals == x1)
  LatIdxc <- which( nc$dim$lat$vals == y1)
  LatIdxf <- which( nc$dim$lat$vals == y2)
  DepIdx <- which( nc$dim$depth$vals == 0)
  
  meanicc <- ncvar_get(nc, "mean", start=c(LonIdxc, LatIdxc, DepIdx), count = c(1,1,1))
  an_cosicc <- ncvar_get(nc,"an_cos", start=c(LonIdxc, LatIdxc, DepIdx), count = c(1,1,1))
  an_sinicc <- ncvar_get(nc,"an_sin", start=c(LonIdxc, LatIdxc, DepIdx), count = c(1,1,1))
  sa_cosicc <- ncvar_get(nc,"sa_cos", start=c(LonIdxc, LatIdxc, DepIdx), count = c(1,1,1))
  sa_sinicc <- ncvar_get(nc,"sa_sin", start=c(LonIdxc, LatIdxc, DepIdx), count = c(1,1,1))
  DOY <- dat$DOY[i]
  t <- 2 * pi * DOY/366
  valueicc <- meanicc + an_cosicc*cos(t) + an_sinicc*sin(t) + sa_cosicc*cos(2*t) + sa_sinicc*sin(2*t)
  
  meaniff <- ncvar_get(nc, "mean", start=c(LonIdxf, LatIdxf, DepIdx), count = c(1,1,1))
  an_cosiff <- ncvar_get(nc,"an_cos", start=c(LonIdxf, LatIdxf, DepIdx), count = c(1,1,1))
  an_siniff <- ncvar_get(nc,"an_sin", start=c(LonIdxf, LatIdxf, DepIdx), count = c(1,1,1))
  sa_cosiff <- ncvar_get(nc,"sa_cos", start=c(LonIdxf, LatIdxf, DepIdx), count = c(1,1,1))
  sa_siniff <- ncvar_get(nc,"sa_sin", start=c(LonIdxf, LatIdxf, DepIdx), count = c(1,1,1))
  valueiff <- meaniff + an_cosiff*cos(t) + an_siniff*sin(t) + sa_cosiff*cos(2*t) + sa_siniff*sin(2*t)
  
  meanifc <- ncvar_get(nc, "mean", start=c(LonIdxf, LatIdxc, DepIdx), count = c(1,1,1))
  an_cosifc <- ncvar_get(nc,"an_cos", start=c(LonIdxf, LatIdxc, DepIdx), count = c(1,1,1))
  an_sinifc <- ncvar_get(nc,"an_sin", start=c(LonIdxf, LatIdxc, DepIdx), count = c(1,1,1))
  sa_cosifc <- ncvar_get(nc,"sa_cos", start=c(LonIdxf, LatIdxc, DepIdx), count = c(1,1,1))
  sa_sinifc <- ncvar_get(nc,"sa_sin", start=c(LonIdxf, LatIdxc, DepIdx), count = c(1,1,1))
  valueifc <- meanifc + an_cosifc*cos(t) + an_sinifc*sin(t) + sa_cosifc*cos(2*t) + sa_sinifc*sin(2*t)

  meanicf <- ncvar_get(nc, "mean", start=c(LonIdxc, LatIdxf, DepIdx), count = c(1,1,1))
  an_cosicf <- ncvar_get(nc,"an_cos", start=c(LonIdxc, LatIdxf, DepIdx), count = c(1,1,1))
  an_sinicf <- ncvar_get(nc,"an_sin", start=c(LonIdxc, LatIdxf, DepIdx), count = c(1,1,1))
  sa_cosicf <- ncvar_get(nc,"sa_cos", start=c(LonIdxc, LatIdxf, DepIdx), count = c(1,1,1))
  sa_sinicf <- ncvar_get(nc,"sa_sin", start=c(LonIdxc, LatIdxf, DepIdx), count = c(1,1,1))
  valueicf <- meanicf + an_cosicf*cos(t) + an_sinicf*sin(t) + sa_cosicf*cos(2*t) + sa_sinicf*sin(2*t)
  
  v1 <- ((x2-x)/(x2-x1))*valueifc + ((x-x1)/(x2-x1)*valueicc)
  v2 <- ((x2-x)/(x2-x1))*valueiff + ((x-x1)/(x2-x1)*valueicf)
  vxy <- ((y2-y)/(y2-y1))*v1 + ((y-y1)/(y2-y1)*v2)
  
  dat[i,12] <- as.numeric(vxy)

}

write.csv(dat, "H://r/data/eReefs_data_carsP.csv", row.names = FALSE)


##################################
#TURBIDITY LAYER FROM NETCDF
##################################
nc <- nc_open("//argos-hba.it.csiro.au/sdode-data/Product407/A20171822017212.L3m_MO_KD490_Kd_490_4km.nc")
nc # Calls the netCDF, which lists all of the variables
dname <- "Kd_490"
Lons <- ncvar_get(nc,"lon")
Lats <- ncvar_get(nc,"lat")

LonIdx <- which(nc$dim$lon$vals > 105 & nc$dim$lon$vals < 177.95 )
LatIdx <- which( nc$dim$lat$vals > -55 & nc$dim$lat$vals < -8 )

tmp <- ncvar_get(nc,dname)[LonIdx, LatIdx]
dlname <- ncatt_get(nc,dname,"long_name")
dunits <- ncatt_get(nc,dname,"units")
fillvalue <- ncatt_get(nc,dname,"_FillValue")
dim(tmp) # Lon, Lat, Depth

tmp <- apply(t(tmp[ , ]), 2, rev) # if no depth profiles
dim(tmp)
tmp <- tmp[c(1128:1),] 

# Plot mean field
e <- c(105, 178, -55, -8) # LonMin, LonMax, LatMin, LatMax
r <- raster(tmp, xmn = e[1], xmx = e[2], ymn = e[3], ymx = e[4]) # Make a raster
r[] <- tmp # Write to raster
cols <- rev(rainbow(20))
plot(r, col = cols) # Check it works!

library(MASS)
write.matrix(tmp, "Turb7.prn")

##################################
#ANG LAYER FROM NETCDF
##################################
nc <- nc_open("~/R/data/A20160012016031.L3m_MO_RRS_angstrom_9km.nc")
nc # Calls the netCDF, which lists all of the variables
dname <- "angstrom"
Lons <- ncvar_get(nc,"lon")
Lats <- ncvar_get(nc,"lat")

LonIdx <- which(nc$dim$lon$vals > 105 & nc$dim$lon$vals < 177.95 )
LatIdx <- which( nc$dim$lat$vals > -55 & nc$dim$lat$vals < -8 )
 
Lons <- ncvar_get(nc,"lon")[LonIdx]
Lats <- ncvar_get(nc,"lat")[LatIdx]
Lats <- Lats[c(564:1)]

tmp <- ncvar_get(nc,dname)[LonIdx, LatIdx]
dlname <- ncatt_get(nc,dname,"long_name")
dunits <- ncatt_get(nc,dname,"units")
fillvalue <- ncatt_get(nc,dname,"_FillValue")
dim(tmp) # Lon, Lat, Depth

tmp <- apply(t(tmp[ , ]), 2, rev) # if no depth profiles
dim(tmp)

fliptv <- tmp[c(564:1),] #as lats need to be increasing in interp2

#get data from satelites
lat <- read.csv("~/R/Satellite_CSV/lat.csv")
long <- read.csv("~/R/Satellite_CSV/lon1.csv")
Lat <- data.matrix(rbind(-8.0208, lat))
long <- data.matrix(long)
nLon <- nrow(long)
nLat <- nrow(Lat)

mesh <- meshgrid(long, Lat)
var2 <- interp2(Lons, Lats, fliptv, mesh$X, mesh$Y, method = "linear")

var2m <- matrix(var2, nrow = length(Lat), ncol = length(long))
dim(var2m)
var2m <- var2m[c(1128:1),]
write.matrix(var2m,"ang.prn")

# Plot mean field
r <- raster(var2m, xmn = e[1], xmx = e[2], ymn = e[3], ymx = e[4]) # Make a raster
r[] <- var2m # Write to raster
cols <- rev(rainbow(20))
plot(r, col = cols) # Check it works!




