fNetcdf_2_Raster <- function(nc, v, x, y) {
  
  library(raster)
  library(ncdf4)
  
  # Extract data from the netCDF file  
  nc <- nc_open(nc)
  dat <- ncvar_get(nc, v) # x, y, year 
  dat[] <- dat
  X <- dim(dat)[1]
  Y <- dim(dat)[2]
  tt <- nc.get.time.series(nc, v = "time", time.dim.name = "time") # from packages ncdf4.helpers&PCICt
  tt <- as.POSIXct(tt)
  tt <- as.Date(tt)
  nc_close(nc)
  rs <- raster(nrow = Y, ncol = X) # Make a raster with the right dims to fill with lat&lon
  
  # Fix orientation of original data
  drs <- data.frame(coordinates(rs))
  
  # Create rasters stacks
  rs_list <- list() # empty list to allocate results
  st <- stack()
  for (i in 1:length(tt)) {
    dt1 <- rasterFromXYZ(cbind(drs, as.vector(dat[,, i])))
    dt1[]<- ifelse(dt1[] <= -2, NA, dt1[]) # for some models that have weird temperatures (-273)
    dt1[]<- ifelse(dt1[] >= 40, NA, dt1[])
    st <- addLayer(st, flip(dt1, 2))
    print(paste0(i, " of ", length(tt)))
  }
  names(st) <- seq(as.Date(paste(from, "1", "1", sep = "/")), as.Date(paste(to, "12", "1", sep = "/")), by = "year")
  crs(st) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  return(st)
}