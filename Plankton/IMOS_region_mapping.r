## Map segment to IMCRA and marine bioregions
## Claire Davies (CSIRO) & Jason Everett (CSIRO / UQ)
## Created 6 Nov 2020

library(sp)
library(maptools)

## Read the IMCRA mesoscale shapefile.

imcra <-
    readShapePoly("shape/imcra4_meso",
                  proj4string=CRS("+proj=longlat +datum=WGS84"))

imcra_det <- imcra@data[,1:3]
write.csv(imcra_det, "imcra_details.csv")

## Read the segments.

segments <- read_csv("segments.csv") # file with columns named Longitude, Latitude

## Convert the segments table to a SpatialPolygonsDataFrame.

coordinates(segments) <- c("Longitude", "Latitude")
proj4string(segments) <- CRS("+proj=longlat +datum=WGS84")

## Perform the point in polygon overlay.

segImcra <- over(segments, imcra)

## Join the segment points to the IMCRA attributes (just by row order because
## the row orders of segments and segImcra are guaranteed to match).

segImcra <- spCbind(segments, segImcra)

## Plot a map for checking.

mesoColor <- rev(rainbow(nlevels(imcra$MESO_ABBR), start=0, end=5 / 6))
plot(imcra, col=mesoColor[imcra$MESO_ABBR])
points(segImcra, pch=19, cex=0.5, col=mesoColor[segImcra$MESO_ABBR])
points(segImcra, pch=2, col="Black")

## Export the results for Oracle.

outTable <- segImcra@data[c("SILK_ID", "SEGMENT_NO" , "MESO_ABBR")]
#write_csv(outTable, "segment_imcra.csv")

########################################
## Read the IMCRA bioregions shapefile.

imcrapb <-
  readShapePoly("shape/imcra4_pb",
                proj4string=CRS("+proj=longlat +datum=WGS84"))

## Read the segments.

segments <- read.csv("segments.csv") # file with columns named Longitude, Latitude

## Convert the segments table to a SpatialPolygonsDataFrame.

coordinates(segments) <- c("Longitude", "Latitude")
proj4string(segments) <- CRS("+proj=longlat +datum=WGS84")

## Perform the point in polygon overlay.

segImcrapb <- over(segments, imcrapb)

## Join the segment points to the IMCRA attributes (just by row order because
## the row orders of segments and segImcra are guaranteed to match).

segImcrapb <- spCbind(segments, segImcrapb)

## Plot a map for checking.

mesoColor <- rev(rainbow(nlevels(imcrapb$PB_NAME), start=0, end=5 / 6))
plot(imcrapb, col=mesoColor[imcrapb$PB_NAME])
points(segImcrapb, pch=19, cex=0.5, col=mesoColor[segImcrapb$PB_NAME])
points(segImcrapb, pch=20, col="Black")

## Export the results for Oracle.

outTable <- segImcrapb@data[c("SILK_ID","SEGMENT_NO",  "PB_NAME")]
#write_csv(outTable,"segment_imcraPB.csv")


########################################################
## Read the Marine planning bioregions shapefile.

mbr <-
  readShapePoly("shape/marine_regions",
                proj4string=CRS("+proj=longlat +datum=WGS84"))

## Read the segments.

segments <- read.csv("segments.csv") # file with columns named Longitude, Latitude
 
## Convert the segments table to a SpatialPolygonsDataFrame.

coordinates(segments) <- c("Longitude", "Latitude")
proj4string(segments) <- CRS("+proj=longlat +datum=WGS84")

## Perform the point in polygon overlay.

segmbr <- over(segments, mbr)

## Join the segment points to the IMCRA attributes (just by row order because
## the row orders of segments and segImcra are guaranteed to match).

segmbr <- spCbind(segments, segmbr)

## Plot a map for checking.

mesoColor <- rev(rainbow(nlevels(mbr$REGION), start=0, end=5 / 6))
plot(mbr, col=mesoColor[mbr$REGION])
points(segmbr, pch=19, cex=0.5, col=mesoColor[segmbr$REGION])
points(segmbr, pch=20, col="Black")

## Export the results for Oracle.

outTablembr <- segmbr@data[c("SILK_ID", "SEGMENT_NO", "REGION")]
#write.csv(outTablembr, file="segment_mbr.csv", na="", row.names=FALSE)

###########################

