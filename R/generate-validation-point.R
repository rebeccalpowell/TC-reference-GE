##-------------------Script to Interactively Create Tree Cover Validation Data Using Google Earth Imagery -----------------##

##THIS SCRIPT CONVERTS MODIS PIXEL SIZED PLOTS TO 90m LANDSAT PLOTS BY + or - 0.00075 TO THE BOUNDING BOX
##INPUT DATA SHOULD BE 250m PLOTS OR SCRIPT MUST BE MODIFIED.

##load in relevant libraries
library(sp)
library(rgdal)
library(tiff)
library(raster)
library(maptools)
library(RgoogleMaps)
library(dismo)

##set working directory
setwd("Y:/Research_current/Serengeti/CalibrationValidation/Data/")
plots <- readOGR("ValidationClusterAddPlotsMODIS", "pltplygns")
##it is important the coordiantes are in lat long in order to access GE imagery
plots <- spTransform(plots, CRS("+proj=longlat +ellps=WGS84 +no_defs"))

##this code can be used to take on the plots in chunks instead of all at once- once you start a loop you cannot safely exit without finishing
##plots <- plots[1:200,]
##plots <- plots[201:400,]
##plots <- plots[401:600,]
##plots <- plots[601:nrow(plots),]

plots@data$LandsatTC <- NA
plots.df <- as.data.frame(plots)
plots.df$LandsatTC <- NA

##--------------------- this code is needed to initiate a polygon file -----------------
i <- 1
##load in relevant shapefiles or create coordinates of bounding box
lonmin <- extent(plots[i,])@xmin + 0.00075
lonmax <- extent(plots[i,])@xmax - 0.00075
latmin <- extent(plots[i,])@ymin + 0.00075
latmax <- extent(plots[i,])@ymax - 0.00075

##generates bounding box
mat <- matrix(c(lonmin, lonmin, lonmax, lonmax, latmax, latmin, latmin, latmax), nrow = 4, ncol = 2)
box <- qbbox(lat=mat[,2], lon = mat[,1])
p <- Polygon(mat)
ps <- Polygons(list(p),1)
x <- SpatialPolygons(list(ps))
proj4string(x) <- CRS("+proj=longlat +datum=WGS84")

##----------------- this part creates IDs if wanted/needed ----------------------
##Creates a list of ID values for data set
##IDs <- sapply(slot(plots, "polygons"), function(i) slot(i, "ID"))
##In order to use IDs, a line must be added where i <- ID[x], and x is used as the for loop iterator.

##------------------ The meat of the code ---------------
##select number of samples
##creates a spatial grid of sample points (9*9=81, in a 90M^2 site samples are ~10m apart)
n <- 100

##code to make a for loop through all data points
for (i in 1:nrow(plots.df)) {
  ##i <- 156
  print(i)
  ##load in relevant shapefiles or create coordinates of bounding box
  lonmin <- extent(plots[i,])@xmin + 0.00075
  lonmax <- extent(plots[i,])@xmax - 0.00075
  latmin <- extent(plots[i,])@ymin + 0.00075
  latmax <- extent(plots[i,])@ymax - 0.00075

  ##generates bounding box
  mat <- matrix(c(lonmin, lonmin, lonmax, lonmax, latmax, latmin, latmin, latmax), nrow = 4, ncol = 2)
  box <- qbbox(lat=mat[,2], lon = mat[,1])
  p <- Polygon(mat)
  ps <- Polygons(list(p),1)
  boundBox <- SpatialPolygons(list(ps))
  proj4string(boundBox) <- CRS("+proj=longlat +datum=WGS84")

  ##finds range of bounding box- used to determine the best available resolution
  lats = c(latmax, latmin)
  lons = c(lonmin, lonmax)

  ##Rgooglemaps function that determines finest available resolution
  zoomFact = min(MaxZoom(range(lats),range(lons)))

  ##acquire static google earth image
  basemap <- gmap(boundBox, lonlat=TRUE, zoom = zoomFact, type = "satellite")

  ##creates sample points along a regular grid, n= determines the number of points included in bounding box
  set.seed(42)
  samplegrid <- spsample(boundBox, n, type= "regular", bb=bbox(boundBox), offset= c(0.5, 0.5))

  ##plot object boundary
  plot(basemap)
  plot(samplegrid, add=TRUE, col="red", pch=1, cex=.85, lwd=2)
  plot(boundBox, add=TRUE, lwd=2)

  ##activate interactive selection, use esc to finish selection
  trees <- identify(samplegrid$x1, samplegrid$x2, labels=row.names(samplegrid))

  ##create independent variable of selected points, record number of points selected
  selected <- samplegrid[trees,]

  ##plot selected points to confirm they are "tree"
  plot(selected, add=TRUE, pch=18, col="blue")

  ##output for number of points included as trees
  plots.df[i, "LandsatTC"] <- 0 ##in case no points are selected
  plots.df[i, "LandsatTC"] <- (summary(selected)$npoints/n)*100

  plots@data[i, "LandsatTC"] <- 0 ##in case no points are selected
  plots@data[i, "LandsatTC"] <- (summary(selected)$npoints/n)*100

  boundBox <- spChFIDs(boundBox, paste("D", i, sep = ""))
  x <- spRbind(x, boundBox)
  print((summary(selected)$npoints/n)*100)

}


##---------------- ONCE YOU ARE FINISHED WITH ALL RESULTS ----------------##


plots@data[101, 1] <- NA


shapebind <- x[paste("D", seq(from=1, to=nrow(plots@data), by =1), sep = ""),]
row.names(plots.df) <- paste("D", seq(from=1, to=nrow(plots@data), by =1), sep = "")

shapebind <- SpatialPolygonsDataFrame(shapebind,data=as.data.frame(plots.df))

writeOGR(shapebind, "TC_VlClstrAddPlts","pltplygns", driver="ESRI Shapefile")

##write out results to CSV
write.csv(plots.df, file="S:/POWELL/Research_current/Serengeti/Validation_Data/landsat4.csv")

