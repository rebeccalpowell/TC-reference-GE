##load in necessary libraries
library(sp)
library(rgdal)
library(tiff)
library(raster)
library(maptools)
library(rgeos)

##set working directory to validation data file
setwd("S:/POWELL/Research_current/Serengeti/CalibrationValidation/")

##load in image scene
modis_test <- raster("S:/POWELL/Research_current/Serengeti/CalibrationValidation/Data/Modis250m.sur_refl_b01.tif")
proj4string(modis_test) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")
##load in tree gradient by which the samples will be stratified
gradient <- readOGR(dsn = "Data_Creation_Package/Serengeti_3Class" , 
                    layer = "Serengeti_3Class")
gradient <- spTransform(gradient, CRS("+proj=longlat +ellps=WGS84 +no_defs"))

##Query each gradient
dense <- gradient[gradient$Class =="Forest",]
sparse <- gradient[gradient$Class =="Sparse",]

##crop raster to AOI to reduce processing time
modis_crop <- crop(modis_test, gradient, snap='in')
##convert raster to polygons
ply_crop <- rasterToPolygons(modis_crop, n=4)

##spatial sample, stratified by gradient area
set.seed(35)
##will repeat for each area to achieve stratified sampling schema
sample_dense <- spsample(dense, 150, type="random")
sample_sparse <- spsample(sparse, 150, type="random")

##extract polygoin OID by points, repeat for each stratification and combine
set1 <- over(sample_dense, ply_crop, returnList=TRUE)
set2 <- over(sample_sparse, ply_crop, returnList=TRUE)

##subset by queary row name/ID
##create blank spatial polygons data frame

##Gradient class 1
plygns1 <- ply_crop[0,]
##loops through all sampled polygons
for(i in 1:length(set1)){
  x <- rownames(set1[[i]])
  ply <- (ply_crop[x,])
  plygns1 <- spRbind(plygns1, ply)
}

##add column to describe strata
plygns1@data$Class <- "Dense"

##gradient class 2
plygns2 <- ply_crop[0,]
##loops through all sampled polygons
for(i in 1:length(set2)){
  x <- rownames(set2[[i]])
  ply <- (ply_crop[x,])
  plygns2 <- spRbind(plygns2, ply)
}
plygns2@data$Class <- "Sparse"

##combine all gradient samples into one variable
final <- spRbind(plygns1, plygns2) 

final$names <- seq(0:101)

##write out shapefiles of selected polygons
writeOGR(final, "Fir_Shps_2","pltplygns",driver="ESRI Shapefile")
