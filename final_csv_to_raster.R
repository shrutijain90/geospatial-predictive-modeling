library(raster)
library(rgdal)

setwd("C:/Users/Shruti Jain/Box Sync/Laptop data/Downloaded/Fraym-Shruti/")

kenya_pred <- read.csv("data/kenya_pred.csv", header = TRUE)

#kenya boundary
boundary <- readOGR("data/boundaries/KEN_adm0.shp")

x <- raster(extent(boundary), res = c(0.008332999999999997,0.008332999999999997), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

kenya_raster <- rasterize(kenya_pred[, c('long','lat')], x, kenya_pred[, 'y_pred'])

writeRaster(kenya_raster, filename="data/kenya_raster.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
