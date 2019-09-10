library(rgeos)
library(raster)
library(rgdal)
library(gdalUtils)
library(MODIS)

setwd("C:/Users/Shruti Jain/Box Sync/Laptop data/Downloaded/Fraym-Shruti")

# data file with y labels
cust <- read.csv("data/customer_profiles_final.csv")

# kenya_centroid
lon <- 37.9062
lat <- -0.0236

#kenya boundary
boundary <- readOGR("data/boundaries/KEN_adm0.shp")
buffer <- raster::buffer(x=boundary, width=0.001)

xy <- cust[,c(4,3)]
spdf <- SpatialPointsDataFrame(coords = xy, data = cust,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
newcrs <- CRS(paste0("+proj=laea + ellps=WGS84 + datum=WGS84 + x_0=0 +y_0=0 +lon_0=", lon, "+lat_0=", lat))
sp_cust <- spTransform(spdf, newcrs)
sp_cust_buffer <- gBuffer(sp_cust, byid = TRUE, width = 500, quadsegs=200)


####### spatial data inputs

#### population, literacy, poverty, and nightlights

## population
pop <- raster("data/population/KEN_popmap15adj_v2b.tif")
pop <- projectRaster(pop, crs=newcrs)
pop[pop<0] <- 0

## literacy
quads <- list.files("data/populationfemale09", pattern = "\\.tif$")
popf09 <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/populationfemale09/",quads[i]))
  popf09 <- stack(popf09,r)
}
popf09 <- stackApply(popf09, indices =  rep(1,nlayers(popf09)), fun = "sum", na.rm = T)
popf09 <- projectRaster(popf09, crs=newcrs)
popf09[popf09<0] <- 0

lit <- raster("data/literacy/KEN_literacy_F.tif") # need female population 09 to aggregate
lit <- projectRaster(lit, crs=newcrs)
lit[lit<0] <- 0
lit[lit>1] <- 1
lit <- projectRaster(lit, popf09)

## poverty
pop09 <- raster("data/population09/ken_ppp_2009.tif")
pop09 <- projectRaster(pop09, crs=newcrs)
pop09[pop09<0] <- 0

pov <- raster("data/poverty/ken08povmpi.tif")  # need population 09 to aggregate
pov <- projectRaster(pov, crs=newcrs)
pov[pov<0] <- 0
pov[pov>1] <- 1
pov <- projectRaster(pov, pop09)

# For weighted sum and aggregated proportions
cust$population = NaN
cust$literacy = NaN
cust$poverty = NaN
for (row in 1:nrow(cust)) {
  df_pop = as.data.frame(extract(pop, sp_cust_buffer[row,], weights=T, normalizeWeights=F))
  cust$population[row] =  sum(df_pop$value * df_pop$weight, na.rm=T)
  
  df_popf09 = as.data.frame(extract(popf09, sp_cust_buffer[row,], weights=T, normalizeWeights=F))
  df_lit = as.data.frame(extract(lit, sp_cust_buffer[row,], weights=T, normalizeWeights=F))
  cust$literacy[row] =  sum(df_popf09$value * df_lit$value * df_lit$weight, na.rm=T) / sum(df_popf09$value * df_lit$weight, na.rm=T)
  
  df_pop09 = as.data.frame(extract(pop09, sp_cust_buffer[row,], weights=T, normalizeWeights=F))
  df_pov = as.data.frame(extract(pov, sp_cust_buffer[row,], weights=T, normalizeWeights=F))
  cust$poverty[row] =  sum(df_pop09$value * df_pov$value * df_pov$weight, na.rm=T) / sum(df_pop09$value * df_pov$weight, na.rm=T)
}

## nightlights
lights <- raster("data/lights/F182013.v4c_web.stable_lights.avg_vis.tif")
lights <- crop(lights, extent(buffer))
lights <- mask(x=lights, mask=buffer)
lights <- projectRaster(lights, crs=newcrs)
lights[lights<0] <- 0
lights[lights>63] <- 63

# For weighted average
cust$nightlights = NaN
cust$nightlights = extract(lights, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

#### temperature, rainfall, and solar

## temperature
quads <- list.files("data/temperature", pattern = "\\.tif$")
temp <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/temperature/",quads[i]))
  r <- crop(r, extent(buffer))
  r <- mask(x=r, mask=buffer)
  temp <- stack(temp,r)
}
temp <- stackApply(temp, indices =  rep(1,nlayers(temp)), fun = "mean", na.rm = T)
temp <- projectRaster(temp, crs=newcrs)

# For weighted average
cust$temperature = NaN
cust$temperature = extract(temp, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

## rainfall
quads <- list.files("data/rainfall", pattern = "\\.tif$")
rain <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/rainfall/",quads[i]))
  r <- crop(r, extent(buffer))
  r <- mask(x=r, mask=buffer)
  rain <- stack(rain,r)
}
rain <- stackApply(rain, indices =  rep(1,nlayers(rain)), fun = "mean", na.rm = T)
rain <- projectRaster(rain, crs=newcrs)

# For weighted average
cust$rainfall = NaN
cust$rainfall = extract(rain, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

## solar
quads <- list.files("data/solar", pattern = "\\.tif$")
sun <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/solar/",quads[i]))
  r <- crop(r, extent(buffer))
  r <- mask(x=r, mask=buffer)
  sun <- stack(sun,r)
}
sun <- stackApply(sun, indices =  rep(1,nlayers(sun)), fun = "mean", na.rm = T)
sun <- projectRaster(sun, crs=newcrs)

# For weighted average
cust$solar = NaN
cust$solar = extract(sun, sp_cust_buffer, weights=T, fun=mean,na.rm=T) 

## solar wb
pvout <- raster("data/solar_LTAy_DailySum/PVOUT.tif")
pvout <- crop(pvout, extent(buffer))
pvout <- mask(x=pvout, mask=buffer)
pvout <- projectRaster(pvout, crs=newcrs)

# For weighted average
cust$pvout = NaN
cust$pvout = extract(pvout, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

#### vegetation, and elevation

##vegetation
quads <- list.files("data/vegetation", pattern = "\\.hdf$")
# Read in a quad 
sds <- get_subdatasets(paste0("data/vegetation/",quads[1]))
gdal_translate(sds[1], dst_dataset = paste0("data/vegetation/NDVI",toString(1),".tif"))
veg <- raster(paste0("data/vegetation/NDVI",toString(1),".tif")) 
# Repeat for other quads 
for (i in 2:length(quads)){
  sds <- get_subdatasets(paste0("data/vegetation/",quads[i]))
  gdal_translate(sds[1], dst_dataset = paste0("data/vegetation/NDVI",toString(i),".tif"))
  r <- raster(paste0("data/vegetation/NDVI",toString(i),".tif"))
  veg <- merge(veg, r)
}
veg = veg*0.00000001

# Clip to national boundaries
veg <- projectRaster(veg, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
veg <- crop(veg, extent(buffer))
veg <- mask(x=veg, mask=buffer)
veg <- projectRaster(veg, crs=newcrs)

# For weighted average
cust$vegetation = NaN
cust$vegetation = extract(veg, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

## elevation
quads <- list.files("data/elevation", pattern = "\\.tif$")
# Read in a quad 
elev <- raster(paste0("data/elevation/",quads[1])) 
# Repeat for other quads 
for (i in 2:length(quads)){
  r <- raster(paste0("data/elevation/",quads[i]))
  elev <- merge(elev, r)
}

# Clip to national boundaries
elev <- projectRaster(elev, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
elev <- crop(elev, extent(buffer))
elev <- mask(x=elev, mask=buffer)
elev <- projectRaster(elev, crs=newcrs)

# For weighted average
cust$elevation = NaN
cust$elevation = extract(elev, sp_cust_buffer, weights=T, fun=mean,na.rm=T)

#write csv file with all spatial data inputs
write.csv(cust, file="data/customer_profiles_inputs.csv")

#### Raster stack for Kenya for prediction

#### population, literacy, poverty, and nightlights

## population
pop <- raster("data/population/KEN_popmap15adj_v2b.tif")
pop <- aggregate(pop, fact=10, fun=sum)
pop <- crop(pop, extent(buffer))
pop <- mask(x=pop, mask=buffer)

## literacy
lit <- raster("data/literacy/KEN_literacy_F.tif") # need female population 09 to aggregate
lit <- crop(lit, extent(buffer))
lit <- mask(x=lit, mask=buffer)
lit <- projectRaster(lit, pop)

## poverty
pov <- raster("data/poverty/ken08povmpi.tif")  # need population 09 to aggregate
pov <- crop(pov, extent(buffer))
pov <- mask(x=pov, mask=buffer)
pov <- projectRaster(pov, pop)

## nightlights
lights <- raster("data/lights/F182013.v4c_web.stable_lights.avg_vis.tif")
lights <- crop(lights, extent(buffer))
lights <- mask(x=lights, mask=buffer)
lights <- projectRaster(lights, pop)

## temperature
quads <- list.files("data/temperature", pattern = "\\.tif$")
temp <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/temperature/",quads[i]))
  r <- crop(r, extent(buffer))
  r <- mask(x=r, mask=buffer)
  temp <- stack(temp,r)
}
temp <- stackApply(temp, indices =  rep(1,nlayers(temp)), fun = "mean", na.rm = T)
temp <- projectRaster(temp, pop)

## rainfall
quads <- list.files("data/rainfall", pattern = "\\.tif$")
rain <- stack()
for (i in 1:length(quads)){
  r <- stack(paste0("data/rainfall/",quads[i]))
  r <- crop(r, extent(buffer))
  r <- mask(x=r, mask=buffer)
  rain <- stack(rain,r)
}
rain <- stackApply(rain, indices =  rep(1,nlayers(rain)), fun = "mean", na.rm = T)
rain <- projectRaster(rain, pop)

## solar wb
pvout <- raster("data/solar_LTAy_DailySum/PVOUT.tif")
pvout <- crop(pvout, extent(buffer))
pvout <- mask(x=pvout, mask=buffer)
pvout <- projectRaster(pvout, pop)

##vegetation
veg <- raster(paste0("data/vegetation/NDVI",toString(1),".tif")) 
# Repeat for other quads 
for (i in 2:4){
  r <- raster(paste0("data/vegetation/NDVI",toString(i),".tif"))
  veg <- merge(veg, r)
}
veg = veg*0.00000001

# Clip to national boundaries
veg <- projectRaster(veg, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
veg <- aggregate(veg, fact=4, fun=mean)
veg <- crop(veg, extent(buffer))
veg <- mask(x=veg, mask=buffer)
veg <- projectRaster(veg, pop)

## elevation
quads <- list.files("data/elevation", pattern = "\\.tif$")
# Read in a quad 
elev <- raster(paste0("data/elevation/",quads[1])) 
# Repeat for other quads 
for (i in 2:length(quads)){
  r <- raster(paste0("data/elevation/",quads[i]))
  elev <- merge(elev, r)
}

# Clip to national boundaries
elev <- aggregate(elev, fact=10, fun=mean)
elev <- crop(elev, extent(buffer))
elev <- mask(x=elev, mask=buffer)
elev <- projectRaster(elev, pop)
elev[elev<0] <- 0

kenya_stack <- stack(pop,lit,pov,lights,temp,rain,pvout,veg,elev)
writeRaster(kenya_stack, filename="data/kenya_stack.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

kenya_df <- as.data.frame(kenya_stack, xy=TRUE, row.names=NULL, optional=FALSE)
write.csv(kenya_df, file="data/kenya_df.csv")
