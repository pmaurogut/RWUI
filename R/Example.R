# TODO: Add comment
# 
# Author: Paco
###############################################################################

library("lidR")
library("raster")
library("sf")
library("rgdal")
library("stars")
## Seccion 1
# buildings<-st_read("Data/Miguel/zonaestudio.shp") # nolint
# buildings<-st_buffer(buildings,1000)

buildings <- st_read("Data/Miguel/42150/Garray.gpkg")
st_write(buildings, "Data/Miguel/Edificios.shp", overwrite = TRUE)

bbox <- st_bbox(buildings)
lidar <- catalog("Data/Miguel")

has_veg <- function(x) {
	x <- x[!is.na(x)]
	if (any(x %in% c(3, 4, 5))) {
		return(list(has_veg=1))
	} else {
		return(list(has_veg=0))
	}
}

vegetation <- pixel_metrics(lidar, has_veg(Classification), 3)
vegetation[vegetation == 0] <- NA
writeRaster(vegetation, "Data/vegetacion.tif", overwrite = TRUE)

buildings$is_building <- 1
buildings <- as(buildings, "Spatial")
buildings_raster <- rasterize(buildings, vegetation, "is_building")
writeRaster(buildings_raster, "Data/edificios.tif", overwrite = TRUE)

vegetation <- raster("Data/vegetacion.tif")
buildings_raster <- raster("Data/edificios.tif")

continuity <- function(vegetation, fact = 10) {
	kernel <- matrix(c(0.5, 1, 0.5, 1, 0, 1, 0.5, 1, 0.5), ncol = 3)
	kernel <- kernel * raster::res(vegetation)[1]
	tmp <- raster::focal(vegetation, kernel, na.rm = TRUE)
	c_ij <- overlay(
			x = vegetation, y = tmp,
			fun = function(x, y) {
				ifelse(is.na(x) | is.na(y), 0, ifelse(x == 0, 0, y))
			}
	)
	result <- aggregate(c_ij, fact = fact)
	result[result <= 0] <- NA
	result
}

friction <- function(vegetation, buildings, fact = 30) {
	kernel <- matrix(c(0.5, 1, 0.5, 1, 1, 1, 0.5, 1, 0.5), ncol = 3)
	kernel <- kernel * raster::res(vegetation)[1]
	tmp <- raster::focal(vegetation, kernel, na.rm = TRUE)
	tmp <- overlay(tmp, buildings, fun = function(x, y) {
				x * y
			})
	c_ij <- overlay(
			x = buildings, y = tmp,
			fun = function(x, y) {
				ifelse(is.na(x) | is.na(y), 0, ifelse(x == 0, 0, y))
			}
	)
	result <- aggregate(c_ij, fact = fact)
	result[result <= 0] <- NA
	result
}

wuix <- function(vegetation, buildings, fact = 30) {
	cont <- continuity(vegetation, fact)
	fric <- friction(vegetation, buildings, fact)
	wuix <- overlay(cont, fric, fun = function(x, y) {
				x * y
			})
	stack(list(wuix = wuix, cont = cont, fric = fric))
}

Coarse<-data.frame(coarse=c(5,10,20,30))


cont <- continuity(vegetation, 5)
writeRaster(cont, "Data/Continuidad.tif", overwrite = TRUE,progress="window")
fric <- friction(vegetation, buildings_raster, 5)
writeRaster(fric, "Data/Friccion.tif", overwrite = TRUE,progress="window")
wuix_5 <- wuix(vegetation, buildings_raster, 5)
writeRaster(wuix[[1]], "Data/Wuix_5.tif", overwrite = TRUE,progress="window")
