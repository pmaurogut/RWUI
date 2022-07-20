# TODO: Add comment
#
# Author: Paco
###############################################################################
library("terra")
library("raster")
library("sf")
setwd("C:/Users/Paco/OneDrive - UVa/CambiumWS/MFE_IFN")
houses <- st_read("Data/CATASTRO/MADRID/MADRID_merged.gpkg", layer = "Building")
houses <- houses[houses$currentUse %in% c("1_residential", "3_industrial", "4_3_publicServices"), ]
houses <- houses[houses$conditionOfConstruction == "functional", ]

forest <- st_read("Data/MFE/mfe_Madrid/mfe_madrid.shp")
forest <- st_transform(forest, crs(houses))
# forest <- st_write(forest, "",delete_layer=TRUE)
is_forest_mfe <- c(
    "Pastizal-Matorral", "Herbazal-Pastizal",
    "Bosque", "Bosque de Plantación", "Bosquetes", "Prados",
    "Galerías arbustivas", "Arbustedos", "Msc desarb/suelo desnudo"
)

is_forest_mfe <- data.frame(lulucf = unique(forest$lulucf), is_forest = NA)

is_forest_mfe <- unique(forest$lulucf)
is_forest_mfe <- is_forest_mfe[is_forest_mfe < 200]
is_forest_mfe <- c(is_forest_mfe, 400)

field <- "lulucf"
resolution <- 150
origin <- 0
categories <- is_forest <- is_forest_mfe
layer <- forest


#' 
#' @param housing_dens 
#' @param is_vegetated 
#' @param is_exposed 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
stewart_wui <- function(housing_dens, is_vegetated, is_exposed, ...) {
    res <- terra::overlay(housing_dens, is_vegetated, is_exposed,
        fun = function(x, y, z) {
            reclass_stewart_wui(x, y, z)
        }
    )
    levels(resolution) <- c(
        "0", "a", "b", "Vegetated no Housing",
        "Dispersed rural", "Intermix", "Interface"
    )

    if (missing(...)) {
        return(resolution)
    } else {
        writeRaster(resolution, ...)
    }
}


#' 
#' @param housing_dens 
#' @param is_vegetated 
#' @param is_exposed 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
reclass_stewart_wui <- function(housing_dens, is_vegetated, is_exposed) {
    cat1 <- ifelse(housing_dens > 0 & housing_dens <= 49.42e-6 &
        is_vegetated %in% c(0, NA) &
        is_exposed %in% c(0, NA), 1, 0)
    cat2 <- ifelse(housing_dens > 49.42e-6 &
        is_vegetated %in% c(0, NA) &
        is_exposed %in% c(0, NA), 2, 0)
    cat3 <- ifelse(housing_dens %in% c(0, NA) &
        is_vegetated == 1, 3, 0)
    cat4 <- ifelse(housing_dens > 0 & housing_dens <= 6.18e-6 &
        is_vegetated == 1, 4, 0)
    cat5 <- ifelse(housing_dens > 6.18e-6 &
        is_vegetated == 1 &
        is_exposed %in% c(0, NA), 5, 0)
    cat6 <- ifelse(housing_dens > 6.18e-6 &
        is_vegetated %in% c(0, NA) & is_exposed == 1, 6, 0)

    return(cat1 + cat2 + cat3 + cat4 + cat5 + cat6)
}

#'
#' @param houses
#' @param aoi
#' @param resolution
#' @param blocks
#' @param origin
#' @param ...
#' @returnType
#' @return
#'
#' @author Paco
#' @export
housing_dens <- function(houses, template, resolution = 150, 
						blocks = 3, origin = 0, classes = FALSE, ...) {
    if (missing(template)) {
        template <- st_as_sf(st_as_sfc(st_bbox(houses)))
        template <- st_buffer(template, blocks * resolution)
        template <- vect(template)
        template <- rast(template, resolution = resolution, crs = crs(houses))
        origin(template) <- origin 
    }
    if (inherits(template, "sf")) {
        template <- st_as_sf(st_as_sfc(st_bbox(template)))
        template <- ext(vect(st_buffer(template, blocks * resolution)))
        template <- rast(template, origin = origin, resolution = resolution)
        if (inherits(template, "try-error")) {
            stop("Wrong template, it cannot be converted to SpatRaster")
        }
    }
    if (inherits(template, "Raster")) {
        template <- try(rast(template))
        if (inherits(template, "try-error")) {
            stop("Wrong template, it cannot be converted to SpatRaster")
        }
    }
    if (!inherits(template, "SpatRaster")) {
        stop("template is not a SpatRaster and cannot be used as template")
    }

    houses <- try({
        houses <- st_as_sf(houses)
    })

    centroids <- try({
        centroids <- st_centroid(houses)
        centroids <- vect(centroids)
        centroids
    })

    density <- terra::rasterize(crds(centroids),
        template,
        fun = length, background = 0
    )
    resolution <- resolution(template)[1]
    density <- terra::focal(density, w = matrix(1, nc = blocks, nr = blocks))
    density <- 1e6 * density / ((blocks * resolution)^2)

    if (missing(...)) {
        return(density)
    } else {
        writeRaster(density, ...)
    }
}


#' 
#' @param fcc 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
delineate_forest <- function(forest, forest_code = 1,to_raster=TRUE, ...){
	
	if(!inherits(forest,"SpatRast")){
		if(inherits(forest,"Raster")){
			forest<- try(rast(forest))
		}
		if(inherits(forest,"try-error")){
			stop("Failed creating template")
		}
	}

	forest[!forest==forest_code]<-NA

	polys <- terra::as.polygons(forest)
	polys <- sf::st_as_sf(polys)
	polys$area <- st_area(polys)
	polys <- polys[,c("area")]
	polys$id <- 1:dim(polys)[1]
	polys <- polys[,c("id","area")]
	
	if(to_raster){
		forest_only<- try(terra::rasterize(vect(polys),forest_only, "exposed"))
		if(missing(...)){
			return(forest_only)
		}else{
			writeRaster(forest_only,...)
		}
	}else{
		if(missing(...)){
			return(polys)
		}else{
			st_write(polys,...)
		}
	}
}


#' 
#' @param forest 
#' @param template 
#' @param forest_code 
#' @param forest_fiel 
#' @param buffer_dist 
#' @param max_ember_dist 
#' @param resoultion 
#' @param origin 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
is_ember_exposed <- function(forest,template,
		forest_field= "lulucf",  forest_code = 1,
		max_ember_dist=2e3,min_patch_size=5e6,
		resoultion=150,origin=0, ...){
	
	if(missing(forest)){
		stop("Need forest layer to calculate exposure")
	}

	if(inherits(forest,c("Spatial","sf","SpatVector"))){
		conversion <- try({
			forest <-st_as_sf(forest)
			forest <- forest[forest[[forest_field]]%in%forest_code,]
			forest$area <- st_area(forest)
			})

		if (!inherits(conversion, "sf")) {
			stop("forest cannot be converted to sf")
		}
	}

	if(inherits(forest,c("raster","SpatRaster"))){
		forest <- delineate_forest(forest,forest_code=forest_code)
	}

	forest <- forest[forest$area>min_patch_size,]
	forest <- st_buffer(forest,max_ember_dist)
	forest$exposed <- 1
	forest <- st_union(forest)

	if(to_raster){

		if(missing(template)){
			template <- try(ext(vect(forest)))
			template <- try(rast(extent= template, resolution = resolution, origin = origin, crs = crs(forest)))
		}

		if(inherits(template,c("try-error"))){
			stop("Cannot create template")
		}

		created <- try({
			is_exposed <-terra::rasterize(forest,template,field="exposed")
			origin((is_exposed))<-origin
		})

		if(inherits(created,c("try-error"))){
			stop("Error rasterizing")
		}

		if(missing(...)){
			return(is_exposed)
		}else{
			writeRaster(is_exposed,...)
		}

	}else{

		if(missing(...)){
			return(forest)
		}else{
			st_write(forest,...)
		}
	}
	
}

#' 
#' @param layer 
#' @param field 
#' @param resolution 
#' @param origin 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
is_vegetated_vect <- function(layer,template,veg_field="fccarb",
		filter= FALSE,forest_field= "lulucf", forest_code = 1,reclass=TRUE,
		resolution = 150, origin = 0, ...){
	
	if (missing(template)) {
		template <- st_as_sf(st_as_sfc(st_bbox(layer)))
		template <- st_buffer(template, blocks * resolution)
		template <- vect(template)
		template <- rast(template, resolution = resolution, crs = crs(layer))
		origin(template) <- origin 
	}
	
	if (inherits(template, c("Spatial","sf","SpatVector"))){
		
		if(inherits(template, c("Spatial","SpatVector"))){
			template <- st_as_sf(template)
		}
		template_created <- try({
				template <- st_as_sf(st_as_sfc(st_bbox(template)))
				template <- ext(vect(st_buffer(template, blocks * resolution)))
				template <- rast(template, origin = origin, resolution = resolution)
				})
		
		if (inherits(template_created, "try-error")) {
			stop("Wrong template, it cannot be converted to SpatRaster")
		}
		
	}
	if (inherits(template, "Raster")) {
		template <- try(rast(template))
		if (inherits(template, "try-error")) {
			stop("Wrong template, it cannot be converted to SpatRaster")
		}
	}
	if (!inherits(template, "SpatRaster")) {
		stop("template is not a SpatRaster and cannot be used as template")
	}
	if(missing(field)){
		stop("No field provided with an sf layer")
	}
	
	if(!inherits(layer,"SpatVector")){
		layer <- try(vect(layer))
		if(inherits(layer,c("try-error"))){
			stop("Error converting layer to SpatVector")
		}
	}
	if(inherits(layer,"SpatVector")){
		
		if(filter){
			layer <- layer[layer[[field]]%in%forest_code,]
		}
		
		created <- try({
					is_veg <-terra::rasterize(layer,template,field=veg_field)
				})
			
		if(inherits(created,c("try-error"))){
			stop("Error rasterizing")
		}
	}else{
		stop("layer could not be converted to SpatVector")
	}
	
	if(reclass){
		reclassed <- try({
					rclmat <- matrix(c(0, 50, 0,
									50, 1000, 1), ncol=3, byrow=TRUE)
					layer <- terra::classify(layer, rclmat, include.lowest=TRUE)
				})
		
		if(inherits(reclassed,c("try-error"))){
			stop("Error reclassifying")
		}
	}
	

	
	if(missing(...)){
		return(layer)
	}else{
		writeRaster(layer,...)
	}
}


forest$is_forest <- 1
is_forest <- is_vegetated_vect(forest,,veg_field="fccarb",
		filter= FALSE,forest_field= "lulucf", forest_code = is_forest_mfe,
		filename="Data/is_forest.tif",overwrite=TRUE)

is_exposed <- is_ember_exposed(forest,field="lulucf",forest_code = is_forest_mfe,
		filename="Data/is_exposed.tif",overwrite=TRUE)

is_vegetated<-get_forest_mfe(forest,field="lulucf",filename="Data/is_vegetated.tif",overwrite=TRUE)

housing_dens<-housing_dens(houses,is_vegetated,filename="Data/housing_density.tif",overwrite=TRUE)
classes <- stewart_wui(housing_dens,is_vegetated,is_exposed,
		filename="Data/stewart.tif",overwrite=TRUE)
plot(classes)


field="lulucf";forest_code = is_forest_mfe;
filename="Data/is_exposed.tif";overwrite=TRUE
max_ember_dist=2e3;min_patch_size=5e6


housing_dens_categories <- function(density,...){
	
	if(missing(density)){
		stop("Missing housing density")
	}else{
		if(!inherits(density,"SpatRaster")){
			density <- try(rast(density))
			if(inherits(density,"try-error")){
				stop("Housing density layer cannot be converted to SpatRaster")
			}
		}
	}
	
	rclmat <- matrix(c(0, 6.18e-6, 0,
					6.18e-6, 49.42e-6, 1,
					6.18e-6,Inf,2), ncol=3, byrow=TRUE)
	density <- classify(density, rclmat, include.lowest=TRUE)
	density <- levels(density,c("very-low","low","medium_high"))
	if(missing(...)){
		return(density)
	}else{
		writeRaster(density,...)
	}
	
	
}


#' 
#' @param layer 
#' @param field 
#' @param resolution 
#' @param origin 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
get_fccarb_mfe <- function(layer, resolution = 150, origin = 0, ...){
	
	if(!inherits(layer,c("Spatial","sf","SpatVect"))){
		stop("layer must be an spatial, sf or SpatVect Object")
	}
	if(inherits(layer,"Spatial")){
		layer <- as(layer,"sf")
		layer <- try(vect(layer))
		if(inherits(template,"try-error")){
			stop("Failed creating template")
		}
	}
	if(inherits(layer,"sf")){
		layer <- try(vect(layer))
		if(inherits(layer,"try-error")){
			stop("Failed converting layer to SpatVector")
		}
	}
	if(inherits(layer,"SpatVector")){
		template <- try(raster::raster(st_as_sf(layer),
						origin = origin,resolution=resolution))
		template <- try(rast(template))
		if(inherits(template,"try-error")){
			stop("Failed creating template")
		}
	}
	
	layer <- try(terra::rasterize(layer,template, "fccarb"))
	if(inherits(layer,"try-error")){
		stop("Failed rasterizing layer")
	}
	
	if(missing(...)){
		return(layer)
	}else{
		writeRaster(layer,...)
	}
	
}

fccarb_mfe <- get_fccarb_mfe(forest,filename="Data/Test_get_fccarb_mfe.tif",overwrite=TRUE)


