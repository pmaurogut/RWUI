# TODO: Add comment
#
# Author: Paco
###############################################################################
library("terra")
library("raster")
library("sf")
setwd("C:/Users/UVa/OneDrive - UVa/CambiumWS/MFE_IFN")


create_template <- function(layer, buffer=0, resolution,origin){
	
	if(missing(layer)){
		stop("Layer needed to create extent")
	}
	
	if(inherits(layer,c("SpatVector","sf"))){
		
		if(inherits(layer,"sf")){
			layer <- vect(layer)
		}
		
	}else if(inherits(layer,c("SpatRaster","Raster"))) {
		if(inherits(layer,"Raster")){
			layer <- rast(layer)
		}
		resolution <- res(layer)
		origin <- origin(layer)
	}
	
	template_created <- try({
		template<-as.polygons(ext(layer))
		if(buffer>0){
			template <- buffer(template)
			template <- ext(template)
		}else{
			template <- ext(template)
		}
		template <- rast(extent= template,
				resolution = resolution, crs = crs(layer))	
		origin(template) <- origin
			})
	if(inherits(template_created,"try-error")){
		stop("Unable to create template")
	}else{
		return(template)
	}
	
} 
 
 

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
stewart_wui <- function(housing_density, is_forest, is_exposed, ...) {
	
	housing_density <- as(housing_density,"Raster")
	is_forest <- as(is_forest,"Raster")
	is_exposed <- as(is_exposed,"Raster")
    res <- raster::overlay(housing_density, is_forest, is_exposed,
            fun=function(x,y,z){reclass_stewart_wui(x,y,z)}
    )
	res <- rast(res)
	levels(res) <- c(
			"0", "a", "b", "Vegetated no Housing",
			"Dispersed rural", "Intermix", "Interface"
			)

    if (missing(...)) {
        return(res)
    } else {
        writeRaster(res, ...)
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
        template <- create_template(houses,buffer=resolution*blocks,
				resolution=resolution,origin=origin)
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
    resolution <- res(template)[1]
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
		forest_only<- try(terra::rasterize(vect(polys),forest_only, "id"))
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
		to_raster=TRUE,resoultion=150,origin=0, ...){
	
	if(missing(forest)){
		stop("Need forest layer to calculate exposure")
	}
	

	if(inherits(forest,c("sf","SpatVector"))){
		
		if(inherits(forest,c("SpatVector"))){
			forest <-st_as_sf(forest)
		}
		conversion <- try({
			forest <- forest[forest[[forest_field]]%in%forest_code,]
			})

		if (inherits(conversion, "try-error")) {
			stop("forest cannot be converted to sf")
		}
	}else 	if(inherits(forest,c("raster","SpatRaster"))){
		forest <- delineate_forest(forest,forest_code=forest_code)
	}

	
	filter_dissolve <- try({
		forest_crs <- crs(forest)
		to_rasterize <- st_sf(exposed=1,st_union(forest))
		to_rasterize$area <- st_area(to_rasterize)
		to_rasterize <- to_rasterize[as.numeric(to_rasterize$area)>min_patch_size,]
		to_rasterize <- st_buffer(to_rasterize,max_ember_dist)
		to_rasterize <- vect(to_rasterize)
			})

	if(inherits(filter_dissolve,"try-error")){
		stop("Could not dissolve forest polygons and filter polygons < min_patch_size")
	}
	
	if(to_raster){

		if (missing(template)) {
			if(inherits(forest,"SpatRaster")){
				template <- forest
			}else{
				template <- create_template(forest,buffer=0,
						resolution=resolution,origin=origin)
			}
		}
		
		if (!inherits(template, "SpatRaster")) {
			stop("template is not a SpatRaster and cannot be used as template")
		}

		created <- try({
			is_exposed <-terra::rasterize(to_rasterize,
					template,field="exposed",crs=crs(template))
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
			return(to_rasterize)
		}else{
			writeVector(to_rasterize,...)
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
		template <- create_template(layer,buffer=0,
				resolution=resolution,origin=origin)
	}
		
	if (!inherits(template, "SpatRaster")) {
		stop("template is not a SpatRaster and cannot be used as template")
	}

	if(!inherits(layer,"SpatVector")){
		layer <- try(vect(layer))
		if(inherits(layer,c("try-error"))){
			stop("Error converting layer to SpatVector")
		}
	}
	
	if(missing(forest_field)){
		stop("No field provided ")
	}
	
	if(filter){
		layer <- layer[layer[[forest_field]][,1]%in%forest_code,]
	}
	
	created <- try({
				is_veg <-terra::rasterize(layer,template,field=veg_field)
			})
	
	if(inherits(created,c("try-error"))){
		stop("Error rasterizing")
	}
	
	if(reclass){
		reclassed <- try({
					rclmat <- matrix(c(0, 50, 0,
									50, 1000, 1), ncol=3, byrow=TRUE)
					layer <- terra::classify(is_veg, rclmat, include.lowest=TRUE)
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


houses <- st_read("Data/CATASTRO/42/Merged_42.gpkg",
		layer = "Building")
houses <- houses[houses$currentUse %in% c("1_residential", "3_industrial", "4_3_publicServices"), ]
houses <- houses[houses$conditionOfConstruction == "functional", ]

forest <- st_read("Data/MFE/mfe_castillayleon/mfe_CastillayLeon.shp")
forest <- forest[forest$prov_nom=="Soria",]
forest <- st_transform(forest, crs(houses))
# forest <- st_write(forest, "",delete_layer=TRUE)


is_forest_mfe <- unique(forest[["lulucf"]])
is_forest_mfe <- is_forest_mfe[is_forest_mfe<=200]
is_forest_mfe <- c(is_forest_mfe, 400)

forest$is_forest <- 1

is_forest <- is_vegetated_vect(forest,veg_field="fccarb",
		filter= TRUE,forest_field= "lulucf", forest_code = is_forest_mfe,
		filename="Data/is_forest.tif",overwrite=TRUE)

is_exposed <- is_ember_exposed(forest,template=is_forest,
		forest_field= "lulucf",  forest_code = is_forest_mfe,
		max_ember_dist=2e3,min_patch_size=5e6,
		resoultion=150,origin=0,
		filename="Data/is_exposed.tif",overwrite=TRUE)

#is_exposed_poly <- is_ember_exposed(forest,template=is_forest,
#		forest_field= "lulucf",  forest_code = is_forest_mfe,
#		max_ember_dist=2e3,min_patch_size=5e6,to_raster=FALSE,
#		resoultion=150,origin=0,
#		filename="Data/is_exposed.gpkg",overwrite=TRUE)
#
#forest_poly <- is_exposed_poly <- is_ember_exposed(forest,template=is_forest,
#		forest_field= "lulucf",  forest_code = is_forest_mfe,
#		max_ember_dist=2e3,min_patch_size=0,to_raster=FALSE,
#		resoultion=150,origin=0,
#		filename="Data/forest.gpkg",overwrite=TRUE)


housing_density<-housing_dens(houses,is_forest,filename="Data/housing_density.tif",overwrite=TRUE)
classes <- stewart_wui(housing_density,is_forest,is_exposed,
		filename="Data/stewart.tif",overwrite=TRUE)
plot(classes)

