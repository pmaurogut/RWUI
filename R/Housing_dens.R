# TODO: Add comment
#
# Author: Paco
###############################################################################
library("terra")
library("raster")
library("sf")
setwd("C:/Users/Paco/OneDrive - UVa/CambiumWS/MFE_IFN")


create_template <- function(layer, buffer=0,origin=0, resolution){
	
	if(missing(layer)){
		stop("Layer needed to create extent")
	}
	
	if(inherits(layer,c("SpatVector","sf"))){
		
		if(missing(resolution)){
			stop("resolution must be passed if layer is a vector")
		}
		if(inherits(layer,"sf")){
			template <- st_as_sf(st_as_sfc(st_bbox(layer)))
		}else{
			template<-st_as_sf(as.polygons(ext(layer)))
		}
		
	}else if(inherits(layer,c("SpatRaster","Raster"))) {
		
		if(inherits(layer,"Raster")){
			layer <- rast(layer)
		}
		resolution <- res(layer)
		origin <- origin(layer)
		template<-st_as_sf(as.polygons(ext(layer)))
	}
	print(class(template))
	tryCatch({
		if(buffer>0){
			template <- st_buffer(template, dist=buffer)
			template <- vect(template)
			template <- ext(template)
			
		}else{
			template <- ext(vect(template))
		}
		print(class(template))
		template <- rast(extent= template,
				resolution = resolution, crs = crs(layer,proj=TRUE)@projargs)	
		origin(template) <- origin
		return(template)
			},error=function(e){
			print("Unable to create template")	
			return(e)
		}
	)
	
} 


reclass_stewart_wui <- function(housing_dens, is_vegetated, is_exposed) {
	
	
#	branch for cover < 50% values in other branch are tagged as NA
	result <- ifelse(is_vegetated==1,
			ifelse(housing_dens==0,3,ifelse(housing_dens<6.18,4,5)),NA)
	
#	branch for cover >=50% tag exposed areas as 10 and non exposed (<2km) as 11
	result <- ifelse(is.na(result),ifelse(is_exposed==1,10,11),result)
	
#	branch for cover >=50% tag exposed areas as 10 and non exposed (<2km) as 11
	result <-ifelse(result==10,ifelse(housing_dens>=6.18,6,7),result)

	result <-ifelse(result==11,ifelse(housing_dens>=49.42,2,1),result)
	
	return(result)
}

stewart_wui <- function(housing_density, is_forest, is_exposed, ...) {
	
	housing_density <- as(housing_density,"Raster")
	is_forest <- as(is_forest,"Raster")
	is_exposed <- as(is_exposed,"Raster")
    res <- raster::overlay(housing_density, is_forest, is_exposed,
            fun=function(x,y,z){reclass_stewart_wui(x,y,z)}
    )
	res <- rast(res)
	levels(res) <- c(
			"Very low & low housing dens", "Medium & high housing dens", 
			"Vegetated no Housing","Dispersed rural", "Intermix",
			"Interface", "Other"
			)

    if (missing(...)) {
        return(res)
    } else {
        writeRaster(res, ...)
    }
}

housing_dens <- function(houses, template, origin = 0,resolution = 150, 
						distance = 225, kernel="circular",  classes = FALSE, ...) {
    
	if (missing(template)|inherits(template,c("SpatVector","sf","Raster"))) {
        template <- create_template(houses,buffer=2*distance,
			origin=origin,resolution=resolution)
    }
    if (!inherits(template, "SpatRaster")) {
        stop("template is not a SpatRaster and it cannot be converted to one")
    }

    density <- tryCatch({
		houses <- st_as_sf(houses)
        centroids <- st_centroid(houses)
        centroids <- vect(centroids)
       	terra::rasterize(crds(centroids),
				template,fun = length, background = 0)
    }, error= function(e){
		print("error rasterizing houses")
		return(e)
	})

	tryCatch({
		resolution <- res(template)[1]
		if(kernel=="circular"){
			base_vect <- -ceiling(distance/resolution):ceiling(distance/resolution)
			w <- outer(base_vect,base_vect,
					function(x,y,res){res*sqrt(x^2+y^2)},res=resolution)
			w[]<-ifelse(w[]>distance,0,1)
#		total area kernel in km2
			den <- sum(w[])*resolution*resolution/1e6
			w[]<-w[]/den
		}else{
			blocks <- 2*ceiling(distance/resolution)+1
			w <- matrix(1, nc = blocks, nr = blocks)
			w <- 1e6*w/(length(blocks) * resolution)^2
		}
		
		return(terra::focal(density,w=w,...))
		
		}, error= function(e){
			print("Error generating housing density from raster of houses counts")
			return(e)
		})

}


delineate_forest <- function(forest, forest_code = 1,to_raster=TRUE, ...){
	
	if(!inherits(forest,"SpatRast")){
		if(inherits(forest,"Raster")){
			forest<- try(rast(forest))
		}
		if(inherits(forest,"try-error")){
			stop("Forest is not a SpatRaster or can be converted to one")
		}
	}

	forest[!forest==forest_code]<-NA

	polys <- terra::as.polygons(forest)
	polys <- sf::st_as_sf(polys)
	polys$id <- 1:dim(polys)[1]
	polys <- polys[,"id"]

	if(to_raster){
		try(terra::rasterize(vect(polys),forest_only, "id",...))
	}else{
		if(missing(...)){
			return(polys)
		}else{
			st_write(polys,...)
		}
	}
}

is_ember_exposed <- function(forest,template,
		forest_field= "lulucf",  forest_code = 1,
		max_ember_dist=2e3,min_patch_size=5e6,
		raster_distance=TRUE, to_raster=TRUE,resoultion=150,origin=0, ...){
	
	if(missing(forest)){
		stop("Need forest layer to calculate exposure")
	}
	

	if(inherits(forest,c("sf","SpatVector"))){
		
		if(inherits(forest,c("SpatVector"))){
			forest <-st_as_sf(forest)
		}
		filtered <- try({
			forest <- forest[forest[[forest_field]]%in%forest_code,]
			})

		if (inherits(filtered, "try-error")) {
			stop("forest cannot filter forest areas")
		}
			
	}else 	if(inherits(forest,c("raster","SpatRaster"))){
		
		
		delineated1 <- try({
					forest <- delineate_forest(forest,forest_code,to_raster=FALSE)
				})
		
		if(inherits(delineated1,c("try-error"))){
			stop("Error delineating forest areas")
		}
		
	}

	
	filter_dissolve <- try({
				forest_crs <- crs(forest)
				filtered <- st_sf(exposed=1,st_union(forest))
				filtered$area <- st_area(filtered)
				filtered <- filtered[as.numeric(filtered$area)>min_patch_size,]
			})
	
	if(inherits(filter_dissolve,"try-error")){
		stop("Could not dissolve forest polygons and filter polygons < min_patch_size")
	}
	
	if(raster_distance){
		
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
		
		filtered <- vect(filtered)
		created <- try({
					filtered$large_forest<-1
					forest <-terra::rasterize(filtered,
							template,field="large_forest",crs=crs(template))
					is_exposed <- terra::distance(forest)
					is_exposed[] <- ifelse(is_exposed[]>max_ember_dist,NA,1)
				})
		
		if(to_raster==FALSE){
			
			tryCatch({
					delineate_forest(is_exposed,1,to_raster=FALSE,...)
				},error=function(e){
					stop("Error delineating exposed areas")
					e
				})
		}else{
			if(missing(...)){
				return(is_exposed)
			}else{
				writeRaster(is_exposed,...)
			}
		}

		
	}else{
		
		buffer <- try({
				is_exposed <- st_buffer(filtered,max_ember_dist)
				is_exposed$is_exposed <- TRUE
				is_exposed <- vect(is_exposed)
				})
		
		if(inherits(buffer,c("try-error"))){
			stop("Error computing buffer polygons")
		}
			
		if(to_raster==FALSE){
			
			if(missing(...)){
				return(is_exposed)
			}else{
				writeVector(is_exposed,...)
			}
			
		}else{
			
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
			
			terra::rasterize(is_exposed,
					template,field="is_exposed",crs=crs(template),...)
			
		}
		
	}
	
}


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
		tryCatch({
			rclmat <- matrix(c(0, 50, 0,
					50, 1000, 1), ncol=3, byrow=TRUE)
			layer <- terra::classify(is_veg, rclmat, include.lowest=TRUE)
				},
			error=function(e){
				print("Error reclassifying")
				e
				})
	}

}


vegetated_cover_vect <- function(layer,template,filter=TRUE,
		forest_field= "lulucf", forest_code = 1,reclass=TRUE,
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
		layer$is_forest <- TRUE
	}
	
	created <- try({
				layer <- st_as_sf(layer)
				template <- raster(template)
				forest_cover <-raster::rasterize(layer,template,
						getCover=TRUE,
						filename="test.tif",progress="window")
				forest_cover <- rast(forest_cover)
			})
	
	if(inherits(created,c("try-error"))){
		stop("Error rasterizing")
	}
	
	if(reclass){
		reclassed <- try({
				rclmat <- matrix(c(0, 0.50, 0,0.50, 2, 1), 
						ncol=3, byrow=TRUE)
				layer <- terra::classify(forest_cover, rclmat, include.lowest=TRUE)
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

vegetated_cover_rast <- function(layer,template,filter=TRUE,
		forest_code = 1,reclass=TRUE,
		resolution = 150, origin = 0, ...){
	
	
	if (missing(layer)) {
		stop("layer is missing")
	}
	
	if (!inherits(template, "SpatRaster")) {
		stop("template is not a SpatRaster and cannot be used as template")
	}
	
	if(inherits(layer,c("raster","SpatRaster"))){
			
		delineated1 <- try({
					forest <- delineate_forest(layer,forest_code,to_raster=FALSE)
					forest$is_forest <- 1
				})
		
		if(inherits(delineated1,c("try-error"))){
			stop("Error delineating forest areas")
		}
		
	}else{
		stop("layer must be a raster or SpatRater object")
	}
	
	if(missing(template)){
		vegetated_cover_vect(forest,filter=filter,forest_field="is_forest",
				forest_code=1,reclass=reclass,resolution=resolution,origin=origin,...)
	}else{
		vegetated_cover_vect(forest,template=template,filter=filter,forest_field="is_forest",
				forest_code=1,reclass=reclass,resolution=resolution,origin=origin,...)
	}
	
}



houses <- st_read("Data/CATASTRO/CATASTRO_CyL.gpkg",
		layer = "Building")
houses <- houses[houses$currentUse %in% c("1_residential", "3_industrial", "4_3_publicServices"), ]
houses <- houses[houses$conditionOfConstruction == "functional", ]

forest <- st_read("Data/MFE/mfe_castillayleon/mfe_CastillayLeon.shp")
forest <- st_transform(forest, crs(houses))
# forest <- st_write(forest, "",delete_layer=TRUE)

is_forest_mfe <- unique(forest[["lulucf"]])
is_forest_mfe <- is_forest_mfe[is_forest_mfe<=200]
is_forest_mfe <- c(is_forest_mfe, 400)

forest$is_forest <- 1

forest2 <- forest[forest$prov_nom=="Soria",]

cover <- vegetated_cover_vect(forest2,forest_field= "lulucf", 
		forest_code = is_forest_mfe,reclass=FALSE,
		filename="Output/forest_cover.tif",overwrite=TRUE)

is_exposed <- is_ember_exposed(forest2,template=cover,
		forest_field= "lulucf",  forest_code = is_forest_mfe,
		max_ember_dist=2e3,min_patch_size=5e6,
		resoultion=150,origin=0,
		filename="Output/is_exposed.tif",overwrite=TRUE)

housing_density<-housing_dens(houses,cover,
		filename="Output/housing_density.tif",overwrite=TRUE)

classes <- stewart_wui(housing_density,cover,is_exposed,
		filename="Data/stewart.tif",overwrite=TRUE,progress=1)
plot(classes)
