# TODO: Add comment
# 
# Author: Paco
###############################################################################

library("terra")

#' 
#' @param x 
#' @param res_ratio 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
normalize_cont <-function(x,res_ratio=30, ...){
	res<-100*x/(6*res_ratio^2)
	if(missing(...)){
		return(res)
	}else{
		writeRaster(res,...)
	}
}

#' 
#' @param x 
#' @param res_ratio 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
normalize_frict <-function(x,res_ratio=30, ...){
	
	if((res_ratio %% 2) == 0){
		factor <- (5/2)*(res_ratio^2)
	}else{
		factor <- (5/2)*(res_ratio^2+1)
	}
	
	res<-100*x/factor
	if(missing(...)){
		return(res)
	}else{
		writeRaster(res,...)
	}
}

#' 
#' @param x 
#' @param res_ratio 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
normalize_WUIX <-function(x,res_ratio=30, ...){
	
	if((res_ratio %% 2) == 0){
		factor <- (5/2)*(res_ratio^2)
	}else{
		factor <- (5/2)*(res_ratio^2+1)
	}
	factor <- factor*(6*res_ratio^2)
	res <- 100*x/factor
	if(missing(...)){
		return(res)
	}else{
		writeRaster(res,...)
	}
}

#' 
#' @param vegetation 
#' @param res_ratio 
#' @param normalize 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
continuity <- function(vegetation, res_ratio = 30,normalize=TRUE, ...) {
	
	if(!inherits(vegetation,"SpatRaster")){
		vegetation<-try(as(vegetation,"SpatRaster"))
		if(inherits(vegetation,"try-error")){
			stop("Vegetation is not a SpatRaster or cannot be coerced to one")
		}
	}
	
	resolution <- res(vegetation)[1]
	kernel <- resolution*matrix(c(0.5, 1, 0.5, 1, 0, 1, 0.5, 1, 0.5), ncol = 3)
	tmp <- terra::focal(vegetation, kernel, na.rm = TRUE)
	c_ij <- overlay(
			x = vegetation, y = tmp,
			fun = function(x, y) {
				ifelse(is.na(x) | is.na(y), 0, ifelse(x == 0, 0, y))
			}
	)
	result <- terra::aggregate(c_ij, res_ratio = res_ratio)
	result[result <= 0] <- NA
	
	if(normalize){
		result <-normalize_cont(result)
	}
	
	if(missing(...)){
		return(result)
	}else{
		writeRaster(result,...)
	}
	
	
}

#' 
#' @param vegetation 
#' @param buildings 
#' @param res_ratio 
#' @param normalize 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
friction <- function(vegetation, buildings, res_ratio = 30,normalize=TRUE, ...) {
	
	if(!inherits(vegetation,"SpatRaster")){
		vegetation<-try(as(vegetation,"SpatRaster"))
		if(inherits(vegetation,"try-error")){
			stop("Vegetation is not a SpatRaster or cannot be coerced to one")
		}
	}
	if(!inherits(buildings,"SpatRaster")){
		buildings<-try(as(buildings,"SpatRaster"))
		if(inherits(buildings,"try-error")){
			stop("Vegetation is not a SpatRaster or cannot be coerced to one")
		}
	}
	
	resolution <- res(vegetation)[1]
	kernel <- resolution*matrix(c(0.5, 1, 0.5, 1, 1, 1, 0.5, 1, 0.5), ncol = 3)
	tmp <- terra::focal(vegetation, kernel, na.rm = TRUE)
	tmp <- terra::overlay(tmp, buildings, fun = function(x, y) {
				x * y
			})
	c_ij <- terra::overlay(
			x = buildings, y = tmp,
			fun = function(x, y) {
				ifelse(is.na(x) | is.na(y), 0, ifelse(x == 0, 0, y))
			}
	)
	result <- aggregate(c_ij, res_ratio = res_ratio)
	result[result <= 0] <- NA
	
	if(normalize){
		result <-normalize_frict(result)
	}
	
	
	if(missing(...)){
		return(result)
	}else{
		writeRaster(result,...)
	}
}

#' 
#' @param vegetation 
#' @param buildings 
#' @param res_ratio 
#' @param normalize 
#' @param get_components 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author Paco
#' @export
wuix <- function(vegetation, buildings, res_ratio = 30,normalize=TRUE,get_components=FALSE, ...) {
	
	cont <- continuity(vegetation, res_ratio,normalize)
	fric <- friction(vegetation, buildings, res_ratio,normalize)
	result <- terra::overlay(cont, fric, fun = function(x, y) {
				x * y
			})
	if(normalize){
		result <-normalize_WUIX(result)
	}
	
	if(get_components){
		result <- rast(list(wuix = result, cont = cont, fric = fric))
	}
	
	if(missing(...)){
		return(result)
	}else{
		writeRaster(result,...)
	}
	
}

