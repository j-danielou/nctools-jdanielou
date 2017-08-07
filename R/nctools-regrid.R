
#' Regrid a variable from a ncdf file
#'
#' @param filename The filename of the original ncdf file.
#' @param varid The name of the variable to extract.
#' @param dim A vector including the names of the dimensions to interpolate,
#' e.g. c("longitude", "latitude")
#' @param new A list including information of the new coordinates, the names
#' should match the ones provides in \code{dim}.
#' @param mask An optional matrix providing the ocean/land mask. The maks will
#' be multiplied by each layer of the new regrided array, so dimensions must be
#' consistent. The mask contains 1s in the area for valid data and NAs for the
#' invalid ones (e.g. 1s in the ocean, NAs in the land).
#' @param output The filename to save the regrided variable.
#'
#' @return
#' @export
#'
#' @examples
nc_regrid = function(filename, varid, dim, new, mask=NULL, output, fill=FALSE,
                     radius=1) {

  nc =  nc_open(filename)
  x = ncvar_get(nc, varid, collapse_degen = FALSE)
  old = list(lon=ncvar_get(nc, dim[1]),
             lat=ncvar_get(nc, dim[2]))

  old$lon = kali::checkLongitude(old$lon)
  new$lon = kali::checkLongitude(new$lon)

  x = kali::regrid(object=x, old=old, new=new, mask=mask)
  if(isTRUE(fill))
    x = kali::fillMap(object=x, mask=new$mask, radius=radius, fill.value=NA)

  newVar = nc$var[[varid]]
  newVar$size = dim(x)
  newVar$chunksizes = NA

  .modifyDim = function(x, dim, index) {
    if(is.null(index[[x]])) return(dim[[x]])
    dim[[x]]$size = length(index[[x]])
    dim[[x]]$len = length(index[[x]])
    dim[[x]]$vals = index[[x]]
    return(dim[[x]])
  }

  newVar$dim = lapply(names(nc$dim), FUN=.modifyDim, dim=nc$dim, index=new)

  ncNew = nc_create(filename=output, vars=newVar)
  ncvar_put(ncNew, varid, x)
  nc_close(ncNew)

  nc_close(nc)

  return(invisible(output))
}

