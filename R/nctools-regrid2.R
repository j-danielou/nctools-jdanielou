
#' Title
#'
#' @param filename
#' @param varid
#' @param dim
#' @param new
#' @param mask
#' @param output
#' @param fill
#' @param radius
#' @param log
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nc_regrid2 = function(filename, varid=NULL, dim, new, mask=NULL, output, fill=FALSE,
                     radius=1, log=TRUE, ...) {

  nc =  nc_open(filename)

  old = list(lon=ncvar_get(nc, dim[1]),
             lat=ncvar_get(nc, dim[2]))

  new$lon = checkLongitude(new$lon,
                           primeMeridian = findPrimeMeridian(old$lon),
                           sort=TRUE)

  .modifyDim = function(x, dim, index) {
    if(is.null(index[[x]])) return(dim[[x]])
    dim[[x]]$size = length(index[[x]])
    dim[[x]]$len = length(index[[x]])
    dim[[x]]$vals = index[[x]]
    return(dim[[x]])
  }

  newVars = list()
  vars = if(is.null(varid)) names(nc$var) else varid

  for(varid in vars) {
    newVar = nc$var[[varid]]
    newVar$size[1:2] = c(length(new$lon), length(new$lat))
    newVar$chunksizes = NA
    newVar$dim = lapply(names(nc$dim), FUN=.modifyDim, dim=nc$dim, index=new)
    newVars[[varid]] = newVar
  }

  ncNew = nc_create(filename=output, vars=newVars)

  for(varid in vars) {

    x = ncvar_get(nc, varid, collapse_degen = FALSE)
    if(isTRUE(log)) x = log(x + 1e-4)

    if(!isTRUE(fill)) {
      x = regrid(object=x, old=old, new=new, mask=mask)
    } else {
      x = regrid2(object=x, old=old, new=new, mask=mask, ...)
    }

    if(isTRUE(log)) x = exp(x) - 1e-4
    ncvar_put(ncNew, varid, x)

  }

  nc_close(ncNew)
  nc_close(nc)

  return(invisible(output))
}

