# TO_DO: check on dimension names more properly
# TO_DO: integrate to the flow of other functions
# TO_DO: check on output argument

#' Regrid a variable from a ncdf file
#'
#' @param filename The filename of the original ncdf file.
#' @param varid What variable to read the data from. Can be a string with the
#' name of the variable or an object of class ncvar4. If set to NA,
#' the function will determine if there is only one variable in the file and,
#' if so, read from that, but if there are multiple variables in the file, an error is generated.
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
nc_regrid = function(filename, varid=NA, dim=1:2, new, mask=NULL, output, extrap=FALSE,
                     log=TRUE, keepRange=FALSE, ...) {

  nc =  nc_open(filename)
  on.exit(nc_close(nc))

  varid = .checkVarid(varid=varid, nc=nc)

  x = ncvar_get(nc, varid, collapse_degen = FALSE)
  if(isTRUE(log)) x = log(x + 1e-4)

  if(length(dim)!=2) stop("Regridding only possible for two dimensions.")

  if(is.numeric(dim)) {
    dim = ncvar_dim(nc, varid)[dim]
  }

  if(is.character(dim)) {
    checkDim = dim %in% ncvar_dim(nc, varid)
    if(!all(checkDim)) stop("Dimension names (dim) not found for 'varid'.")
    checkDim = dim %in% names(new)
    if(!all(checkDim)) stop("Dimension names (dim) not found in 'new'.")
  }

  old = list(lon=ncvar_get(nc, dim[1]),
             lat=ncvar_get(nc, dim[2]))

  new$lon = new[[dim[1]]]
  new$lat = new[[dim[2]]]

  new$lon = checkLongitude(new$lon,
                           primeMeridian = findPrimeMeridian(old$lon),
                           sort=FALSE)

  sameGrid = setequal(old$lon, new$lon) & setequal(old$lat, new$lat)

  nc_close(nc)
  on.exit()

  if(isTRUE(sameGrid)) {

    message("Old and new grids are the same, nothing to do. Copying file to output destination.")
    file.copy(from=filename, to=output)
    return(invisible(output))

  }

  if(!isTRUE(extrap)) {
    x = regrid(object=x, old=old, new=new, mask=mask, keepRange=keepRange)
  } else {
    x = regrid2(object=x, old=old, new=new, mask=mask, linear=FALSE,
                extrap=extrap, keepRange=keepRange, ...)
  }

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

  if(!dir.exists(dirname(output)))
    dir.create(dirname(output), recursive = TRUE)

  ncNew = nc_create(filename=output, vars=newVar)

  if(isTRUE(log)) x = exp(x) - 1e-4
  ncvar_put(ncNew, varid, x)
  nc_close(ncNew)

  return(invisible(output))

}


# regrid ------------------------------------------------------------------

regrid = function(object, old, new, mask, keepRange, ...) {
  UseMethod("regrid")
}

#' @export
regrid.matrix = function(object, old, new, mask=NULL, keepRange=FALSE, ...) {

  if(is.null(mask)&!is.null(new$mask)) mask=new$mask
  stopifnot(exists("lat", where=old), exists("lon", where=old))
  stopifnot(length(old$lat)==ncol(object), length(old$lon)==nrow(object))

  stopifnot(exists("lat", where=new), exists("lon", where=new))
  stopifnot(is.numeric(new$lat), !is.matrix(new$lat),
            is.numeric(new$lon), !is.matrix(new$lon))

  old$x = old$lon
  old$y = old$lat
  old$z = object

  new$x = new$lon
  new$y = new$lat

  newp = fields::interp.surface.grid(obj=old, grid.list=new, ...)$z
  newmap = if(!is.null(mask)) newp*mask else newp

  if(keepRange) {
    rr = range(object, na.rm=TRUE)
    newmap[newmap<rr[1]] = rr[1]
    newmap[newmap>rr[2]] = rr[2]
  }

  return(newmap)
}


#' @export
regrid.array = function(object, old, new, mask=NULL, keepRange=FALSE, ...) {
  # new grid
  if(exists("LAT", where=new) & exists("LON", where=new)) {
    stopifnot(is.matrix(new$LAT), is.matrix(new$LON), dim(new$LAT)==dim(new$LON))
    nLAT = new$LAT
    nLON = new$LON
    nlat = as.numeric(nLAT)
    nlon = as.numeric(nLON)
  } else {
    if(exists("lat", where=new) & exists("lon", where=new)) {
      stopifnot(is.numeric(new$lat), !is.matrix(new$lat), is.numeric(new$lon), !is.matrix(new$lon))
      nLAT = matrix(new$lat, ncol=length(new$lat), nrow=length(new$lon), byrow=TRUE)
      nLON = matrix(new$lon, ncol=length(new$lat), nrow=length(new$lon))
      nlat = as.numeric(nLAT)
      nlon = as.numeric(nLON)
    } else {
      stop("'new' must contain latitude and longitude information.")
    }
  }

  ndim = seq_along(dim(object))[-c(1,2)]
  newmap = apply(object, ndim, regrid, old=old, new=new, mask=mask, keepRange=keepRange, ...)
  dim(newmap) = c(dim(nLAT), dim(object)[-c(1,2)])
  return(newmap)
}


# regrid2 -----------------------------------------------------------------


regrid2 = function(object, old, new, mask, linear, extrap, keepRange, ...) {
  UseMethod("regrid2")
}


#' @export
regrid2.matrix = function(object, old, new, mask=NULL, linear=TRUE, extrap=FALSE, keepRange=FALSE, ...) {

  if(is.null(mask)&!is.null(new$mask)) mask=new$mask
  stopifnot(exists("lat", where=old), exists("lon", where=old))
  stopifnot(length(old$lat)==ncol(object), length(old$lon)==nrow(object))

  stopifnot(exists("lat", where=new), exists("lon", where=new))
  stopifnot(is.numeric(new$lat), !is.matrix(new$lat),
            is.numeric(new$lon), !is.matrix(new$lon))

  new$x = new$lon
  new$y = new$lat

  old$x = rep(old$lon, ncol(object))
  old$y = rep(old$lat, each=nrow(object))
  old$z = as.numeric(object)
  old = old[c("x","y","z")]
  old = as.data.frame(old)

  naobj = is.na(object)
  xrow = sum(apply(naobj, 1, all, na.rm=TRUE))
  xcol = sum(apply(naobj, 2, all, na.rm=TRUE))

  ox = 1

  if(xrow>0 | xcol>0) {

    omask = 0 + !is.na(object)
    oo    = data.frame(x=old$x, y=old$y, z=as.numeric(omask))
    oo    = oo[complete.cases(oo), ]
    nmask = akima::interp(x=oo$x, y=oo$y, z=oo$z, xo=new$x, yo=new$y,
                          linear=linear, extrap=extrap, ...)$z
    thr = 0.2
    nmask =  (nmask < thr) & !is.na(mask)
    nmask[is.na(mask)] = NA
    irow = apply(nmask, 1, all, na.rm=TRUE)
    icol = apply(nmask, 2, all, na.rm=TRUE)
    ox = !outer(irow, icol, FUN = "|")
    ox[!ox] = NA
  }

  old = old[complete.cases(old), ]
  newp = akima::interp(x=old$x, y=old$y, z=old$z, xo=new$x, yo=new$y,
                linear=linear, extrap=extrap, ...)$z
  newmap = if(!is.null(mask)) newp*mask else newp
  newmap = newmap*ox

  if(keepRange) {
    rr = range(object, na.rm=TRUE)
    newmap[newmap<rr[1]] = rr[1]
    newmap[newmap>rr[2]] = rr[2]
  }

  return(newmap)
}


#' @export
regrid2.array = function(object, old, new, mask=NULL, linear=TRUE, extrap=FALSE, keepRange=FALSE, ...) {
  # new grid
  if(exists("LAT", where=new) & exists("LON", where=new)) {
    stopifnot(is.matrix(new$LAT), is.matrix(new$LON), dim(new$LAT)==dim(new$LON))
    nLAT = new$LAT
    nLON = new$LON
    nlat = as.numeric(nLAT)
    nlon = as.numeric(nLON)
  } else {
    if(exists("lat", where=new) & exists("lon", where=new)) {
      stopifnot(is.numeric(new$lat), !is.matrix(new$lat),
                is.numeric(new$lon), !is.matrix(new$lon))
      nLAT = matrix(new$lat, ncol=length(new$lat), nrow=length(new$lon), byrow=TRUE)
      nLON = matrix(new$lon, ncol=length(new$lat), nrow=length(new$lon))
      nlat = as.numeric(nLAT)
      nlon = as.numeric(nLON)
    } else {
      stop("'new' must contain latitude and longitude information.")
    }
  }

  ndim = seq_along(dim(object))[-c(1,2)]
  newmap = apply(object, ndim, regrid2, old=old, new=new, mask=mask,
                 linear=linear, extrap=extrap, keepRange=keepRange, ...)
  dim(newmap) = c(dim(nLAT), dim(object)[-c(1,2)])
  return(newmap)
}

