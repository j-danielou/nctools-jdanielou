

#' Change the Prime Meridian of a variable in a ncdf file.
#'
#' @param filename The filename of the original ncdf file.
#' @param output The name of the output file.
#' @param varid The name of the variable to extract.
#' @param primeMeridian position of the Prime Meridian. Use "center" for [-180,180]
#' range values and "left" for [0,360] range (Pacific centered).
#' @param verbose If TRUE, run verbosely.
#' @param overwrite overwrite output file if already exists?
#'
#' @return
#' @export
#'
#' @examples
nc_changePrimeMeridian = function(filename, output, varid=NA, primeMeridian="center",
                                  verbose=FALSE, overwrite=FALSE, compression=NA) {

  if(missing(output) & !isTRUE(overwrite))
    stop("output file is missing. Set 'overwrite' to TRUE to make changes in the original file.")

  if(missing(output)) output = filename

  if(file.exists(output) & !isTRUE(overwrite))
    stop("output file already exists. Set 'overwrite' to TRUE.")

  tmp = paste(output, ".temp", sep="")

  nc = nc_open(filename)
  on.exit(nc_close(nc))

  if(is.na(varid) & (length(nc$var)>1)) stop("More than one variable, you must specify 'varid'.")
  if(is.na(varid)) varid = names(nc$var)

  ivar = nc$var[[varid]]
  lon = ivar$dim[[1]]$vals
  ndim = length(ivar$dim)
  pm = findPrimeMeridian(lon)

  pmCheck = is.null(pm) | identical(pm, primeMeridian)

  if(isTRUE(pmCheck)) {
    warning("Longitud values are correct, nothing to do.")
    nc_close(nc)
    on.exit()
    if(!identical(filename, output))
      file.copy(filename, output, overwrite=TRUE)
    return(invisible(output))
  }

  newlon = checkLongitude(lon, primeMeridian = primeMeridian)
  ind = sort(newlon, index.return=TRUE)$ix
  ivar$dim[[1]]$vals = newlon[ind]
  if(!is.na(compression)) ivar$compression = compression
  ivar$chunksizes = NA

  newvar = c(list(x=ncvar_get(nc, varid, collapse_degen=FALSE),
                  drop=FALSE, i=ind), rep(TRUE, ndim-1))
  newvar = do.call('[', newvar)
  invisible(gc())
  ncNew = nc_create(filename=tmp, vars=ivar)
  on.exit(if(file.exists(tmp)) file.remove(tmp))
  ncvar_put(ncNew, varid, newvar)
  nc_close(ncNew)
  nc_close(nc)

  if(file.exists(output)) file.remove(output)
  file.rename(tmp, output)

  return(invisible(output))

}


#' Extract a mask from a ncdf file.
#'
#' @param filename The filename of the original ncdf file.
#' @param output The file to write the mask. If NULL, the default, a list with the mask is returned.
#'
#' @return A ncdf file with the mask if output is not NULL, and a list with the mask information.
#' @export
#'
#' @examples
nc_mask = function(filename, output=NULL) {

  nc = nc_open(filename)
  dims = ncvar_dim(nc)
  dimCheck = all(sapply(dims, identical, y=dims[[1]]))
  if(!dimCheck) {
    stop("Variables dimension don't match, cannot extract the grid.")
  }
  dims = dims[[1]]
  count = rep(1, length(dims))
  count[1:2] = -1
  idims = nc$dim[nc$var[[1]]$dimids[1:2]+1]
  ivar = ncvar_get(nc, varid=names(nc$var)[1], count=count)
  mask = !is.na(ivar)
  mask[!mask] = NA
  mask = 0 + mask
  storage.mode(mask) = "integer"

  out = ncvar_dim(nc, value=TRUE)[[1]][1:2]
  out$mask = mask

  if(!is.null(output)) {
    iVar = ncvar_def(name="mask", units="0/1", dim=idims, prec="integer",
                     missval=-9999, longname="grid mask")
    ncNew = nc_create(filename=output, vars=iVar)
    ncvar_put(ncNew, varid=iVar, vals=mask)
    nc_close(ncNew)
    return(invisible(out))
  }


  return(out)

}
