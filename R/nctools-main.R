
#' Create a new file with one variable from a ncdf file
#'
#' @param filename The filename of the original ncdf file.
#' @param varid The name of the variable to extract.
#' @param output The name of the output file.
#'
#' @return NULL, it creates the new file in the disk.
#' @export
#'
#' @examples
nc_extract = function(filename, varid, output) {
  nc = nc_open(filename)
  ncNew = nc_create(filename=output, vars=nc$var[[var]])
  ncvar_put(ncNew, var, ncvar_get(nc, var, collapse_degen=FALSE))
  nc_close(ncNew)
  return(invisible(output))
}

#' Concatenate records of the same variable from different ncdf files
#'
#' @param filenames A vector with the file names.
#' @param varid The name of the variable to concatenate.
#' @param output The name of the output file.
#'
#' @return NULL, it creates the new file in the disk.
#' @export
#'
#' @examples
nc_rcat = function(filenames, varid, output) {
  # add function validation
  # check for unlim
  for(i in seq_along(filenames)) {
    nc = nc_open(filenames[i])
    if(!any(ncdim_isUnlim(nc))) stop("Files don't have an unlimited dimension.")
    if(i==1) {
      isUnlim = ncdim_isUnlim(nc)[ncvar_dim(nc)[[varid]]]
      ncNew = nc_create(filename=output, vars=nc$var[[varid]])
      start = rep(1, length(isUnlim))
    }
    ncSize = nc$var[[varid]]$size
    count = ncSize*isUnlim -1*!isUnlim
    ncvar_put(ncNew, varid, ncvar_get(nc, varid, collapse_degen=FALSE),
              start=start, count=count)
    start = start + ncSize*isUnlim
    nc_close(nc)
  }
  nc_close(ncNew)
  return(invisible(output))
}


#' Subset a ncdf variable
#'
#' @param filename
#' @param varid
#' @param output
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nc_subset = function(filename, varid, output, ...) {

  bounds = list(...)
  nc = nc_open(filename)
  x  = ncvar_get(nc, varid, collapse_degen=FALSE)
  dims = ncvar_dim(nc, varid, value=TRUE)
  dimnames(x) = dims

  .isInside = function(x, bound) {
    # longitude in a torus
    if(is.null(bound)) return(TRUE)
    if(diff(bound)<0) stop("Upper bound is lower than lower bound.")
    out = which((x>=bound[1]) & (x<=bound[2]))
    return(out)
  }

  index = setNames(lapply(names(dims),
                          function(x) .isInside(dims[[x]], bounds[[x]])),
                   names(dims))

  x = do.call("[", c(list(x), index))

  newVar = nc$var[[varid]]
  newVar$size = dim(x)
  newVar$chunksizes = NA

  .modifyDim = function(x, dim, index) {
    if(isTRUE(index[[x]])) return(dim[[x]])
    dim[[x]]$size = length(index[[x]])
    dim[[x]]$len = length(index[[x]])
    dim[[x]]$vals = dim[[x]]$vals[index[[x]]]
    return(dim[[x]])
  }

  newVar$dim = lapply(names(nc$dim), FUN=.modifyDim, dim=nc$dim, index=index)

  ncNew = nc_create(filename=output, vars=newVar)
  ncvar_put(ncNew, varid, x)
  nc_close(ncNew)

  nc_close(nc)

  return(invisible(output))
}



#' Make a dimension unlimited
#'
#' @param filename
#' @param unlim Name of the variable to set as unlimited
#' @param output Name of the output file. If NULL,
#' replace the original value
#'
#' @return
#' @export
#'
#' @examples
ncdim_unlim = function(filename, unlim, output=NULL) {
  # open ncdf connection
  if(is.null(output)) output = filename
  outputTemp = paste(output, ".temp", sep="")
  nc = nc_open(filename)

  .makeUnlim = function(x, unlim) {
    names(x$dim) = sapply(x$dim, "[[", "name")
    if(is.null(x$dim[[unlim]])) return(x)
    x$dim[[unlim]]$unlim = TRUE
    x$unlim = TRUE
    return(x)
  }

  # new variables with unlimited dimension
  newVars = lapply(nc$var, FUN=.makeUnlim, unlim=unlim)

  ncNew = nc_create(filename=outputTemp, vars=newVars)

  for(iVar in names(newVars))
    ncvar_put(ncNew, iVar, ncvar_get(nc, iVar, collapse_degen=FALSE))

  nc_close(ncNew)
  nc_close(nc)

  renameFlag = file.rename(outputTemp, output)

  return(invisible(newVars))

}



# Extra tools -------------------------------------------------------------

#' Data output in ncdf format
#'
#' @param x
#' @param filename
#' @param varid
#' @param dim
#' @param longname
#' @param units
#' @param prec
#' @param missval
#' @param compression
#' @param chunksizes
#' @param verbose
#' @param dim.units
#' @param dim.longname
#' @param unlim
#'
#' @return
#' @export
#'
#' @examples
write.ncdf = function(x, filename, varid, dim, longname, units, prec="float",
                      missval=-9999, compression=9, chunksizes=NA, verbose=FALSE,
                      dim.units, dim.longname, unlim=FALSE) {

  if(missing(dim)) dim = lapply(base::dim(x), seq_len)

  if(length(dim)!=length(dim(x)))
    stop("dim argument does not match data dimension.")

  if(is.null(names(dim)))
    dim = setNames(dim, paste("dim", seq_along(dim(x)), sep=""))

  if(missing(longname)) longname = ""
  if(missing(units))    units    = ""

  if(missing(dim.units)) dim.units = rep("", length(dim))
  if(length(dim.units)!=length(dim))
    stop("dim units provided are not equal to dimension size.")

  if(missing(dim.longname)) dim.longname = rep("", length(dim))
  if(length(dim.longname)!=length(dim))
    stop("dim longnames provided are not equal to dimension size.")

  dims = list()
  for(i in seq_along(dim))
    dims[[names(dim)[i]]] =
    ncdim_def(name=names(dim)[i], units=dim.units[i], vals=dim[[names(dim)[i]]],
              unlim=names(dim)[i]==unlim, longname=dim.longname[i])

  iVar = ncvar_def(name=varid, units=units, dim=dims, prec=prec ,missval=missval, longname=longname,
                   compression=compression, chunksizes=chunksizes, verbose=verbose)

  ncNew = nc_create(filename=filename, vars=iVar, verbose=verbose)

  ncvar_put(ncNew, varid=iVar, vals=x, verbose=verbose)

  nc_close(ncNew)

  return(invisible(filename))

}

