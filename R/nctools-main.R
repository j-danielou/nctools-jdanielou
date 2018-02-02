
# Main functions for ncdf files -------------------------------------------


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
  on.exit(nc_close(nc))
  ncNew = nc_create(filename=output, vars=nc$var[[varid]])
  on.exit(try(nc_close(ncNew), silent = TRUE), add=TRUE)
  ncvar_put(ncNew, varid, ncvar_get(nc, varid, collapse_degen=FALSE))
  nc_close(ncNew)
  return(invisible(output))
}

#' Renaming variable and dimensions in a netCDF File
#'
#' @param filename The filename of the original ncdf file.
#' @param oldnames A string vector containing the names of the
#' variable or dimensions in the file that are to be renamed.
#' @param newnames A string vector containing the new names of
#' the variables or dimensions.
#' @param output Optional, the output file with the changes. By default,
#' it will overwrite the old file.
#' @param verbose If TRUE, run verbosely.
#'
#' @return
#' @export
#'
#' @examples
nc_rename = function(filename, oldnames, newnames, output, verbose=FALSE) {

  ncc = nc_open(filename)
  on.exit(try(nc_close(ncc), silent = TRUE))

  vars = names(nc$var)
  dims = names(nc$dim)

  gv = which(oldnames %in% vars)
  gd = which(oldnames %in% dims)
  nm = which(!(oldnames %in% c(vars, dims)))

  msgV = paste(sQuote(oldnames[gv]), sQuote(newnames[gv]), sep=" -> ", collapse="\n")
  msgD = paste(sQuote(oldnames[gd]), sQuote(newnames[gd]), sep=" -> ", collapse="\n")
  msgN = sprintf("Some of the 'oldnames' (%s) were not found in variables or dimensions.",
                 paste(sQuote(oldnames[nm]), collapse=", "))

  if(length(nm)>0) warning(msgN)

  old_varname = oldnames[gv]
  new_varname = newnames[gv]
  old_dimname = oldnames[gd]
  new_dimname = newnames[gd]

  nc_close(ncc)

  if(length(old_varname)>0) {

    if(isTRUE(verbose)) cat("Changing variable names:\n", msgV, "\n",sep="")

    filename = .nc_renameVar(filename=filename, oldname=old_varname,
                             newname=new_varname, output=output, verbose=verbose)

  }

  if(length(old_dimname)>0) {

    if(isTRUE(verbose)) cat("Changing dimension names:\n", msgD, "\n",sep="")

    output  = .nc_renameDim(filename=filename, oldname=old_dimname,
                            newname=new_dimname, output=output, verbose=verbose)

  }

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
#' @param filename The name of the ncdf file to subset.
#' @param varid The name of the variable to subset. If missing and only one variable in the file, that one is used.
#' @param output The name of the ncdf output file to create.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nc_subset = function(filename, varid, output, ...) {

  bounds = list(...)
  nc = nc_open(filename)

  if(missing(varid)) {
    msg = sprintf("ncdf file has more than one variable (%s), argument 'varid' must be specified.",
                  paste(sQuote(names(nc$var)), collapse=", "))
    if(length(nc$var)>1) stop(msg)
    varid = names(nc$var)[1]
  }

  dims = ncvar_dim(nc, varid, value=TRUE)

  .getIndex = function(x, bound, FUN, default=1) {
    FUN = match.fun(FUN)
    # longitude in a torus
    if(is.null(bound)) return(default)
    if(diff(bound)<0) stop("Upper bound is lower than lower bound.")
    out = which((x>=bound[1]) & (x<=bound[2]))
    if(length(out)==0)
      stop(sprintf("All index are out of bounds: (%s).",
                   paste(bound, collapse=", ")))
    return(FUN(out))
  }

  index = setNames(lapply(names(dims),
                          function(x) .getIndex(dims[[x]], bounds[[x]], FUN=identity, default=TRUE)),
                   names(dims))

  start = setNames(sapply(names(dims),
                          function(x) .getIndex(dims[[x]], bounds[[x]], FUN=min, default=1)),
                   names(dims))

  count = setNames(sapply(names(dims),
                          function(x) .getIndex(dims[[x]], bounds[[x]], FUN=length, default=-1)),
                   names(dims))

  x  = ncvar_get(nc, varid, collapse_degen=FALSE, start=start, count=count)

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
nc_unlim = function(filename, unlim, output=NULL) {
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
