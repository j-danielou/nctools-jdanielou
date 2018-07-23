
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
#' @param overwrite overwrite output file if already exists?
#'
#' @return
#' @export
#'
#' @examples
nc_rename = function(filename, oldnames, newnames, output, verbose=FALSE, overwrite=FALSE) {

  if(missing(output) & !isTRUE(overwrite))
    stop("output file is missing. Set 'overwrite' to TRUE to make changes in the original file.")

  if(missing(output)) output = filename

  if(file.exists(output) & !isTRUE(overwrite))
    stop("output file already exists. Set 'overwrite' to TRUE.")

  tmp = paste(output, ".temp", sep="")
  on.exit(if(file.exists(tmp)) file.remove(tmp))

  ncc = nc_open(filename)

  vars = names(ncc$var)
  dims = names(ncc$dim)

  nc_close(ncc)

  gv = which(oldnames %in% vars)
  gd = which(oldnames %in% dims)

  if(length(gv)==0 & length(gd)==0) {
    message("Nothing to change. Exiting...")
    return(invisible())
  }

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


  if(length(old_dimname)>0) {

    if(isTRUE(verbose)) cat("Changing dimension names:\n", msgD, "\n",sep="")

    filename  = .nc_renameDim(filename=filename, oldname=old_dimname,
                            newname=new_dimname, output=tmp, verbose=FALSE)

  }

  if(length(old_varname)>0) {

    if(isTRUE(verbose)) cat("Changing variable names:\n", msgV, "\n",sep="")

    filename = .nc_renameVar(filename=filename, oldname=old_varname,
                             newname=new_varname, output=tmp, verbose=FALSE)

  }

  if(file.exists(output)) file.remove(output)
  file.rename(tmp, output)

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
      isUnlim  = ncdim_isUnlim(nc)[ncvar_dim(nc)[[varid]]]
      unlimDim = names(nc$dim)[which(ncdim_isUnlim(nc))]
      ncNew = nc_create(filename=output, vars=nc$var[[varid]])
      start = rep(1, length(isUnlim))
      refSize = nc$var[[varid]]$size[which(!isUnlim)]
    }
    ncSize = nc$var[[varid]]$size
    if(!identical(ncSize[which(!isUnlim)], refSize))
      stop("File dimensions doesn't match.")
    count = ncSize*isUnlim -1*!isUnlim
    # add values to varid
    ncvar_put(ncNew, varid, ncvar_get(nc, varid, collapse_degen=FALSE),
              start=start, count=ncSize)
    # add values to unlimited dimension
    ncvar_put(ncNew, varid=unlimDim, ncvar_get(nc, varid=unlimDim),
              start=start[which(isUnlim)], count=ncSize[which(isUnlim)])
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
#' @param newvarid New name for varid in the output file.
#' @param compression Compression level for the new variable (forces ncdf v4).
#' @param force_v4 Logical. Should the resulting file be ncdf v4?
#' @param ... the dimensions and bounds of values to subset.
#' @param ignore.case Logical. Ignore case when matching the dimensions?
#'
#' @return
#' @export
#'
#' @examples
nc_subset = function(filename, varid, output, newvarid, compression,
                     force_v4=FALSE, ..., ignore.case=FALSE) {

  bounds = list(...)
  if(isTRUE(ignore.case)) names(bounds) = tolower(names(bounds))
  nc = nc_open(filename)

  if(missing(varid)) {
    msg = sprintf("ncdf file has more than one variable (%s), argument 'varid' must be specified.",
                  paste(sQuote(names(nc$var)), collapse=", "))
    if(length(nc$var)>1) stop(msg)
    varid = names(nc$var)[1]
  }

  if(missing(newvarid)) newvarid = varid

  dims = ncvar_dim(nc, varid, value=TRUE)
  dimNames = names(dims)

  if(isTRUE(ignore.case)) names(dims) = tolower(names(dims))

  check = any(names(bounds) %in% names(dims), na.rm=TRUE)

  if(!isTRUE(check)) {
    warning("Dimensions to subset don't match, nothing to do.")
    return(invisible())
  }

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
                   dimNames) # keep original names

  start = setNames(sapply(names(dims),
                          function(x) .getIndex(dims[[x]], bounds[[x]], FUN=min, default=1)),
                   names(dims))

  count = setNames(sapply(names(dims),
                          function(x) .getIndex(dims[[x]], bounds[[x]], FUN=length, default=-1)),
                   names(dims))

  x  = ncvar_get(nc, varid, collapse_degen=FALSE, start=start, count=count)

  newVar = nc$var[[varid]]
  newVar$size = dim(x)
  newVar$name = newvarid
  if(!missing(compression)) newVar$compression = compression
  newVar$chunksizes = NA

  .modifyDim = function(x, dim, index) {
    if(isTRUE(index[[x]])) return(dim[[x]])
    dim[[x]]$size = length(index[[x]])
    dim[[x]]$len = length(index[[x]])
    dim[[x]]$vals = dim[[x]]$vals[index[[x]]]
    return(dim[[x]])
  }

  newVar$dim = lapply(names(nc$dim), FUN=.modifyDim, dim=nc$dim, index=index)
  newVar$dim = newVar$dim[newVar$dimids + 1]

  ncNew = nc_create(filename=output, vars=newVar, force_v4=force_v4)
  ncvar_put(ncNew, newvarid, x)

  globalAtt = ncatt_get(nc, varid=0)

  xcall = paste(gsub(x=gsub(x=capture.output(match.call()),
                     pattern="^[ ]*", replacement=""), pattern="\"",
               replacement="'"), collapse="")

  globalAtt$history = sprintf("%s: %s [nctools version %s, %s]",
                              date(), xcall, packageVersion("nctools"), R.version.string)

  # copy global attributes from original nc file.
  ncatt_put_all(ncNew, varid=0, attval=globalAtt)

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
  on.exit(nc_close(nc))

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

  renameFlag = file.rename(outputTemp, output)

  return(invisible(newVars))

}


#' Apply Functions Over Dimensions of a netCDF variable.
#'
#' @param filename Name of the existing netCDF file to be opened.
#' @param varid What variable to read the data from. Can be a string with the
#' name of the variable or an object of class ncvar4. If set to NA,
#' the function will determine if there is only one variable in the file and,
#' if so, read from that, but if there are multiple variables in the file, an error is generated.
#' @param MARGIN a vector giving the dimensions which the function will be applied over.
#' It can be a character vector selecting dimension names.
#' @param FUN the function to be applied
#' @param ... optional arguments to FUN.
#' @param output Name of the file to save results.
#' @param drop Logical. Drop degenered dimensions (i.e. dimensions of length 1)? Not implemented.
#' @param newdim the values to be assigned a the dimension resulting from the
#' application of FUN.
#' @param name new name of the resulting variable. If NULL (by default), the original name is kept.
#' @param longname long name of the resulting variable.
#' @param units units of the resulting variable.
#' @param compression If set to an integer between 1 (least compression) and 9 (most compression), this enables compression for the variable as it is written to the file. Turning compression on forces the created file to be in netcdf version 4 format, which will not be compatible with older software that only reads netcdf version 3 files.
#' @param verbose Print debugging information.
#' @param force_v4 If TRUE, then the created output file will always be in netcdf-4 format (which supports more features, but cannot be read by version 3 of the netcdf library). If FALSE, then the file is created in netcdf version 3 format UNLESS the user has requested features that require version 4. Deafult is TRUE.
#' @param ignore.case If TRUE, ignore case in matching dimension names and MARGIN. Default is FALSE.
#'
#' @return
#' @export
#'
#' @examples
nc_apply = function(filename, varid, MARGIN, FUN, ..., output=NULL, drop=TRUE,
                    newdim = NULL, name=NULL, longname=NULL, units=NULL,
                    compression=NA, verbose=FALSE, force_v4=TRUE,
                    ignore.case=FALSE) {


  funName = deparse(substitute(FUN))

  if(is.null(output)) stop("You must specify an 'output'file.")

  if(file.exists(output)) {
    oTest = file.remove(output)
    if(!oTest) stop(sprintf("Cannot write on %s.", output))
  }

  FUN = match.fun(FUN)

  nc = nc_open(filename)
  on.exit(nc_close(nc))

  if(is.na(varid)) {
    if(length(nc$var)==1) varid = nc$var[[1]]$name
    msg = sprintf("Several variables found in %s, must specify 'varid'.", filename)
    if(length(nc$var)>1) stop(msg)
  }

  if(inherits(varid, "ncvar4")) varid = varid$name

  X = ncvar_get(nc, varid, collapse_degen = FALSE)

  dn = ncvar_dim(nc, varid, value=TRUE)
  dnn = if(isTRUE(ignore.case)) tolower(names(dn)) else names(dn)

  if (is.character(MARGIN)) {
    if(isTRUE(ignore.case)) MARGIN = tolower(MARGIN)
    MARGIN <- match(MARGIN, dnn)
    if(anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  }

  Y = apply(X=X, MARGIN = MARGIN, FUN = FUN, ...)

  lans = length(Y)/prod(dim(X)[MARGIN]) # length of the answer (FUN)

  if(is.null(newdim)) {

    if(lans>1) {
      vals = suppressWarnings(as.numeric(dimnames(Y)[[1]]))
      newdim = if(all(!is.na(vals))) vals else seq_len(lans)
    } else {
      newdim = seq_len(lans)
    }

  } # make newdim from answer

  if(length(newdim)!=lans) {
    msg = sprintf("The length of newdim (%s) doesn't match answer length (%s), ignoring newdim.",
                  length(newdim), lans)
    warning(msg)
    newdim = seq_len(lans)
  }
  if(!is.numeric(newdim)) {
    warning("The argument newdim must be numeric, ignoring newdim.")
    newdim = seq_len(lans)
  }

  oldVar = nc$var[[varid]]
  newDim = nc$dim[(oldVar$dimids + 1)[MARGIN]]

  if(!is.na(compression)) oldVar$compression = compression

  if(lans>1) {

    Y = aperm(Y, perm = c(seq_along(dim(Y))[-1], 1))

    ansDim = ncdim_def(name=funName, units="", vals=newdim)
    newDim = c(newDim, list(ansDim))
    names(newDim)[length(MARGIN)+1] = funName

  }

  # drop for tomorrow

  xlongname = sprintf("%s of %s", funName, oldVar$longname)

  varLongname = if(is.null(longname)) xlongname else longname
  varName  = if(is.null(name)) oldVar$name else name
  varUnits = if(is.null(units)) oldVar$units else units

  newVar = ncvar_def(name=varName, units = varUnits,
                     missval = oldVar$missval, dim = newDim,
                     longname = varLongname, prec = oldVar$prec,
                     compression = oldVar$compression)

  ncNew = nc_create(filename=output, vars=newVar)
  ncvar_put(nc=ncNew, varid=varName, vals=Y)
  nc_close(ncNew)

  return(invisible(output))

}


# Extra tools -------------------------------------------------------------

#' Data output in ncdf format
#'
#' @param x An array to write to a ncdf file.
#' @param filename The file to write.
#' @param varid The name of the variable in the ncdf file.
#' @param dim A list with the values of the dimensions. Names are taken from the list.
#' @param longname The longname for the variable to be created.
#' @param units The units for the variable to be created.
#' @param prec The precision for the variable to be created.
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

