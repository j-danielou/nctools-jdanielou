#' @export
write_ncdf.default = function(x, filename, varid, dim, longname, units, prec="float",
                              missval=-9999, compression=9, chunksizes=NA, verbose=FALSE,
                              dim.units, dim.longname, unlim=FALSE, global=list(),
                              force_v4=FALSE, ...) {

  if(!is.list(global)) stop("'global' must be a list")

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

  ncNew = nc_create(filename=filename, vars=iVar, force_v4=force_v4, verbose=verbose)
  on.exit(nc_close(ncNew))

  ncvar_put(ncNew, varid=iVar, vals=x, verbose=verbose)

  xcall = paste(gsub(x=gsub(x=capture.output(match.call()),
                            pattern="^[ ]*", replacement=""), pattern="\"",
                     replacement="'"), collapse="")
  globalAtt = global
  globalAtt$history = sprintf("File create on %s: %s [nctools version %s, %s]",
                       date(), xcall, packageVersion("nctools"), R.version.string)
  # create global attributes.
  ncatt_put_all(ncNew, varid=0, attval=globalAtt)

  return(invisible(filename))

}

#' @export
write_ncdf.list = function(x, filename, varid, dim, longname, units, prec="float",
                           missval=-9999, compression=9, chunksizes=NA, verbose=FALSE,
                           dim.units, dim.longname, unlim=FALSE, global=list(),
                           force_v4=FALSE, ...) {

  if(!is.list(global)) stop("'global' must be a list")

  nvar = length(x)

  if(missing(varid)) varid = names(x)
  if(length(varid)!=nvar) stop("One 'varid' per variable must be provided")

  if(missing(dim)) dim = lapply(base::dim(x[[1]]), seq_len)

  if(length(dim)!=length(dim(x[[1]])))
    stop("dim argument does not match data dimension.")

  ind = lapply(x, dim)
  ind = sapply(ind, FUN = identical, y=ind[[1]])
  if(!all(ind)) stop("All arrays to be added to the ncdf file must have the same dimension.")

  if(is.null(names(dim)))
    dim = setNames(dim, paste("dim", seq_along(dim(x[[1]])), sep=""))

  if(missing(longname)) longname = rep("", nvar)
  if(length(longname)==1) longname = rep(longname, nvar)
  if(length(longname)!=nvar) stop("One longname value per variable must be provided.")

  if(missing(units))    units    = rep("", nvar)
  if(length(units)==1) units = rep(units, nvar)
  if(length(units)!=nvar) stop("One units value per variable must be provided.")

  if(length(prec)==1) prec = rep(prec, nvar)
  if(length(prec)!=nvar) stop("One precision value per variable must be provided.")

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

  iVar = list()

  for(i in seq_along(x)) {

    iVar[[i]] = ncvar_def(name=varid[i], units=units[i], dim=dims, prec=prec[i], missval=missval,
                          longname=longname[i], compression=compression, chunksizes=chunksizes, verbose=verbose)

  }

  ncNew = nc_create(filename=filename, vars=iVar, force_v4=force_v4, verbose=verbose)
  on.exit(nc_close(ncNew))

  for(i in seq_along(x)) ncvar_put(ncNew, varid=iVar[[i]], vals=x[[i]], verbose=verbose)

  xcall = paste(gsub(x=gsub(x=capture.output(match.call()),
                            pattern="^[ ]*", replacement=""), pattern="\"",
                     replacement="'"), collapse="")

  globalAtt = global
  globalAtt$history = sprintf("File create on %s: %s [nctools version %s, %s]",
                              date(), xcall, packageVersion("nctools"), R.version.string)
  # create global attributes.
  ncatt_put_all(ncNew, varid=0, attval=globalAtt)


  return(invisible(filename))

}
