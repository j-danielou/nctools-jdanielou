

.replaceInDim = function(x, dim, id, value) {
  .repInDim = function(x, dim, id, value) {
    if(x$name != dim) return(x)
    x[[id]] = value
    return(x)
  }
  x$dim = lapply(x$dim, FUN=.repInDim, dim=dim,
                 id=id, value=value)
  return(x)
}

.nc_renameDim = function(filename, oldname, newname, output, verbose=FALSE) {

  tmp = paste(output, "temp", sep=".")
  on.exit(if(file.exists(tmp)) file.remove(tmp))

  if(length(oldname)!=length(newname))
    stop("oldname and newname must have the same length.")

  nc = nc_open(filename)
  is_v4 = grepl(x=nc$format, pattern="NETCDF4")

  if(is_v4) {
    varids = names(nc$var)
    for(varid in varids) {
      if(is.na(ncvar_compression(nc, varid)))
        nc = ncvar_change_compression(nc, varid, compression = 9)
    }
  }

  x = nc$var

  for(i in seq_along(oldname)) {
    x = lapply(x, FUN=.replaceInDim, dim=oldname[i],
               id="name", value=newname[i])
  }


  ncNew = nc_create(filename = tmp, vars = x, verbose=verbose,
                    force_v4 = is_v4)

  for(varid in names(x))
    ncvar_put(ncNew, varid=varid, vals=ncvar_get(nc, varid=varid),
              verbose=verbose)

  # copy global attributes from original nc file.
  ncatt_put_all(ncNew, varid=0, attval=ncatt_get(nc, varid=0))

  nc_close(ncNew)
  nc_close(nc)

  renamed = file.rename(tmp, output)

  if(!renamed) {
    file.remove(output)
    renamed = file.rename(tmp, output)
    if(!renamed) {
      file.remove(tmp)
      stop(sprintf("Couldn't write %s.", output))
    }
  }

  return(invisible(output))

}

.nc_renameVar = function(filename, oldname, newname, output, verbose=FALSE) {

  if(length(oldname)!=length(newname))
    stop("oldname and newname must have the same length.")

  if(output!=filename) {
    file.copy(from=filename, to=output, overwrite = TRUE)
    filename = output
  }

  nc = nc_open(filename, write=TRUE)
  on.exit(try(nc_close(nc), silent = TRUE))

  for(i in seq_along(oldname)) {

    nc = ncvar_rename(nc, old_varname = oldname[i], new_varname = newname[i],
                      verbose=verbose)

  }

  # nc_close(nc)

  return(invisible(output))

}

.getCompression = function(x) return(x$compression)

.setCompression = function(x, compression) {
  x$compression = compression
  x$chunksizes = NA
  return(x)
}



# Argument checking -------------------------------------------------------


.checkVarid = function(varid, nc) {

  if(missing(varid)) varid = NA

  if(is.na(varid)) {
    if(length(nc$var)==1) varid = nc$var[[1]]$name
    msg = sprintf("Several variables found in %s, must specify 'varid'.", nc$filename)
    if(length(nc$var)>1) stop(msg)
  }

  if(inherits(varid, "ncvar4")) varid = varid$name

  if(!is.character(varid))
    stop("varid must be a string or an object of class 'ncvar4'.")

  varExists = varid %in% names(nc$var)

  msg = sprintf("Variable '%s' not found in '%s'.", varid, nc$filename)
  if(!varExists) stop(msg)

  return(varid)

}

# match and array with a dimsize vector.
.getDimensions = function(x, dimsize) {
  mydim = dim(x)
  out = NA_real_*numeric(length(mydim))
  for(i in seq_along(mydim)) {
    ind = which(dimsize %in% mydim[i])
    if(length(ind)<1) stop("Array incompatible with dimension units.")
    if(length(ind)>1) ind = min(ind)
    out[i] = ind
    dimsize[seq_len(ind)] = -1
  }
  return(out)
}
