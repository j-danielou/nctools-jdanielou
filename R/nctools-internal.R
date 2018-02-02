

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

  if(length(oldname)!=length(newname))
    stop("oldname and newname must have the same length.")

  nc = nc_open(filename)
  on.exit(try(nc_close(nc), silent=TRUE))

  if(missing(output)) output = filename

  for(i in seq_along(oldname)) {

    nc$var = lapply(nc$var, FUN=.replaceInDim, dim=oldname[i],
               id="name", value=newname[i])

  }

  tmp = paste(output, "temp", sep=".")

  ncNew = nc_create(filename = tmp, vars = nc$var, verbose=verbose)
  on.exit(try(nc_close(ncNew), silent=TRUE), add=TRUE)
  for(varid in names(nc$var))
    ncvar_put(ncNew, varid=varid, vals=ncvar_get(nc, varid=varid),
              verbose=verbose)

  nc_close(ncNew)
  nc_close(nc)

  file.rename(tmp, output)
  on.exit(try(file.remove(tmp), silent=TRUE), add=TRUE)

  return(invisible(output))

}

.nc_renameVar = function(filename, oldname, newname, output, verbose=FALSE) {

  if(length(oldname)!=length(newname))
    stop("oldname and newname must have the same length.")

  if(missing(output)) output = filename

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

  nc_close(nc)

  return(invisible(output))

}
