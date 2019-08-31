#' @export
write.ncdf = function(x, filename, varid, dim, longname, units, prec="float",
                      missval=-9999, compression=9, chunksizes=NA, verbose=FALSE,
                      dim.units, dim.longname, unlim=FALSE, ...) {
  .Deprecated("write_ncdf")
  write_ncdf.default(x=x, filename=filename, varid=varid, dim=dim, longname=longname,
                     units=units, prec=prec, missval=missval, compression=compression,
                     chunksizes=chunksizes, verbose=verbose, dim.units=dim.units,
                     dim.longname=dim.longname, unlim=unlim, ...)
}


