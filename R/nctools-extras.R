

#' Change the Prime Meridian of a variable in a ncdf file.
#'
#' @param filename The filename of the original ncdf file.
#' @param output The name of the output file.
#' @param varid The name of the variable to extract.
#' @param primeMeridian position of the Prime Meridian. Use "center" for [-180,180]
#' range values and "left" for [0,360] range (Pacific centered).
#' @param verbose If TRUE, run verbosely.
#' @param overwrite overwrite output file if already exists?
#' @param compression
#' @param mem.limit Maximum RAM size to be used per iteration.
#' If the total size of the variable exceeds this value, the reshape is done iteratively.
#' By default is 3072MB (3GB).
#'
#' @return
#' @export
#'
#' @examples
nc_changePrimeMeridian = function(filename, output, varid=NA, MARGIN=1, primeMeridian="center",
                                  verbose=FALSE, overwrite=FALSE, compression=NA,
                                  mem.limit=3072, ignore.case=FALSE) {

  if(missing(output) & !isTRUE(overwrite))
    stop("output file is missing. Set 'overwrite' to TRUE to make changes in the original file.")

  if(missing(output)) output = filename

  if(file.exists(output) & !isTRUE(overwrite))
    stop("output file already exists. Set 'overwrite' to TRUE.")

  tmp = paste(output, ".temp", sep="")

  nc = nc_open(filename)
  on.exit(nc_close(nc))

  varid = .checkVarid(varid=varid, nc=nc)

  dn = ncvar_dim(nc, varid, value=TRUE)
  dnn = if(isTRUE(ignore.case)) tolower(names(dn)) else names(dn)

  if(length(MARGIN)!=1) stop("MARGIN must be of length one (dim of 'longitude').")

  if (is.character(MARGIN)) {
    if(isTRUE(ignore.case)) MARGIN = tolower(MARGIN)
    MARGIN = match(MARGIN, dnn)
    if(anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  }

  ivar = nc$var[[varid]]
  lon = ivar$dim[[MARGIN]]$vals
  ndim = length(ivar$dim)
  if(ndim<2) stop("Data must have at least two dimensions!")

  cellLimit = (mem.limit*2^20)/8

  bigData = (prod(ivar$size) > cellLimit) # 3GB by default

  if(bigData) {

    npiece = floor(prod(ivar$size)/cellLimit)
    useDim = setdiff(which(ivar$size >= npiece), MARGIN)
    useDim = which.min(ivar$size[useDim])
    itDim  = max(ceiling(ivar$size[useDim]/npiece), 1)

    starts = seq(from=1, to=ivar$size[useDim], by=itDim)
    counts = diff(c(starts-1, ivar$size[useDim]))
    npiece = length(starts)

    start = rep(1, ndim)
    count = rep(-1, ndim)

  }

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
  ivar$dim[[MARGIN]]$vals = newlon[ind]
  if(!is.na(compression)) ivar$compression = compression
  ivar$chunksizes = NA

  ncNew = nc_create(filename=tmp, vars=ivar)
  on.exit(if(file.exists(tmp)) file.remove(tmp))

  if(!bigData) {

    newvar = c(list(x=ncvar_get(nc, varid, collapse_degen=FALSE),
                    drop=FALSE), rep(TRUE, ndim))
    newvar[[MARGIN+2]] = ind
    newvar = do.call('[', newvar)
    ncvar_put(ncNew, varid, newvar)

  } else {

    message("Using big data method.")
    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, 0)


    for(i in seq_len(npiece)) {

      start[useDim] = starts[i]
      count[useDim] = counts[i]

      newvar = c(list(x=ncvar_get(nc, varid, collapse_degen=FALSE,
                                  start=start, count=count),
                      drop=FALSE), rep(TRUE, ndim))
      newvar[[MARGIN+2]] = ind
      newvar = do.call('[', newvar)
      invisible(gc())
      ncvar_put(ncNew, varid, newvar, start=start, count=count)
      nc_sync(ncNew)

      pb = txtProgressBar(style=3)
      setTxtProgressBar(pb, i/npiece)

    }

  }

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
