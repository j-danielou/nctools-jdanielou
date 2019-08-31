
#' @export
nc_loc = function(filename, varid, MARGIN, loc, output=NULL, drop=TRUE,
                  newdim = NULL, name=NULL, longname=NULL, units=NULL,
                  compression=NA, verbose=FALSE, force_v4=TRUE,
                  ignore.case=FALSE) {


  if(is.null(output)) stop("You must specify an 'output'file.")

  if(file.exists(output)) {
    oTest = file.remove(output)
    if(!oTest) stop(sprintf("Cannot write on %s.", output))
  }

  nc = nc_open(filename)
  on.exit(nc_close(nc))

  varid = .checkVarid(varid=varid, nc=nc)

  dn = ncvar_dim(nc, varid, value=TRUE)
  dnn = if(isTRUE(ignore.case)) tolower(names(dn)) else names(dn)

  if(length(MARGIN)!=1) stop("MARGIN must be of length one.")

  if (is.character(MARGIN)) {
    if(isTRUE(ignore.case)) MARGIN = tolower(MARGIN)
    MARGIN = match(MARGIN, dnn)
    if(anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  }

  depth = ncvar_get(nc, varid=dnn[MARGIN])

  dims = nc$var[[varid]]$size
  dims[MARGIN] = dims[MARGIN]-1
  ind0 = c(list(X=ncvar_get(nc, varid, collapse_degen = FALSE) - loc, drop=FALSE),
           lapply(dims, seq_len))

  x0 = do.call("[", ind0)
  ind0[[MARGIN+2]] = ind0[[MARGIN+2]]+1L
  x1 = do.call("[", ind0)
  x0 = sign(x0)*sign(x1)

  ind = apply(x0, -MARGIN, FUN=function(x) which(x<1)[1])

  D1 = depth[ind]
  D2 = depth[ind+1]

  iList = append(lapply(dims[-MARGIN], seq_len), MARGIN-1, values=NA_integer_)

  index = as.matrix(do.call(expand.grid, iList))
  index[, MARGIN] = ind
  x1 = ind0$X[index]
  index[, MARGIN] = ind + 1
  x2 = ind0$X[index]

  Dx = (-x1*(D2-D1))/(x2-x1) + D1
  dim(Dx) = dims[-MARGIN]

  oldVar = nc$var[[varid]]
  newDim = nc$dim[(oldVar$dimids + 1)[-MARGIN]]
  thisDim = nc$dim[(oldVar$dimids + 1)[MARGIN]][[1]]

  if(!is.na(compression)) oldVar$compression = compression

  xlongname = sprintf("%s of %s=%s %s", dnn[MARGIN], oldVar$longname, loc, oldVar$units)

  varLongname = if(is.null(longname)) xlongname else longname
  varName  = if(is.null(name)) dnn[MARGIN] else name
  varUnits = thisDim$units

  newVar = ncvar_def(name=varName, units = varUnits,
                     missval = oldVar$missval, dim = newDim,
                     longname = varLongname, prec = oldVar$prec,
                     compression = oldVar$compression)

  ncNew = nc_create(filename=output, vars=newVar)
  ncvar_put(ncNew, varName, Dx)
  nc_close(ncNew)

  return(invisible(output))

}
