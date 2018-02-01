#' #' Subset a ncdf variable
#' #'
#' #' @param filename The name of the ncdf file to subset.
#' #' @param varid The name of the variable to subset. If missing and only one variable in the file, that one is used.
#' #' @param output The name of the ncdf output file to create.
#' #' @param ...
#' #'
#' #' @return
#' #'
#' #' @examples
#' nc_subset2 = function(filename, varid, output, ...) {
#'
#'   bounds = list(...)
#'   nc = nc_open(filename)
#'
#'   if(missing(varid)|is.na(varid)) {
#'     msg = sprintf("ncdf file has more than one variable (%s), argument 'varid' must be specified.",
#'                   paste(sQuote(names(nc$var)), collapse=", "))
#'     if(length(nc$var)>1) stop(msg)
#'     varid = names(nc$var)[1]
#'   }
#'   x  = ncvar_get(nc, varid, collapse_degen=FALSE)
#'   dims = ncvar_dim(nc, varid, value=TRUE)
#'   dimnames(x) = dims
#'
#'   .isInside = function(x, bound) {
#'     # longitude in a torus
#'     if(is.null(bound)) return(TRUE)
#'     if(diff(bound)<0) stop("Upper bound is lower than lower bound.")
#'     out = which((x>=bound[1]) & (x<=bound[2]))
#'     return(out)
#'   }
#'
#'   index = setNames(lapply(names(dims),
#'                           function(x) .isInside(dims[[x]], bounds[[x]])),
#'                    names(dims))
#'
#'   x = do.call("[", c(list(x), index))
#'
#'   newVar = nc$var[[varid]]
#'   newVar$size = dim(x)
#'   newVar$chunksizes = NA
#'
#'   .modifyDim = function(x, dim, index) {
#'     if(isTRUE(index[[x]])) return(dim[[x]])
#'     dim[[x]]$size = length(index[[x]])
#'     dim[[x]]$len = length(index[[x]])
#'     dim[[x]]$vals = dim[[x]]$vals[index[[x]]]
#'     return(dim[[x]])
#'   }
#'
#'   newVar$dim = lapply(names(nc$dim), FUN=.modifyDim, dim=nc$dim, index=index)
#'
#'   ncNew = nc_create(filename=output, vars=newVar)
#'   ncvar_put(ncNew, varid, x)
#'   nc_close(ncNew)
#'
#'   nc_close(nc)
#'
#'   return(invisible(output))
#' }
