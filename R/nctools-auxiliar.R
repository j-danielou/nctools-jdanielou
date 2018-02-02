
# Generic ncdf auxiliar functions -----------------------------------------


# Variables ---------------------------------------------------------------


#' Get dimension names for each variable in a ncdf file
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param varid
#' @param value
#'
#' @return
#' @export
#'
#' @examples
ncvar_dim = function(nc, varid=NULL, value=FALSE) {

  if (!inherits(nc, "ncdf4"))
    stop("first argument (nc) is not of class ncdf4!")

  if(isTRUE(value)) {
    .getDimValues = function(x) {
      vals = stats::setNames(lapply(x, function(x) x$vals),
                             sapply(x, function(x) x$name))
    }

    out = lapply(nc$var, function(x) x$dim)
    out = lapply(out, .getDimValues)
  } else {
    out = lapply(nc$var, function(x) names(nc$dim)[x$dimids+1])
  }
  if(!is.null(varid)) out = out[[varid]]
  return(out)
}


#' Get the size of the variables in a ncdf object.
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param varid
#'
#' @return
#' @export
#'
#' @examples
ncvar_size = function(nc, varid=NULL) {

  if (!inherits(nc, "ncdf4"))
    stop("first argument (nc) is not of class ncdf4!")

  out = lapply(nc$var, function(x) x$size)
  if(!is.null(out)) return(out[[varid]])
  return(out)
}



# Dimensions --------------------------------------------------------------


#' Get the dimensions of the variables in a ncdf object.
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#'
#' @return A named list with the dimensions of the variables.
#' @export
#'
#' @examples
ncdim_size = function(nc) {

  if (!inherits(nc, "ncdf4"))
    stop("first argument (nc) is not of class ncdf4!")

  lapply(nc$dim, function(x) x$len)
}


#' Is the dimention unlimited?
#'
#' @param nc  An open connection to a netCDF file as in nc_open(file).
#'
#' @return
#' @export
#'
#' @examples
ncdim_isUnlim = function(nc) {

  if (!inherits(nc, "ncdf4"))
    stop("first argument (nc) is not of class ncdf4!")

  sapply(nc$dim, function(x) x$unlim)
}



# Attributes --------------------------------------------------------------

#' Get all variable's attributes from a ncdf object
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param type A character to choose attributes for variables ("var")
#' or dimensions ("dim").
#' @return A list with all the atributes.
#' @export
#'
#' @examples
ncatt_get_all = function(nc, type=c("var", "dim")) {

  if (!inherits(nc, "ncdf4"))
    stop("first argument (nc) is not of class ncdf4!")

  type = match.arg(type)
  vars = names(nc[[type]])
  names(vars) = vars
  atts = lapply(vars, FUN=function(var) ncatt_get(nc, var))
  return(atts)
}

