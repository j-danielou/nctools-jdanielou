
# Generic ncdf auxiliar functions -----------------------------------------


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
  type = match.arg(type)
  vars = names(nc[[type]])
  names(vars) = vars
  atts = lapply(vars, FUN=function(var) ncatt_get(nc, var))
  return(atts)
}


# Size of the variables and dimensions ------------------------------------

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
  out = lapply(nc$var, function(x) x$size)
  if(!is.null(out)) return(out[[varid]])
  return(out)
}

#' Get the dimensions of the variables in a ncdf object.
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#'
#' @return A named list with the dimensions of the variables.
#' @export
#'
#' @examples
ncdim_size = function(nc) {
  lapply(nc$dim, function(x) x$len)
}


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

#' Is the dimention unlimited?
#'
#' @param nc  An open connection to a netCDF file as in nc_open(file).
#'
#' @return
#' @export
#'
#' @examples
ncdim_isUnlim = function(nc) {
  sapply(nc$dim, function(x) x$unlim)
}

