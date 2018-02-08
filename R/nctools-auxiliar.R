
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


#' Get the compression of the variables in a ncdf object.
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param varid Either the name of the variable or an ncvar object.
#'
#' @return
#' @export
#'
#' @examples
ncvar_compression = function(nc, varid=NA) {
  if(is.na(varid)) return(sapply(nc$var, FUN=.getCompression))
  if(inherits(varid, "ncvar4")) varid = varid$name
  if(!is.character(varid))
    stop("varid must be a string or an object of class 'ncvar4'.")
  return(.getCompression(nc$var[[varid]]))
}

#' Change the compression of the variables in a ncdf object.
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param varid Either the name of the variable or an ncvar object indicating whose compression value will be changed.
#' If NA, all the variables will be changed to the new compression.
#' @param compression The new compression value.
#'
#' @return
#' @export
#'
#' @examples
ncvar_change_compression = function(nc, varid=NA, compression) {
  if(is.na(varid)) {
    nc$var = lapply(nc$var, FUN=.setCompression, compression=compression)
    return(nc)
  }
  if(inherits(varid, "ncvar4")) varid = varid$name
  if(!is.character(varid))
    stop("varid must be a string or an object of class 'ncvar4'.")
  nc$var[[varid]] = .setCompression(nc$var[[varid]], compression=compression)
  return(nc)
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

#' Put several attributes to a ncdf object
#'
#' @param nc An open connection to a netCDF file as in nc_open(file).
#' @param varid The variable whose attribute is to be written. Can be a
#' character string with the variable's name, an object of class ncvar4,
#' or an id contained in the "id" field of a ncvar object. As a special case,
#' if varid==0, then a global attribute is written instead of a variable's
#' attribute.
#' @param attname The names of the attributes.
#' @param attval The values of the attributes.
#' @param prec Precision to write the attribute.
#' @param verbose Can be set to TRUE if additional information is desired while
#' the attribute is being created.
#' @param definemode See 'definemode' in ncatt_put for details.
#'
#' @return
#' @export
#'
#' @examples
ncatt_put_all = function(nc, varid, attname, attval,
                         prec=NA, verbose=FALSE, definemode=FALSE) {


  if(missing(attname) & missing(attval))
    stop("You must provide values and names for the attributes.")

  if(missing(attname) & !is.null(names(attval)))
    attname = names(attval)

  if(missing(attval) & is.null(names(attname)))
    stop("You must provide values and names for the attributes.")

  if(missing(attval) & !is.null(names(attname))) {
    attval = attname
    attname = names(attname)
  }

  if(any(is.na(attname)))
    stop("Missing names for attributes are not allowed.")

  if(length(attval)!=length(attname))
    stop("An equal number of names and values for attributes must be provided.")

  names(attval) = attname

  .ncatt_put = function(attname, attval, nc, varid, prec, verbose, definemode) {
    ncatt_put(nc, varid=varid, attname=attname, attval=attval[[attname]],
              prec=prec, verbose=verbose, definemode=definemode)
  }

  lapply(attname, FUN=.ncatt_put, attval=attval, nc=nc, varid=varid,
         prec=prec, verbose=verbose, definemode=definemode)

  return(invisible())

}
