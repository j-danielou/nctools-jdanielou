
#' Calculate sample quantiles over the dimensions of a ncdf variable
#'
#' @param filename Name of the existing netCDF file to be opened.
#' @param varid What variable to read the data from. Can be a string with the
#' name of the variable or an object of class ncvar4. If set to NA,
#' the function will determine if there is only one variable in the file and,
#' if so, read from that, but if there are multiple variables in the file, an error is generated.
#' @param MARGIN a vector giving the dimensions which the function will be applied over.
#' It can be a character vector selecting dimension names. By default, the first two dimensions are taken.
#' @param na.rm ogical; if TRUE, any NA and NaN's are removed before the quantiles are computed.
#' @param probs numeric vector of probabilities with values in [0,1]. Passed to stats::quantile.
#' @param output Name of the file to save results.
#' @param drop Logical. Drop degenered dimensions (i.e. dimensions of length 1)? Not implemented.
#' @param compression If set to an integer between 1 (least compression) and 9 (most compression), this enables compression for the variable as it is written to the file. Turning compression on forces the created file to be in netcdf version 4 format, which will not be compatible with older software that only reads netcdf version 3 files.
#' @param verbose Print debugging information.
#' @param force_v4 If TRUE, then the created output file will always be in netcdf-4 format (which supports more features, but cannot be read by version 3 of the netcdf library). If FALSE, then the file is created in netcdf version 3 format UNLESS the user has requested features that require version 4. Deafult is TRUE.
#' @param ignore.case If TRUE, ignore case in matching dimension names and MARGIN. Default is FALSE.
#'
#' @return
#' @export
#'
#' @examples
nc_quantile = function(filename, varid, MARGIN=c(1,2), na.rm=TRUE, probs=c(0, 0.5, 1),
                       output=NULL, drop=TRUE, compression=NA, verbose=FALSE,
                       force_v4=TRUE, ignore.case=FALSE) {

  if(!is.numeric(probs)) stop("Argument 'probs' must be numeric.")

  nc_apply(filename=filename, varid=varid, MARGIN=MARGIN, FUN=quantile,
           na.rm=na.rm, probs=probs, output=output, drop=drop, newdim = probs,
           name=NULL, longname=NULL, units=NULL, compression=compression, verbose=verbose,
           force_v4=force_v4, ignore.case=ignore.case)

}


