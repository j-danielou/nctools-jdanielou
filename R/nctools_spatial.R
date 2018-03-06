
#' Check and correct Longitude values given the position of the Prime Meridian (International Reference Meridian)
#'
#' @param x numeric vector of longitude values.
#' @param primeMeridian position of the Prime Meridian. Use "center" for [-180,180]
#' range values and "left" for [0,360] range (Pacific centered).
#' @param sort sort longitude values?
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
checkLongitude = function(x, primeMeridian="center", sort=FALSE, ...) {

  .longitude2Center = function(x, ...) {
    if (!any(x > 180))
      return(x)
    x[x > 180] = x[x > 180] - 360
    return(x)
  }

  .longitude2Left = function(x, ...) {
    if (!any(x < 0))
      return(x)
    x[x < 0] = x[x < 0] + 360
    return(x)
  }

  out = switch(primeMeridian,
               center = .longitude2Center(x, ...),
               left = .longitude2Left(x, ...))

  if(isTRUE(sort)) out = sort(out)

  return(out)
}


#' Find the Prime Meridian
#'
#' @param x numeric vector of longitude values.
#'
#' @return
#' @export
#'
#' @examples
findPrimeMeridian = function(x) {
  if(any(x<0)) return("center")
  if(any(x>180)) return("left")
  warning("Indeterminate Prime Meridian from longitude values.")
  return(NULL)
}

