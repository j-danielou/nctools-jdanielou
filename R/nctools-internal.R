

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

