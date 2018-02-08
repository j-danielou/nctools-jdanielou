library(nctools)

ovar = c("O2_ZS", "DEPTH1_1_bnds")
nvar = c("o2", "toDelete")
odim = c("DEPTH1_1", "LATITUDE", "LALA", "TIME", "time")
ndim = c("depth", "lat", "pooo", "time", "time")

filename = "test.nc"
file.copy("master.nc", filename)

# test 1
cat("Test 1a ---------------------- \n\n")
nc_rename(filename, ovar, nvar, "test1a.nc", verbose=TRUE) # other file
nc_rename(filename, ovar, nvar, "test1a.nc", verbose=TRUE, overwrite = TRUE) # other file
cat("Test 1b ---------------------- \n\n")
file.copy("master.nc", "test1b.nc", overwrite=TRUE)
nc_rename("test1b.nc", ovar, nvar, verbose=TRUE) # same file
nc_rename("test1b.nc", ovar, nvar, verbose=TRUE, overwrite = TRUE) # same file

# test 2
cat("Test 2a ---------------------- \n\n")
nc_rename(filename, odim, ndim, "test2a.nc", verbose=TRUE) # other file
nc_rename(filename, odim, ndim, "test2a.nc", verbose=TRUE, overwrite = TRUE) # other file
cat("Test 2b ---------------------- \n\n")
file.copy(filename, "test2b.nc", overwrite = TRUE)
nc_rename("test2b.nc", odim, ndim, verbose=TRUE) # same file
nc_rename("test2b.nc", odim, ndim, verbose=TRUE, overwrite = TRUE) # same file

# test 3
cat("Test 3a ---------------------- \n\n")
nc_rename(filename, c(ovar, odim), c(nvar, ndim), "test3a.nc", verbose=TRUE) # other file
nc_rename(filename, c(ovar, odim), c(nvar, ndim), "test3a.nc", verbose=TRUE, overwrite = TRUE) # other file
cat("Test 3b ---------------------- \n\n")
file.copy(filename, "test3b.nc", overwrite = TRUE)
nc_rename("test3b.nc", c(ovar, odim), c(nvar, ndim), verbose=TRUE) # same file
nc_rename("test3b.nc", c(ovar, odim), c(nvar, ndim), verbose=TRUE, overwrite = TRUE) # same file


