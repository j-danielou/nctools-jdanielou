## TO_DO

### URGENT
one handler of varid (NA, exceptions) for all the functions. 
uniform handling of files (temporal files, overwrite, etc.)
  check when filename==output, so temp file is written. Always?
  nctools::ncdim_unlim: it failes when files are open?
  ncrcat is closing files?

progress bar for nc_rcat and all functions

conversion from HDF4 (if hdf2nc installed.)

### NOT URGENT
change dimension values

newfunctions:

- nc_sd
- nc_fivenums (save 5 variables, varid_mxx)
- group all them in one help page! 

- add basic arithmetics?

nc_subset: example with depth=1 level.

check on regrid

## PIPES

How to use pipes: vectorize functions, return output file %>%

### a simple wrapper to NCO (another package)

nco("ncrcat", files=files, output=, ...)
ncks(...) = function(...) nco("ncks", ...)

input == output, check for --no_tmp_fl

