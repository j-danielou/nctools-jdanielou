## TO_DO

### URGENT
uniform handling of files (temporal files, overwrite, etc.)
  check when filename==output, so temp file is written. Always?
  nctools::ncdim_unlim: it failes when files are open?
  ncrcat is closing files?

check attributes are copied

consistent check for 'output' argument, specially removing if already in existence.

for nc_subset and related, check if something will be extracted, error otherwise.
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

