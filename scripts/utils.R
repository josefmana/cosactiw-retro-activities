# This script includes utility functions used throughout analyses for this project.
# Ought to be loaded before analysis starts.

# return rounded number as a character
rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )

# return CI as [lower_CI, upper_CI]
ciprint <- function(x, d = 2) paste0( "[", paste( rprint(x, d), collapse = ", "), "]" )

# strip leading zero and return result, print < .001 if the result is too small a number
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )

# calculate and return mean ± SD ignoring NAs
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )
