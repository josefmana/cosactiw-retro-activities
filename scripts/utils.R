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

# plot posterior predictions
ppc_plot <-
  
  function( fit,
            data,
            y, x,
            labs = list(NULL, NULL),
            meth = "ppc_dens_overlay_grouped",
            draws = sample(1:4000,100),
            stat = "mean"
          ) {
    
    pp_check(
      
      object = data[ , y],
      yrep = posterior_predict(fit, newdata = data)[draws, ],
      fun = meth,
      stat = stat,
      group = data[ , x]
      
    ) +
      
      labs(
        title = labs[[1]],
        subtitle = labs[[2]]
      ) +
      
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5)
      )
    
  }

  
