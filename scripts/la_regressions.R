# Use this script to describe retrospectively recalled leisure activities with respect to cognitive SA via Bayesian regressions.

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

library(here)
library(tidyverse)
library(brms)
library(bayesplot)

theme_set( theme_bw(base_size = 12) )
source( here("scripts","utils.R") ) # in-house functions

sapply( c("figures","tables","models"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders


# DATA ----

# read them
for ( i in names( readRDS("_data.rds") ) ) assign( i, readRDS("_data.rds")[[i]] )

# prepare a wide data set with all activities and their counts disregarding intensity and seasoness
d1 <- sapply(
  
  map$activity,
  function(i)
    
    sapply(
      unique(act$ID),
      function(j)
        if_else(i %in% act[ act$ID == j, "Activity"] , 1, 0)
    )

) %>%
  
  as.data.frame() %>%
  mutate(
    activities = rowSums( across( everything() ) ), # sum all activities per participant
    !!!setNames(rep( NA, length( unique(map$type) ) ), unique( map$type) ),
    !!!setNames(rep( NA, length( unique(map$category) ) ), unique( map$category) ),
    across( unique(map$type), ~ rowSums( across( all_of( map[map$type == cur_column(), "activity"] ) ) ) ),
    across( unique(map$category), ~ rowSums( across( all_of( map[map$category == cur_column(), "activity"] ) ) ) )
  ) %>%
  rownames_to_column("ID") %>%
    
  # join other data sets
  left_join(cog, by = "ID") %>%
  left_join(dem, by = "ID")


# FREQUENCIES ----

# prepare data
d2 <- act %>%
  
  left_join( dem[ , c("ID","Age_years")] ) %>%
  filter( complete.cases(Intensity) ) %>%
  left_join(cog) %>%
  mutate(
    
    # re-format binary variables
    Seasonal = factor(Seasonal),
    Activity_type = factor(Activity_type),
    Activity = factor(Activity),
    ID = factor(ID),
    
    # intensity/frequency (outcome) measures
    Intensity = as.integer( if_else(Intensity == "do not know", NA, Intensity) ),
    logIntensity = log( as.numeric(Intensity)-1 ),
    
    # time (predictor) variables
    Time_past = Time_bin - Age_years,
    Time_bin = factor(Time_bin, levels = unique(act$Time_bin), ordered = T),
    Time_num = as.numeric(Time_bin),
    #Time_log = log(Time_num),
    
    # PPC categories
    Type_Season_Time = paste(Activity_type, Seasonal, Time_bin, sep = "_"),
    SA_Type_Time = paste(SA, Activity_type, Time_bin, sep = "_"),
    SA_Season_Time = paste(SA, Seasonal, Time_bin, sep = "_")
    
  ) %>%
  within( . , {
    contrasts(SA) <- -contr.sum(2)/2 # nonSA = -0.5, SA = 0.5
    contrasts(Seasonal) <- -contr.sum(2)/2 # non-seasonal = -0.5, seasonal = 0.5
    contrasts(Activity_type) <- -contr.sum(2)/2 # mental = -0.5, physical = 0.5
  } )


## MODEL FITTING ----

# set-up prior
# (will be same for seasonal and non-seasonal models)
prior <-  c(
  prior(normal(1, 3), class = "Intercept"),
  prior(normal(0, 3), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma"),
  prior(lkj(2), class = "cor")
  #prior(dirichlet(1), class = "simo") # default Dirichlet(1) for monotonic terms in the ordered time model
)

### ---- non-seasonal activities (physical vs mental) ----

# set-up model formulas
form_n <- list(

  linear_time = bf( logIntensity ~ 1 + Time_num * SA * Activity_type + (1 + Time_num | ID) + (1 + Time_num | Category) ),
  #log_time = bf( logIntensity ~ 1 + Time_log * SA * Activity_type + (1 + Time_log | ID) + (1 + Time_log | Category) ),
  reverse_time = bf( logIntensity ~ 1 + Time_past * SA * Activity_type + (1 + Time_past | ID) + (1 + Time_past | Category) ),
  ordered_time = bf( logIntensity ~ 1 + mo(Time_bin) * SA * Activity_type + (1 + mo(Time_bin) | ID) + (1 + mo(Time_bin) | Category) )
  
)

# fit the models in a loop
fit_n <- lapply(
  
  setNames( names(form_n), names(form_n) ),
  function(i)
    
    brm(
      formula = form_n[[i]],
      family = gaussian(link = "identity"),
      prior = prior,
      data = subset(d2, Seasonal == F),
      file = here( "models", paste0(i,"_nonseasonal.rds") ),
      file_refit = "on_change"
    )
  
)


### ---- physical activities (seasonal vs non-seasonal) ----

# set-up model formulas
form_p <- list(
  
  linear_time = bf( logIntensity ~ 1 + Time_num * SA * Seasonal + (1 + Time_num | ID) + (1 + Time_num | Category) ),
  #log_time = bf( logIntensity ~ 1 + Time_log * SA * Seasonal + (1 + Time_log | ID) + (1 + Time_log | Category) ),
  reverse_time = bf( logIntensity ~ 1 + Time_past * SA * Seasonal + (1 + Time_past | ID) + (1 + Time_past | Category) ),
  ordered_time = bf( logIntensity ~ 1 + mo(Time_bin) * SA * Seasonal + (1 + mo(Time_bin) | ID) + (1 + mo(Time_bin) | Category) )
  
)

# fit the models in a loop
fit_p <- lapply(
  
  setNames( names(form_p), names(form_p) ),
  function(i)
    
    brm(
      formula = form_p[[i]],
      family = gaussian(link = "identity"),
      prior = prior,
      data = subset(d2, Activity_type == "physical"),
      file = here( "models", paste0(i,"_physical.rds") ),
      file_refit = "on_change"
    )
  
)


## CONVERGENCE SOFT CHECKS  ----

for ( i in names(fit_n) ) {
  
  print(i)
  plot(fit_n[[i]], ask = F)
  plot(fit_p[[i]], ask = F)
  
}

## MODEL COMPARISONS ----

# compare expected out-of-sample predictive performance approximated via PSIS-LOO
loo <- list(
  
  nonseasonal = with( fit_n, loo(linear_time, reverse_time, ordered_time) ),
  physical = with( fit_p, loo(linear_time, reverse_time, ordered_time) )
  
)

# save results of loo comparisons
write.table(
  
  x = lapply(
    
    names(loo),
    function(i)
      
      loo[[i]]$diffs %>%
      as.data.frame() %>%
      mutate(type = i, .before = 1) %>%
      rownames_to_column("model")
    
  ) %>%
    
    do.call( rbind.data.frame, . ),
  
  file = here("tables","loo_model_comparisons.csv"),
  row.names = F,
  quote = F,
  sep = ","

)

# Selecting ordered time in both cases and will continue with them only.

