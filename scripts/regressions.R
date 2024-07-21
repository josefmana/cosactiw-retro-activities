# Use this script to describe retrospectively remembered activities with respect to cognitive SA

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

library(here)
library(tidyverse)
library(brms)
library(bayesplot)

color_scheme_set("viridisA")
theme_set( theme_bw(base_size = 12) )

sapply( c("figures","tables","models"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders


# UTILS ----

rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )

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


## COUNTS ----

# prepare data for a GLMM
d2 <- d1 %>%
  select( ID, SA, all_of(map$category) ) %>%
  pivot_longer(
    cols = all_of(map$category),
    values_to = "Count",
    names_to = "Category"
  ) %>%
  mutate(
    Type =
      unlist(
        sapply( 1:nrow(.), function(i) with( map, unique( type[category == Category[i]] ) ) ),
        use.names = F
      ),
    .before = Category
  )

# Poisson model
fit0poiss <- brm(
  formula = Count ~ 1 + (1 | ID) + (1 + SA || Type / Category),
  family = poisson(link = "log"),
  prior = NULL,
  data = d2,
  control = list(adapt_delta = .9),
  file = here("models","count_poisson.rds"),
  file_refit = "on_change"
)

# Negative-Binomial model
fit0negbin <- brm(
  formula = Count ~ 1 + (1 | ID) + (1 + SA || Type / Category),
  family = negbinomial(link = "log", link_shape = "log"),
  prior = NULL,
  data = d2,
  control = list(adapt_delta = .9),
  file = here("models","count_negbinomial.rds"),
  file_refit = "on_change"
)


# FREQUENCIES ----

# prepare data
d3 <- act %>%
  filter( complete.cases(Intensity) ) %>%
  mutate(
    Intensity = as.integer( if_else(Intensity == "do not know", NA, Intensity) ),
    logIntensity = log( as.numeric(Intensity)-1 ),
    Time_bin = factor(Time_bin, levels = unique(act$Time_bin), ordered = T),
    Time_num = as.numeric(Time_bin),
    Time_Type = paste0(Activity_type,"_",Time_bin),
    Type_Season_Time = paste(Activity_type, Seasonal, Time_bin, sep = "_"),
    Time_log = log(Time_num)
  ) %>%
  left_join(cog)


## EXPLORATION ----

### GAUSSIAN ----

# fit a Gaussian model on log intensities and explore its posterior predictions w.r.t. data
fit0 <- brm(
  
  formula = logIntensity ~ 1 + (1 | ID) + (1 | Activity_type / Category) + (1 | Time_bin),
  family = gaussian(link = "identity"),
  prior = NULL,
  data = d3,
  file = here("models","explore_gauss.rds")
  
)

# convergence soft check
plot(fit0)

# some posterior checks per activity category
brms::pp_check(fit0, type = "dens_overlay_grouped", group = "Activity_type:Category", ndraws = 100) # missing mark
brms::pp_check(fit0, type = "stat_grouped", group = "Activity_type:Category", stat = "mean") # good
brms::pp_check(fit0, type = "stat_grouped", group = "Activity_type:Category", stat = "median") # missed most of them
brms::pp_check(fit0, type = "stat_grouped", group = "Activity_type:Category", stat = "sd") # missed most of them

# the model misses mean of logIntensity for seasonal activities in a variable-specific way
pp_check(
  object = c( na.omit(d3$logIntensity) ),
  yrep = posterior_predict( fit0, newdata = na.omit(d3) ),
  fun = ppc_stat_grouped,
  stat = "mean",
  group = na.omit(d3)$Seasonal
)

# on the other hand, the model get mean of logIntensity for SA vs non-SA subjects similarly well
pp_check(
  object = c( na.omit(d3$logIntensity) ),
  yrep = posterior_predict( fit0, newdata = na.omit(d3) ),
  fun = ppc_stat_grouped,
  stat = "mean",
  group = na.omit(d3)$SA
)



fit2 <- brm(
  
  formula = logIntensity ~ 1 + Time_num * (Activity_type + Seasonal) + (1 + Time_num | ID) + (1 + Time_num | Category),
  family = gaussian(link = "identity"),
  prior = c(
    prior(normal(1, 3), class = "Intercept"),
    prior(normal(0, 3), class = "b"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma"),
    prior(lkj(2), class = "cor")
  ),
  data = d3,
  file = here("models","explore_gauss2.rds"),
  file_refit = "on_change"
  
)



fit3 <- brm(
  
  formula = logIntensity ~ 1 + Time_num * (SA + Seasonal) + (1 + Time_num | ID) + (1 + Time_num | Activity_type / Category),
  family = gaussian(link = "identity"),
  prior = c(
    prior(normal(1, 3), class = "Intercept"),
    prior(normal(0, 3), class = "b"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma"),
    prior(lkj(2), class = "cor")
  ),
  data = d3,
  file = here("models","explore_gauss3.rds"),
  file_refit = "on_change"
  
)


### ORDERED LOGIT ----

# fit a Gaussian model on log intensities and explore its posterior predictions w.r.t. data
fit2 <- brm(
  
  formula = Intensity ~ 1 + (1 | ID) + (1 | Activity_type / Category) + (1 | Time_bin),
  family = cumulative(link = "logit"),
  prior = NULL,
  data = d3,
  file = here("models","explore_logit.rds")
  
)

