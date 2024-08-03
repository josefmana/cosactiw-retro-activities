# Use this script to describe retrospectively recalled activities with respect to cognitive SA

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

library(here)
library(tidyverse)
library(brms)
library(bayesplot)
library(bayestestR)
library(patchwork)

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
    SA =
      factor(
        SA,
        levels = c("nonSA","SA"),
        ordered = T
      ),
    Type =
      factor(
        unlist(
          sapply( 1:nrow(.), function(i) with( map, unique( type[category == Category[i]] ) ) ),
          use.names = F
        )
      ),
    .before = Category
  ) %>%
  within( . , {
    contrasts(SA) <- -contr.sum(2)/2 # nonSA = -0.5, SA = 0.5
    contrasts(Type) <- -contr.sum(2)/2 # mental = -0.5, physical = 0.5
  } )

# Poisson model
#fit0poiss <- brm(
#  formula = Count ~ 1 + SA * Type + (1 | ID) + (1 | Category),
#  family = poisson(link = "log"),
#  prior = c(
#    prior("normal(-2.3,2.5)", class = "Intercept"),
#    prior("normal(0,1)", class = "b"),
#    prior("exponential(1)", class = "sd")
#  ),
#  data = d2,
#  control = list(adapt_delta = .9),
#  file = here("models","count_poisson.rds"),
#  file_refit = "on_change"
#)


# FREQUENCIES ----

# prepare data
d3 <- act %>%
  
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
      data = subset(d3, Seasonal == F),
      file = here( "models", paste0(i,"_nonseasonal.rds") ),
      file_refit = "on_change"
    )
  
)

# posterior predictive check for conditional means
for ( i in names(fit_n) ) {
  
  print(
    pp_check(
      object = c(na.omit( subset(d3, Seasonal == F) )$logIntensity),
      yrep = posterior_predict( fit_n[[i]], newdata = na.omit( subset(d3, Seasonal == F) ) ),
      fun = ppc_stat_grouped,
      stat = "mean",
      group = na.omit( subset(d3, Seasonal == F) )$SA_Type_Time
    )
  )
  
}


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
      data = subset(d3, Activity_type == "physical"),
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
with( fit_n, loo(linear_time, reverse_time, ordered_time) )
with( fit_p, loo(linear_time, reverse_time, ordered_time) )

# Selecting ordered time in both cases and will continue with them only.
fit <- list(nonseasonal = fit_n$ordered_time, physical = fit_p$ordered_time)
rm(fit_n, fit_p, form_n, form_p, prior, i) # remove the old model objects
gc() # clean memory if possible

## POSTERIOR PREDICTIVE CHECKS ----

# density plots
lapply(
  
  names(fit),
  function(i) {
    
    # prepare parameters for plotting
    if (i == "nonseasonal") {
      
      data <- na.omit( subset(d3, Seasonal == F) )
      x <- "SA_Type_Time"
      txt <- "activity type"
      
      color_scheme_set("viridisA")
      
    } else if (i == "physical") {
      
      data <- na.omit( subset(d3, Activity_type == "physical") )
      x <- "SA_Season_Time"
      txt <- "seasoness"
      
      color_scheme_set("blue")
      
    }
    
    ### ---- density overlay ----
    
    # plot it
    ppc_plot(
      fit = fit[[i]],
      data = data,
      y = "logIntensity",
      x = x,
      labs = list(
        "Observed (thick lines) and predicted (thin lines) distributions of log(Intensity) of recalled leisure activities",
        paste0("Cells represent different combinations of SA status, ",txt,", and time bin.")
      ),
      meth = "dens_overlay_grouped",
      draws = sample(1:4e3, 1e2),
      stat = NULL
    )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("ppc_density_",i,"_ordered_time.jpg") ),
      dpi = 300,
      width = 12.6,
      height = 13.3
    )
    
    ### ---- stat plots ----
    
    lapply(
      
      c("mean","sd"),
      function(j) {
        
        # plot it
        ppc_plot(
          fit = fit[[i]],
          data = data,
          y = "logIntensity",
          x = x,
          labs = list(
            paste0("Observed (thick bars) and predicted (histograms) ",j," of log(Intensity) of recalled leisure activities"),
            paste0("Cells represent different combinations of SA status, ",txt,", and time bin.")
          ),
          meth = "stat_grouped",
          draws = 1:4e3,
          stat = j
        )
        
        # save it
        ggsave(
          plot = last_plot(),
          filename = here( "figures", paste0("ppc_stat_",j,"_",i,"_ordered_time.jpg") ),
          dpi = 300,
          width = 12.6,
          height = 13.3
        )
        
      }
    )
    
  }
)


## POST-PROCESSING ----

# prepare a data frame for posterior inference
# will be used for both nonseasonal and physical only models
d_seq <- with(

  d3,
  expand.grid(
    ID = NA,
    Category = NA,
    Activity_type = c( levels(Activity_type), NA ),
    Seasonal = c( levels(Seasonal), NA ),
    SA = c( levels(SA), NA ),
    Time_bin = levels(Time_bin)

  )
) %>%
  
  mutate( Time_bin = as.ordered(Time_bin) ) %>%
  within( . , {
    contrasts(SA) <- -contr.sum(2)/2
    contrasts(Seasonal) <- -contr.sum(2)/2
    contrasts(Activity_type) <- -contr.sum(2)/2
  } )

# extract posterior draws for each combination of variables in each model
ppred <- lapply(
  
  setNames( names(fit), names(fit) ),
  function(i) {
    
    # prepare new data and column names for posterior E(X)s
    if (i == "nonseasonal") {
      
      ndat <- subset(d_seq, Seasonal == F)
      cnms <- with(ndat, paste(SA, Activity_type, Time_bin, sep = "_") )
      
    } else if (i == "physical") {
      
      ndat <- subset(d_seq, Activity_type == "physical")
      cnms <- with(ndat, paste(SA, Seasonal, Time_bin, sep = "_") )
      
    }
    
    # compute it
    # re_formula NA to ignore participant-level terms
    return( posterior_epred(object = fit[[i]], newdata = ndat, re_formula = NA) %>% `colnames<-`(cnms) )

  }
)

# add median and 95% PPIs for each prediction to the d_seq
post_sum <- lapply(
  
  setNames( names(ppred), names(ppred) ),
  function(i) {
    
    # prepare new data and column names for posterior E(X)s
    if (i == "nonseasonal") d <- subset(d_seq, Seasonal == F) else if (i == "physical") d <- subset(d_seq, Activity_type == "physical")
    
    # add summaries
    cbind(

      d,
      apply(
        ppred[[i]],
        2,
        function(x)
          c( logxEstimate = median(x),
             logxETI_low = quantile(x, prob = .025, names = F),
             logxETI_high = quantile(x, prob = .975, names = F),
             rawxEstimate = median( exp(x) ),
             rawxETI_low = quantile( exp(x), prob = .025, names = F),
             rawxETI_high = quantile( exp(x), prob = .975, names = F)
            )
      ) %>%
        t()

    ) %>%
      
      as.data.frame() %>%
      pivot_longer(
        cols = starts_with("log") | starts_with("raw"),
        names_to = c("scale","term"),
        names_pattern = "(.*)x(.*)",
        values_to = "value"
      ) %>%
      pivot_wider(
        names_from = "term",
        values_from = "value",
        id_cols = c("Activity_type","Seasonal","SA","Time_bin","scale")
      )
  }
)


### ---- interaction plots ----

#### ---- nonseasonal model ----

# prepare
fig_nonseasonal <- lapply(
  
  setNames( c("log","raw"), c("log","raw") ),
  function(i) {
    
    # prepare data
    data <-
      post_sum$nonseasonal %>%
      filter( !is.na(Activity_type) ) %>%
      filter( !is.na(SA) ) %>%
      filter( scale == i )
    
    # plot both ways
    list(
      
      data %>%
        ggplot() +
        aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high, colour = Activity_type, group = Activity_type) +
        geom_point( position = position_dodge(width = .5), size = 4 ) +
        geom_linerange( position = position_dodge(width = .5), size = 1.5 ) +
        geom_line(linetype = "dashed", linewidth = .8) +
        scale_colour_manual( values = c("#56B4E9","#E69F00") ) +
        facet_wrap( ~ SA, nrow = 2 ) +
        theme(panel.grid = element_blank(), legend.position = "bottom"),
      
      data %>%
        ggplot() +
        aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high, colour = SA, group = SA) +
        geom_point( position = position_dodge(width = .5), size = 4 ) +
        geom_linerange( position = position_dodge(width = .5), size = 1.5 ) +
        geom_line(linetype = "dashed", linewidth = .8) +
        scale_colour_manual( values = c("#CC79A7","#999999") ) +
        facet_wrap( ~ Activity_type, nrow = 2 ) +
        theme(panel.grid = element_blank(), legend.position = "bottom")
      
    )
    
  }
)

# plot and save
lapply(
  
  names(fig_nonseasonal),
  function(i) {
    
    # plot it
    with( fig_nonseasonal, get(i)[[1]] | get(i)[[2]] ) +
      plot_layout(axis_titles = "collect") +
      plot_annotation(
        title = "Three-way interaction between time bin, SA, and activity type",
        subtitle = paste0("Left and right sets of panels are representation of the same interaction on a ",i," scale"),
        theme = theme( plot.title = element_text(hjust = .5, face = "bold"), plot.subtitle = element_text(hjust = .5) )
      )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("conditional_means_nonseasonal_",i,"_scale.jpg") ),
      dpi = 300,
      width = 12.6,
      height = 8.31
    )
    
  }
)

#### ---- physical model ----

# prepare it
fig_physical <- lapply(
  
  setNames( c("log","raw"), c("log","raw") ),
  function(i) {
    
    # prepare data
    data <-
      post_sum$physical %>%
      filter( !is.na(Seasonal) ) %>%
      filter( !is.na(SA) ) %>%
      filter( scale == i )
    
    # plot both ways
    list(
      
      data %>%
        ggplot() +
        aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high, colour = Seasonal, group = Seasonal) +
        geom_point( position = position_dodge(width = .5), size = 4 ) +
        geom_linerange( position = position_dodge(width = .5), size = 1.5 ) +
        geom_line(linetype = "dashed", linewidth = .8) +
        scale_colour_manual( values = c("#F0E442","#009E73") ) +
        facet_wrap( ~ SA, nrow = 2 ) +
        theme(panel.grid = element_blank(), legend.position = "bottom"),
      
      data %>%
        ggplot() +
        aes(x = Time_bin, y = Estimate, ymin = ETI_low, ymax = ETI_high, colour = SA, group = SA) +
        geom_point( position = position_dodge(width = .5), size = 4 ) +
        geom_linerange( position = position_dodge(width = .5), size = 1.5 ) +
        geom_line(linetype = "dashed", linewidth = .8) +
        scale_colour_manual( values = c("#CC79A7","#999999") ) +
        facet_wrap( ~ Seasonal, nrow = 2 ) +
        theme(panel.grid = element_blank(), legend.position = "bottom")
      
    )
    
  }
)

# plot and save it
lapply(
  
  names(fig_physical),
  function(i) {
    
    # plot it
    with( fig_physical, get(i)[[1]] | get(i)[[2]] ) +
      plot_layout(axis_titles = "collect") +
      plot_annotation(
        title = "Three-way interaction between time bin, SA, and seasoness",
        subtitle = paste0("Left and right sets of panels are representation of the same interaction on a ",i," scale"),
        theme = theme( plot.title = element_text(hjust = .5, face = "bold"), plot.subtitle = element_text(hjust = .5) )
      )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figures", paste0("conditional_means_physical_",i,"_scale.jpg") ),
      dpi = 300,
      width = 12.6,
      height = 8.31
    )
    
  }
)

