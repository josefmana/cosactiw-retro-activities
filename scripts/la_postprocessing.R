# Use this script to describe selected models of retrospectively recalled activities with respect to cognitive SA.

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

sapply( c("figures","tables"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders


# DATA ----

# read them
for ( i in names( readRDS("_data.rds") ) ) assign( i, readRDS("_data.rds")[[i]] )

# prepare data frame for analyses
df <- act %>%
  
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

# models
fit <- list(
  
  nonseasonal = readRDS( here("models","ordered_time_nonseasonal.rds") ),
  physical = readRDS( here("models","ordered_time_physical.rds") )
  
)


# POSTERIOR PREDICTIVE CHECKS ----

# density plots
lapply(
  
  names(fit),
  function(i) {
    
    # prepare parameters for plotting
    if (i == "nonseasonal") {
      
      data <- na.omit( subset(df, Seasonal == F) )
      x <- "SA_Type_Time"
      txt <- "activity type"
      
      color_scheme_set("viridisA")
      
    } else if (i == "physical") {
      
      data <- na.omit( subset(df, Activity_type == "physical") )
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
      filename = here( "figures", paste0("_ppc_density_",i,"_ordered_time.jpg") ),
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
          filename = here( "figures", paste0("_ppc_stat_",j,"_",i,"_ordered_time.jpg") ),
          dpi = 300,
          width = 12.6,
          height = 13.3
        )
        
      }
    )
    
  }
)


# STATISTICAL INFERENCE ----

# prepare a data frame for posterior inference
# will be used for both nonseasonal and physical only models
d_seq <- with(
  
  df,
  expand.grid(
    ID = NA,
    Category = NA,
    Activity_type = c(levels(Activity_type), NA),
    Seasonal = c(levels(Seasonal), NA),
    SA = c(levels(SA), NA),
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
             explogxEstimate = median( exp(x) ),
             explogxETI_low = quantile( exp(x), prob = .025, names = F),
             explogxETI_high = quantile( exp(x), prob = .975, names = F)
          )
      ) %>%
        t()
      
    ) %>%
      
      as.data.frame() %>%
      pivot_longer(
        cols = starts_with("log") | starts_with("explog"),
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

# save the results
lapply(
  
  names(post_sum),
  function(i)
    
    write.table(
      x = post_sum[[i]],
      file = here( "tables", paste0("_la_",i,"_conditional_means.csv") ),
      sep = ",",
      row.names = F,
      quote = F
    )
  
)


## ---- interaction plots ----

### ---- nonseasonal model ----

# prepare
fig_nonseasonal <- lapply(
  
  setNames( c("log","explog"), c("log","back-transformed raw") ),
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
    
    # prepare filename
    fn <- paste0("_conditional_means_nonseasonal_", gsub("-| ","_", i),"_scale.jpg")
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here("figures",fn),
      dpi = 300,
      width = 12.6,
      height = 8.31
    )
    
  }
)

### ---- physical model ----

# prepare it
fig_physical <- lapply(
  
  setNames( c("log","explog"), c("log","back-transformed raw") ),
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
    
    # prepare filename
    fn <- paste0("_conditional_means_physical_", gsub("-| ","_", i),"_scale.jpg")
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here("figures",fn),
      dpi = 300,
      width = 12.6,
      height = 8.31
    )
    
  }
)


## ---- contrasts/pairwise comparisons ----

# prepare a summary data frame
post_comp <- list(
  
  nonseasonal = list(
    
    # activity type and SA comparisons
    expand_grid(
      diff = c("physical - mental", "SA - nonSA", "physical_SA - mental_SA", "physical_nonSA - mental_nonSA"),
      diff2 = NA,
      timebin = levels(d_seq$Time_bin),
      actype = NA,
      SA = NA
    ),
    
    # two-way interactions of interest
    expand_grid(
      diff = "physical - mental",
      diff2 = "SA - nonSA",
      timebin = levels(d_seq$Time_bin),
      actype = NA,
      SA = NA
    ),
    
    # time comparisons
    expand_grid(
      diff = data.frame( v2 = seq(30,80,5) ) %>%
        mutate(v1 = v2 + 5) %>%
        mutate(diff = paste0(v1," - ",v2) ) %>%
        select(diff) %>%
        unlist(use.names = F) %>%
        c("85 - 30"),
      diff2 = NA,
      timebin = NA,
      actype = c("mental","physical",NA),
      SA = c("SA","nonSA",NA)
    ),
    
    # two-way interactions with time
    expand_grid(diff = "85 - 30", diff2 = "SA - nonSA", timebin = NA, actype = c("mental","physical",NA), SA = NA),
    expand_grid(diff = "85 - 30", diff2 = "physical - mental", timebin = NA, actype = NA, SA = c("SA","nonSA",NA) )
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    mutate( # re-code
      SA = case_when(
        diff == "physical_SA - mental_SA" ~ "SA",
        diff == "physical_nonSA - mental_nonSA" ~ "nonSA",
        .default = SA
      ),
      diff = gsub("_SA|_nonSA", "", diff)
    ),
  
  physical = list(
    
    # activity type and SA comparisons
    expand_grid(
      diff = c("FALSE - TRUE", "SA - nonSA", "FALSE_SA - TRUE_SA", "FALSE_nonSA - TRUE_nonSA"),
      diff2 = NA,
      timebin = levels(d_seq$Time_bin),
      season = NA,
      SA = NA
    ),
    
    # two-way interactions of interest
    expand_grid(
      diff = "FALSE - TRUE",
      diff2 = "SA - nonSA",
      timebin = levels(d_seq$Time_bin),
      season = NA,
      SA = NA
    ),
    
    # time comparisons
    expand_grid(
      diff = data.frame( v2 = seq(30,80,5) ) %>%
        mutate(v1 = v2 + 5) %>%
        mutate(diff = paste0(v1," - ",v2) ) %>%
        select(diff) %>%
        unlist(use.names = F) %>%
        c("85 - 30"),
      diff2 = NA,
      timebin = NA,
      season = c("FALSE","TRUE",NA),
      SA = c("SA","nonSA",NA)
    ),
    
    # two-way interactions with time
    expand_grid(diff = "85 - 30", diff2 = "SA - nonSA", timebin = NA, season = c("FALSE","TRUE",NA), SA = NA),
    expand_grid(diff = "85 - 30", diff2 = "FALSE - TRUE", timebin = NA, season = NA, SA = c("SA","nonSA",NA) )
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    mutate( # re-code
      SA = case_when(
        diff == "FALSE_SA - TRUE_SA" ~ "SA",
        diff == "FALSE_nonSA - TRUE_nonSA" ~ "nonSA",
        .default = SA
      ),
      diff = gsub("_SA|_nonSA", "", diff)
    )
  
)

### ---- save results of the nonseasonal model ----

write.table(
  
  x = lapply(
    
    setNames( c("log","logexp"), c("log","logexp") ),
    function(s)
      
      with(
        
        post_comp$nonseasonal,
        # extracting
        sapply(
          
          1:nrow(post_comp$nonseasonal),
          function(i) {
            
            # extract terms
            term1 <- strsplit(diff[i], " - ")[[1]][1]
            term2 <- strsplit(diff[i], " - ")[[1]][2]
            
            # find out which variable is being compared
            # and extract column name appropriately
            if ( term1 %in% na.omit( unique(actype) ) ) {
              
              col1 <- paste(SA[i], term1, timebin[i], sep = "_")
              col2 <- paste(SA[i], term2, timebin[i], sep = "_")
              
            } else if ( term1 %in% na.omit( unique(SA) ) ) {
              
              col1 <- paste(term1, actype[i], timebin[i], sep = "_")
              col2 <- paste(term2, actype[i], timebin[i], sep = "_")
              
            } else if ( term1 %in% na.omit( unique(timebin) ) ) {
              
              col1 <- paste(SA[i], actype[i], term1, sep = "_")
              col2 <- paste(SA[i], actype[i], term2, sep = "_")
              
            }
            
            # extract the resulting posterior
            if (s == "logexp") dif <- exp( ppred$nonseasonal[ , col1] ) - exp( ppred$nonseasonal[ , col2] )
            else if (s == "log") dif <- ppred$nonseasonal[ , col1] - ppred$nonseasonal[ , col2]
            
            # compute interaction if relevant
            if( !is.na(diff2[i]) ) {
              
              # extract second terms
              term21 <- strsplit(diff2[i], " - ")[[1]][1]
              term22 <- strsplit(diff2[i], " - ")[[1]][2]
              
              # can be only SA or activity_type
              # find out columns to be compared
              if ( term1 %in% na.omit( unique(actype) ) ) {
                
                # term2 is SA
                col11 <- paste(term21, term1, timebin[i], sep = "_")
                col12 <- paste(term21, term2, timebin[i], sep = "_")
                col21 <- paste(term22, term1, timebin[i], sep = "_")
                col22 <- paste(term22, term2, timebin[i], sep = "_")
                
              }  else if ( term1 %in% na.omit( unique(timebin) ) ) {
                
                # term2 can be SA or activity_type
                if ( term21 %in% na.omit( unique(actype) ) ) {
                  
                  # term2 is activity_type
                  col11 <- paste(SA[i], term21, term1, sep = "_")
                  col12 <- paste(SA[i], term21, term2, sep = "_")
                  col21 <- paste(SA[i], term22, term1, sep = "_")
                  col22 <- paste(SA[i], term22, term2, sep = "_")
                  
                } else if ( term21 %in% na.omit( unique(SA) ) ) {
                  
                  # term2 is SA
                  col11 <- paste(term21, actype[i], term1, sep = "_")
                  col12 <- paste(term21, actype[i], term2, sep = "_")
                  col21 <- paste(term22, actype[i], term1, sep = "_")
                  col22 <- paste(term22, actype[i], term2, sep = "_")
                  
                }
                
              }
              
              # extract the resulting posterior
              if (s == "logexp") dif <- ( exp( ppred$nonseasonal[ , col11] ) - exp( ppred$nonseasonal[ , col12] ) ) - ( exp( ppred$nonseasonal[ , col21] ) - exp( ppred$nonseasonal[ , col22] ) )
              else if (s == "log") dif <- ( ppred$nonseasonal[ , col11] - ppred$nonseasonal[ , col12] ) - ( ppred$nonseasonal[ , col21] - ppred$nonseasonal[ , col22] )
              
            }
            
            # extract final results
            est <- paste0( rprint(median(dif), 2)," ", ciprint( unlist( eti(dif) )[-1] ) )
            pd <- paste0( rprint( 100 * c( p_direction(dif) )$pd, 2 ), "%" )
            return( c(Estimate = est, pd = pd) )
            
          }
        )
        
      ) %>%
      
      t() %>%
      cbind.data.frame( post_comp$nonseasonal, . ) %>%
      mutate(scale = s, .before = Estimate)
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    pivot_wider( names_from = scale, values_from = c(Estimate, pd) ) %>%
    relocate(pd_log, .after = Estimate_log),
  
  file = here("tables","_la_nonseasonal_pairwise_comparisons.csv"),
  row.names = F,
  quote = F,
  sep = ";"

)


### ---- save results of the physical model ----

write.table(
  
  x = lapply(
    
    setNames( c("log","logexp"), c("log","logexp") ),
    function(s)
      
      with(
        
        post_comp$physical,
        # extracting
        sapply(
          
          1:nrow(post_comp$physical),
          function(i) {
            
            # extract terms
            term1 <- strsplit(diff[i], " - ")[[1]][1]
            term2 <- strsplit(diff[i], " - ")[[1]][2]
            
            # find out which variable is being compared
            # and extract column name appropriately
            if ( term1 %in% na.omit( unique(season) ) ) {
              
              col1 <- paste(SA[i], term1, timebin[i], sep = "_")
              col2 <- paste(SA[i], term2, timebin[i], sep = "_")
              
            } else if ( term1 %in% na.omit( unique(SA) ) ) {
              
              col1 <- paste(term1, season[i], timebin[i], sep = "_")
              col2 <- paste(term2, season[i], timebin[i], sep = "_")
              
            } else if ( term1 %in% na.omit( unique(timebin) ) ) {
              
              col1 <- paste(SA[i], season[i], term1, sep = "_")
              col2 <- paste(SA[i], season[i], term2, sep = "_")
              
            }
            
            # extract the resulting posterior
            if (s == "logexp") dif <- exp( ppred$physical[ , col1] ) - exp( ppred$physical[ , col2] )
            else if (s == "log") dif <- ppred$physical[ , col1] - ppred$physical[ , col2]
            
            # compute interaction if relevant
            if( !is.na(diff2[i]) ) {
              
              # extract second terms
              term21 <- strsplit(diff2[i], " - ")[[1]][1]
              term22 <- strsplit(diff2[i], " - ")[[1]][2]
              
              # can be only SA or season
              # find out columns to be compared
              if ( term1 %in% na.omit( unique(season) ) ) {
                
                # term2 is SA
                col11 <- paste(term21, term1, timebin[i], sep = "_")
                col12 <- paste(term21, term2, timebin[i], sep = "_")
                col21 <- paste(term22, term1, timebin[i], sep = "_")
                col22 <- paste(term22, term2, timebin[i], sep = "_")
                
              }  else if ( term1 %in% na.omit( unique(timebin) ) ) {
                
                # term2 can be SA or season
                if ( term21 %in% na.omit( unique(season) ) ) {
                  
                  # term2 is season
                  col11 <- paste(SA[i], term21, term1, sep = "_")
                  col12 <- paste(SA[i], term21, term2, sep = "_")
                  col21 <- paste(SA[i], term22, term1, sep = "_")
                  col22 <- paste(SA[i], term22, term2, sep = "_")
                  
                } else if ( term21 %in% na.omit( unique(SA) ) ) {
                  
                  # term2 is SA
                  col11 <- paste(term21, season[i], term1, sep = "_")
                  col12 <- paste(term21, season[i], term2, sep = "_")
                  col21 <- paste(term22, season[i], term1, sep = "_")
                  col22 <- paste(term22, season[i], term2, sep = "_")
                  
                }
                
              }
              
              # extract the resulting posterior
              if (s == "logexp") dif <- ( exp( ppred$physical[ , col11] ) - exp( ppred$physical[ , col12] ) ) - ( exp( ppred$physical[ , col21] ) - exp( ppred$physical[ , col22] ) )
              else if (s == "log") dif <- ( ppred$physical[ , col11] - ppred$physical[ , col12] ) - ( ppred$physical[ , col21] - ppred$physical[ , col22] )
              
            }
            
            # extract final results
            est <- paste0( rprint(median(dif), 2)," ", ciprint( unlist( eti(dif) )[-1] ) )
            pd <- paste0( rprint( 100 * c( p_direction(dif) )$pd, 2 ), "%" )
            return( c(Estimate = est, pd = pd) )
            
          }
        )
        
      ) %>%
      
      t() %>%
      cbind.data.frame( post_comp$physical, . ) %>%
      mutate(scale = s, .before = Estimate)
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    pivot_wider( names_from = scale, values_from = c(Estimate, pd) ) %>%
    relocate(pd_log, .after = Estimate_log),
  
  file = here("tables","_la_physical_pairwise_comparisons.csv"),
  row.names = F,
  quote = F,
  sep = ";"
  
)
