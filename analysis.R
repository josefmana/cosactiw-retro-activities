# Use this script to address research questions regarding retrospectively remembered activities with respect to cognitive SA.

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

library(here)
library(openxlsx)
library(tidyverse)
library(effsize)
library(brms)
library(bayesplot)

color_scheme_set("viridisA")
theme_set( theme_bw(base_size = 12) )

sapply( c("figs","tabs","mods"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders


# UTILS ----

rprint <- function(x, d = 2) sprintf( paste0("%.",d,"f"), round(x, d) )
ciprint <- function(x, d = 2) paste0( "[", paste( rprint(x, d), collapse = ", "), "]" )
zerolead <- function(x, d = 3) ifelse( x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = T) )
msd <- function(x, d = 2) paste0( rprint( mean(x, na.rm = T), d ), " ± ", rprint( sd(x, na.rm = T), d ) )


# DATA ----

# read them
for ( i in names( readRDS("_data.rds") ) ) assign( i, readRDS("_data.rds")[[i]] )

# prepare a wide data set with all activites and their counts disregarding intensity and seasoness
d1 <- lapply(
  
  setNames( c("all","non-seasonal","seasonal"), c("all","non-seasonal","seasonal") ),
  function(x)
    
    sapply(
      
      map$activity,
      function(i)
        
        sapply(
          unique(act$ID),
          function(j) {
            if (x == "all") if_else(i %in% act[ act$ID == j, "Activity"] , 1, 0)
            else if(x == "seasonal") if_else(i %in% with( subset(act, Seasonal == T), Activity[ID == j] ), 1, 0)
            else if(x == "non-seasonal") if_else(i %in% with( subset(act, Seasonal == F), Activity[ID == j] ), 1, 0)
          }
        )
    ) %>%
    
    as.data.frame() %>%
    mutate(
      activities = rowSums( across( everything() ) ), # sum all activities per patients
      !!!setNames(rep( NA, length( unique(map$type) ) ), unique( map$type) ),
      !!!setNames(rep( NA, length( unique(map$category) ) ), unique( map$category) ),
      across( unique(map$type), ~ rowSums( across( all_of( map[map$type == cur_column(), "activity"] ) ) ) ),
      across( unique(map$category), ~ rowSums( across( all_of( map[map$category == cur_column(), "activity"] ) ) ) )
    ) %>%
    rownames_to_column("ID") %>%
    
    # join other data sets
    left_join(cog, by = "ID") %>%
    left_join(dem, by = "ID")

)

# check participants who reported the same activity as both seasonal and non-seasonal
cbind.data.frame(d1$all$ID, d1$all$activities, d1$seasonal$activities, d1$`non-seasonal`$activities) %>%
  `colnames<-`( c("ID","all","seas","nseas") ) %>%
  mutate(ctrl = seas + nseas) %>%
  mutate(ctrl0 = all - ctrl) %>%
  filter(ctrl0 != 0) %>%
  select(ID) %>%
  unlist(use.names = F)

# extract number of SA and nonSA individuals
n <- sapply( c("SA","nonSA"), function(i) nrow( subset(d1$all, SA == i) ) )


# FREQUENCIES ----

## TABLES ----

### table for activity frequencies per SA status ----
tab1 <- lapply(
  
  c( "activities", unique(map$type), unique(map$category) ),
  function(y)
    
    sapply(
      
      c("nonSA","SA"),
      function(i)
        
        sapply(
          names(d1),
          function(j)
            if_else( sum( subset(d1[[j]], SA == i)[ , y] ) == 0, "-", msd( subset(d1[[j]], SA == i)[ , y] ) )
        )
      
    ) %>%
    
    as.data.frame() %>%
    rownames_to_column("season") %>%
    mutate(category = y, .before = 1)
  
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  filter( !( ( category %in% c( "mental", unique( subset(map, type == "mental")$category ) ) ) & ( grepl("season",season) ) ) ) %>%
  mutate(
    type = sapply(
      1:nrow(.),
      function(i)
        case_when(
          category[i] == "activities" ~ "all",
          category[i] == "mental" | category[i] %in% unique( subset(map, type == "mental")$category ) ~ "mental",
          category[i] == "physical" | category[i] %in% unique( subset(map, type == "physical")$category ) ~ "physical"
        )
    ),
    category = case_when(
      category %in% c("activities","mental","physical") ~ "all",
      .default = category
    ),
    .before = 1
  ) %>%
  arrange(type)

# add statistical tests and effect sizes
tab1 <- sapply(
  
  1:nrow(tab1),
  function(i) {
    
    # extract variables to be used
    s <- tab1$season[i]
    v <- with( tab1, case_when( type[i] == "all" & category[i] == "all" ~ "activities", category[i] == "all" ~ type[i], .default = category[i]) )
    
    # calculate t-test
    t <- t.test( as.formula( paste0(v," ~ SA") ), data = d1[[s]] )
    d <- cohen.d( as.formula( paste0(v," ~ SA") ), data = d1[[s]] )
    
    # calculate Mann-Whiteny test
    wilcox <- wilcox.test( as.formula( paste0(v," ~ SA") ), data = d1[[s]] )
    vda <- VD.A( as.formula( paste0(v," ~ SA") ), data = d1[[s]] )
    
    # return it
    return( c(
      cohens_d = paste0( rprint(d$estimate), " ", ciprint(d$conf.int) ),
      t_stat = rprint(t$statistic),
      df = rprint(t$parameter),
      p_ttest = zerolead(t$p.value),
      VD.A = rprint(vda$estimate),
      W = rprint(wilcox$statistic),
      p_mannwhitney = zerolead(wilcox$p.value)
    ) )
  }
  
) %>%
  
  t() %>%
  as.data.frame() %>%
  cbind.data.frame( tab1, . ) %>%
  mutate( across( everything(), ~ case_when( grepl("NaN",.x) ~ "-", is.na(.x) ~ "-", .default = .x ) ) )

# save the table as .csv
write.table(x = tab1, file = here("tabs","activities_frequencies.csv"), sep = ";", row.names = F, quote = F)


### tables of frequencies of participants reporting at least one activity in a category ----

# for all of them, Chi-square based p > .05
lapply(
  
  names(d1),
  function(x)
    
    d1[[x]] %>%
    select( SA, all_of(map$category) ) %>%
    mutate( across( all_of(map$category), ~ if_else(.x > 0, T, F) ) ) %>%
    group_by(SA) %>%
    summarise_all( list(sum) ) %>%
    column_to_rownames("SA") %>%
    t() %>%
    as.data.frame() %>%
    mutate( Type = unlist(sapply( rownames(.), function(i) with( map, unique( type[category == i] ) ) ), use.names = F), .before = 1 ) %>%
    rownames_to_column("Category") %>%
    write.table(
      file = here( "tabs", paste0(x,"_activity_categories_reports.csv") ),
      sep = ",",
      row.names = F,
      quote = F
    )
  
)


## FIGURES ----

### plot individual-level frequencies per seasonal/non-seasonal categories ----
lapply(
  
  names(d1),
  function(x) {
    
    # plot it
    lapply(
      
      c( "activities", unique(map$type), unique(map$category) ),
      function(i)
        
        d1[[x]] %>%
        count( get(i), SA ) %>%
        rename( "Frequency" = "get(i)" ) %>%
        complete( Frequency, SA, fill = list(n = 0) ) %>%
        mutate( category = if_else(i == "activities", "all-in", i), .before = 1 )
    ) %>%
      
      do.call( rbind.data.frame, . ) %>%
      mutate(
        `Successful Aging: ` = if_else(SA == "SA", T, F),
        category = factor(category, levels = c( "all-in", unique(map$type), unique(map$category) ), ordered = T),
        n = as.integer(n)
      ) %>%
      
      ggplot() +
      aes(x = Frequency, y = n, fill = `Successful Aging: `) +
      geom_bar( colour = "black", stat = "identity", width = .7, position = position_dodge(width = .7) ) +
      scale_x_continuous( labels = seq(0,15,1), breaks = seq(0,15,1) ) +
      facet_wrap( ~ category, scales = "free", ncol = 3 ) +
      theme_bw( base_size = 14 ) +
      labs(
        y = NULL,
        title = paste0("Per-participant leisure activity categories report (", x, " activities)"),
        subtitle = paste0(
          "Bars show frequency (y-axis) of the number of activities (x-axis) reported by ", n["SA"], " SA and ", n["nonSA"], " non-SA participants."
        )
      ) +
      theme(
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5)
      )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figs", paste0(x,"_activity_categories_frequencies.jpg") ),
      dpi = 300,
      width = 12.6,
      height = 13.3
    )

  }
)


### plot sample-level frequencies per seasonal/non-seasonal categories ----
lapply(
  
  names(d1),
  function(x) {
    
    # plot it
    d1[[x]] %>%
      select( SA, all_of(map$category) ) %>%
      mutate( across( all_of(map$category), ~ if_else(.x > 0, T, F) ) ) %>%
      group_by(SA) %>%
      summarise_all( list(sum) ) %>%
      pivot_longer(cols = -SA, names_to = "Category", values_to = "Frequency") %>%
      mutate(
        Type = unlist(sapply( 1:nrow(.), function(i) with( map, unique( type[category == Category[i]] ) ) ), use.names = F),
        `Successful Aging: ` = if_else(SA == "SA", T, F)
      ) %>%
      ggplot() +
      aes(x = reorder(Category, Frequency), y = Frequency, fill = `Successful Aging: `) +
      geom_bar( stat = "identity", width = .7, position = position_dodge(width = .7), colour = "black" ) +
      geom_text( aes(label = Frequency), position = position_dodge(width = .7), hjust = -.33 ) +
      scale_y_continuous( limits = c(0,130) ) +
      facet_wrap( ~ Type, scales = "free_y" ) +
      coord_flip() +
      labs(
        x = NULL,
        title = paste0("Leisure activity categories (", x, " activities)"),
        subtitle = paste0(
          "Bars show number of participants (x-axis) out of ", n["SA"], " SA and ",
          n["nonSA"], " non-SA participants\nthat reported at least one activity from each respective category (y-axis)."
        )
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5)
      )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figs", paste0(x,"_activitiy_categories_reports.jpg") ),
      dpi = 300,
      width = 10,
      height = 10.6
    )

  }
)


### plot sample-level frequencies per seasonal/non-seasonal activities ----
lapply(
  
  names(d1),
  function(x) {
    
    # plot it
    d1[[x]] %>%
      select( SA, all_of(map$activity) ) %>%
      mutate( across( all_of(map$activity), ~ if_else(.x > 0, T, F) ) ) %>%
      group_by(SA) %>%
      summarise_all( list(sum) ) %>%
      pivot_longer(cols = -SA, names_to = "Activity", values_to = "Frequency") %>%
      mutate(
        Category = factor(
          unlist(sapply( 1:nrow(.), function(i) with( map, category[activity == Activity[i]] ) ), use.names = F),
          levels = c("cognitively_demanding","aerobic_exercise","community_activities","flexibility_health_exercise","creative_productive","house_and_garden","mixed_mental","mixed_exercise","receptive","strengthening_resistance_exercise"),
          ordered = T
        ),
        Type = unlist(sapply( 1:nrow(.), function(i) with( map, type[activity == Activity[i]] ) ), use.names = F),
        `Successful Aging: ` = if_else(SA == "SA", T, F)
      ) %>%
      
      ggplot() +
      aes(x = reorder(Activity, Frequency), y = Frequency, fill = `Successful Aging: `) +
      geom_bar( stat = "identity", width = .7, position = position_dodge(width = .7), colour = "black" ) +
      scale_y_continuous( limits = c(0,130) ) +
      facet_wrap( ~ Category * Type, scales = "free_y", ncol = 2 ) +
      coord_flip() +
      labs(
        x = NULL,
        title = paste0("Leisure activities per types and categories (", x, " activities)"),
        subtitle = paste0(
          "Bars show number of participants (x-axis) out of ", n["SA"], " SA and ",
          n["nonSA"], " non-SA participants\nthat reported at least one respective activity (y-axis)."
        )
      ) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5)
      )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figs", paste0("_",x,"_single_activities_reports.jpg") ),
      dpi = 300,
      width = 10,
      height = 30.6
    )
    
  }
)


## MODELS ----

# prepare data for a GLMM
d2 <- d1$all %>%
  select( ID, SA, all_of(map$category) ) %>%
  pivot_longer(
    cols = all_of(map$category),
    values_to = "Frequency",
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

# zero-inflated Poisson model
fit0poiss <- brm(
  formula = Frequency ~ 1 + (1 | ID) + (1 + SA || Type / Category),
  family = zero_inflated_poisson(link = "log", link_zi = "logit"),
  prior = NULL,
  data = d2,
  file = here("mods","frequencies_zero_inflated_poisson.rds")
)

# zero-inflated Negative-Binomial model
fit0negbin <- brm(
  formula = Frequency ~ 1 + (1 | ID) + (1 + SA || Type / Category),
  family = zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit"),
  prior = NULL,
  data = d2,
  file = here("mods","frequencies_zero_inflated_negbin.rds")
)


# INTENSITIES ----

# prepare data
d3 <- act %>%
  filter( complete.cases(Intensity) ) %>%
  mutate(
    Intensity = as.integer( if_else(Intensity == "do not know", NA, Intensity) ),
    logIntensity = log( as.numeric(Intensity)-1 ),
    Time_bin = factor(Time_bin, levels = unique(act$Time_bin), ordered = T)
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
  file = here("mods","explore_gauss.rds")
  
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


### ORDERED LOGIT ----

# fit a Gaussian model on log intensities and explore its posterior predictions w.r.t. data
fit1 <- brm(

  formula = Intensity ~ 1 + (1 | ID) + (1 | Activity_type / Category) + (1 | Time_bin),
  family = cumulative(link = "logit"),
  prior = NULL,
  data = d3,
  file = here("mods","explore_logit.rds")
  
)
