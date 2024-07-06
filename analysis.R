# Use this script to address research questions regarding retrospectively remembered activities with respect to cognitive SA.

rm( list = ls() ) # clear environment

library(here)
library(openxlsx)
library(tidyverse)
library(effsize)

sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders


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


# FREQUENCIES ----

# prepare table for activity frequencies per SA status
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

# plot the frequencies per seasonal/non-seasonal activities
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
        mutate( category = i, .before = 1 )
    ) %>%
      
      do.call( rbind.data.frame, . ) %>%
      mutate(
        `Successful Aging: ` = if_else(SA == "SA", T, F),
        category = factor(category, levels = c( "activities", unique(map$type), unique(map$category) ), ordered = T),
        n = as.integer(n)
      ) %>%
      
      ggplot() +
      aes(x = Frequency, y = n, fill = `Successful Aging: `) +
      geom_bar( colour = "black", stat = "identity", width = .7, position = position_dodge(width = .7) ) +
      #geom_text( aes(label = n), position = position_dodge(width = .7), vjust = -.6 ) +
      #scale_y_continuous(name = NULL, labels = NULL) +
      scale_x_continuous( labels = seq(0,15,1), breaks = seq(0,15,1) ) +
      labs(y = NULL) +
      facet_wrap( ~ category, scales = "free", ncol = 3 ) +
      theme_bw( base_size = 14 ) +
      theme( axis.ticks.y = element_blank(), legend.position = "bottom" )
    
    # save it
    ggsave(
      plot = last_plot(),
      filename = here( "figs", paste0(x,"_activities_frequencies.jpg") ),
      dpi = 300,
      width = 12.6,
      height = 13.3
    )

  }
  
)


# INTENSITIES ----

# go with GLMMs for the intensity, can use two types: (i) log-responses, (ii) ordered-logit
