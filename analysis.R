# Use this script to address research questions regarding retrospectively remembered activities with respect to cognitive SA.

rm( list = ls() ) # clear environment

library(here)
library(openxlsx)
library(tidyverse)


# DATA ----

# read them
for ( i in names( readRDS("_data.rds") ) ) assign( i, readRDS("_data.rds")[[i]] )

# prepare a wide data set with all activites and their counts disregarding intensity and seasoness
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
    activites = rowSums( across( everything() ) ), # sum all activities per patients
    !!!setNames(rep( NA, length( unique(map$type) ) ), unique( map$type) ),
    !!!setNames(rep( NA, length( unique(map$category) ) ), unique( map$category) ),
    across( unique(map$type), ~ rowSums( across( all_of( map[map$type == cur_column(), "activity"] ) ) ) ),
    across( unique(map$category), ~ rowSums( across( all_of( map[map$category == cur_column(), "activity"] ) ) ) )
  ) %>%
  rownames_to_column("ID") %>%
  
  # join other data sets
  left_join(cog, by = "ID") %>%
  left_join(dem, by = "ID")
