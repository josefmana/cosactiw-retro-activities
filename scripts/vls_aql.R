# Use this script to analyse cognitive SA effect on social activities

rm( list = ls() ) # clear environment

library(here)
library(tidyverse)
library(psych)
library(lavaan)
library(effsize)

sapply( c("figures","tables","models"), function(i) if( !dir.exists(i) ) dir.create(i) ) # prepare folders
source( here("scripts","utils.R") ) # in-house functions


# DATA ----

for ( i in names( readRDS("_data.rds") ) ) assign( i, readRDS("_data.rds")[[i]] ) # read them

# list item-scale mapping
scls <- list(
  
  Total = unique(vls$item),
  Private = unique(vls$item)[1:6],
  Public = unique(vls$item)[7:11]
  
)

# wide data set for dimensionality, internal consistency and sum score differences
d1 <-
  vls %>%
  pivot_wider(names_from = item, values_from = response) %>%
  mutate(
 
    # SA status
    SA = factor( case_when(ID %in% subset(cog, SA == "SA")$ID ~ 1, ID %in% subset(cog, SA == "nonSA")$ID ~ 0) ),
    
    # demographics
    Age_years = unlist(sapply( 1:nrow(.), function(i) dem[dem$ID == ID[i], "Age_years"] ), use.names = F),
    Edu_years = unlist(sapply( 1:nrow(.), function(i) dem[dem$ID == ID[i], "Education_years"] ), use.names = F),
    
    # VLS-ALQ scores
    VLS_total = rowSums( across( all_of(scls$Total) ) ),
    VLS_private = rowSums( across( all_of(scls$Private) ) ),
    VLS_public = rowSums( across( all_of(scls$Public) ) ),
    
    # add before item scores
    .after = ID

  )

# long data set for IRT
d2 <-
  vls %>%
  mutate(
    SA = case_when(ID %in% subset(cog, SA == "SA")$ID ~ 1, ID %in% subset(cog, SA == "nonSA")$ID ~ 0),
    Age_years = unlist(sapply( 1:nrow(.), function(i) dem[dem$ID == ID[i], "Age_years"] ), use.names = F),
    Edu_years = unlist(sapply( 1:nrow(.), function(i) dem[dem$ID == ID[i], "Education_years"] ), use.names = F),
    .after = ID
  )


# RELIABILITY ----

## ---- CONFIRMATORY FACTOR ANALYSES ----

# set-up models for lavaan
mods <- list(
  
  Total = paste0( "VLS_total =~ ", paste(scls$Total, collapse = " + ") ),
  Multidimensional = "
  VLS_private =~ VLS_I_go_out_with_friends + VLS_I_visit_friends_or_relatives + VLS_I_attend_parties_or_celebrations + VLS_I_talk_to_a_person_on_the_phone_skype + VLS_I_invite_friends_over_to_my_house_for_dinner + VLS_I_go_out_to_eat_at_a_restaurant
  VLS_public =~ VLS_I_get_involved_in_political_activities + VLS_I_give_a_public_speech + VLS_I_attend_a_meeting_of_the_club + VLS_I_participate_in_organised_public_events + VLS_I_volunteer
  VLS_religious =~ VLS_I_attend_church_services
  ",
  `Second order` = "
  VLS_private =~ VLS_I_go_out_with_friends + VLS_I_visit_friends_or_relatives + VLS_I_attend_parties_or_celebrations + VLS_I_talk_to_a_person_on_the_phone_skype + VLS_I_invite_friends_over_to_my_house_for_dinner + VLS_I_go_out_to_eat_at_a_restaurant
  VLS_public =~ VLS_I_get_involved_in_political_activities + VLS_I_give_a_public_speech + VLS_I_attend_a_meeting_of_the_club + VLS_I_participate_in_organised_public_events + VLS_I_volunteer
  VLS_religious =~ VLS_I_attend_church_services
  VLS_total =~ VLS_private + VLS_public + VLS_religious
  ",
  Private = paste0( "VLS_private =~ ", paste(scls$Private, collapse = " + ") ),
  Public = paste0( "VLS_public =~ ", paste(scls$Public, collapse = " + ") )
  
)

# compute the CFAs
CFAs <- lapply(
  
  set_names( names(mods) ),
  function(i)
    cfa( mods[[i]], data = d1[ , scls$Total], estimator = "MLR", missing = "ML" )
  
)

# check convergence
sapply( names(CFAs), function(i) summary(CFAs[[i]])$header$optim.converged )

# extract CFA metrics of interest
cfind <- sapply(
  
  names(CFAs),
  function(i) with(
    
    summary(CFAs[[i]], fit.measures = T, standardized = T),
    c( n = data$nobs, fit[c( "npar", paste0( c("cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper"), ".robust") ) ] )
    
  )
  
) %>%
  
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Scale")


## ---- INTERNAL CONSISTENCY ----

# calculate alphas for each scale
alphs <- lapply( set_names( names(scls) ), function(i) alpha( d1[ , scls[[i]] ], check.keys = T ) )

# prepare a table with internal consistency measures per (sub)scale
icons <-

  sapply( names(alphs), function(i) alphs[[i]]$total[ c("raw_alpha","std.alpha","G6(smc)","average_r") ] ) %>%
  t() %>%
  as.data.frame() %>%
  mutate( across( everything(), ~ unlist(.x, use.names = F) ) ) %>%
  rownames_to_column("Scale")


## --- RELIABILITY TABLE ----

write.table(
  x = left_join(cfind, icons, by = "Scale"),
  file = here("tables","vls_reliabilities.csv"),
  sep = ",",
  row.names = F,
  quote = F
)


# SUM SCORE DIFFERENCES ----

write.table(
  
  x = d1 %>%
    group_by(SA) %>%
    summarise( across( starts_with("VLS"), msd ) ) %>%
    pivot_longer(cols = -SA) %>%
    pivot_wider(names_from = SA, values_from = value) %>%
    rename("Scale" = "name") %>%
    
    left_join(
      
      sapply(
        
        names(d1)[ grepl( "VLS", names(d1) ) ],
        function(y) {
          
          # calculate t-test
          t <- t.test(as.formula( paste0(y," ~ SA") ), data = d1)
          d <- cohen.d(as.formula( paste0(y," ~ SA") ), data = d1)
          
          # calculate Mann-Whiteny test
          wilcox <- wilcox.test(as.formula( paste0(y," ~ SA") ), data = d1)
          vda <- VD.A(as.formula( paste0(y," ~ SA") ), data = d1)
          
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
        rownames_to_column("Scale")
      
    ),
  
  file = here("tables","vls_sum_scores.csv"),
  sep = ";",
  row.names = F,
  quote = F

)

