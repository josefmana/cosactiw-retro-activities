# Run this script first, it prepares data for further analyses.

rm( list = ls() ) # clear environment

library(here)
library(openxlsx)
library(tidyverse)

mfile <- here("_raw","COSACTIW_DCERA-aktivity2.xlsx") # file with mapping
dfile <- here("_raw","COSACTIW_NANOK_pro-jamovi.xlsx") # file with data

# activities mapping
map <-
  read.xlsx( mfile, sheet = "Sheet 1" ) %>%
  select(type, category, activity) %>%
  add_row(type = "physical", category = "flexibility_health_exercise", activity = "chi_kung") %>%
  arrange(type, category, activity)

# demography
dem <-
  read.xlsx( dfile, sheet = "1ID_DATE_AGE_SCREENING_Demogr" ) %>%
  #filter( `Was_the_examination_valid?` == 1 ) %>%
  rename(
    "Age_years" = "Age",
    "Education_years" = "Number_of_years_of_study",
    "Education_level" = "Highest_education_level",
    "Subjective_age_years" = "Subjective_age_number",
    "Valid" = "Was_the_examination_valid?"
  ) %>%
  mutate(
    ID = gsub(" ", "", ID), # no spaces in ID allowed
    PA = factor(
      PA,
      levels = 1:6,
      labels = c(
        "I'm an athletic person",
        "I enjoy movement/exercise",
        "I exercised at least 3 times a week",
        "I don't avoid movement/exercise",
        "I'm not an athletic person",
        "I had to stop doing sports (Injury)"
      )
    ),
    Education_level = factor(
      Education_level,
      levels = 1:4,
      labels = c(
        "Primary school",
        "Vocational school (OU)",
        "Secondary school",
        "College/university"
      ),
      ordered = T
    )
  ) %>%
  select(ID, Valid, Age_years, Subjective_age_years, Education_years, Education_level, PA, FAQ, MMSE)

# cognition
cog <-
  read.xlsx( dfile, sheet = "3COGNITIVE_TESTS" ) %>%
  mutate(
    ID = gsub(" ", "", ID), # no spaces in ID allowed,
    SA = factor(
      `SA_New-BNT-TMT`,
      levels = 1:2,
      labels = c(
        "SA",
        "nonSA"
      )
    ),
    WHO_PA = factor(
      `WHO-PA`,
      levels = 0:1,
      labels = c(
        "nonPA",
        "PA"
      )
    )
  ) %>%
  select(ID, SA, WHO_PA, `RAVLT_1-5`, RAVLT_delayed_recall, TMT_A_time, TMT_B_time, Spon_sem, VF_animals)

# activities
act <-
  read.xlsx( dfile, sheet = "5RETROS-LEISURE_ACTIVITIES" ) %>%
  mutate(
    ID = gsub(".0", "", gsub(" ", "", ID), fixed = T), # no spaces or decimals in ID allowed
    Activity = if_else(Activity == "theathre", "theatre", Activity),
    Activity_type = case_when(
      Activity == "reading" ~ "mental",
      Activity == "theatre" ~ "mental",
      Activity == "self_education" ~ "mental",
      .default = Activity_type
    ),
    Seasonal = if_else(Intensity >= 10, T, F),
    Intensity = factor(
      Intensity %% 10, # modulo 10 to pool seasonal
      levels = c(0:5),
      labels = c(
        "do not know",
        "occasionally or not at all",
        "several times a month",
        "once a week",
        "several times a week",
        "every day or nearly every day"
      )
    ),
    Category = unlist( sapply( 1:nrow(.), function(i) map[map$activity == Activity[i], "category"] ), use.names = F )
  )

# save the outcomes as .rds
saveRDS( list(dem = dem, cog = cog, act = act, map = map), file = here("_data.rds") )
