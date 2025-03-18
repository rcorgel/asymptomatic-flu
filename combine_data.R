################################################################################
# File Name: combine_data                                                      #
#                                                                              #
# Purpose:   Combine county-level claims data with external data sources.      #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Format data (impute, filter, collapse)                         #
#            4. Merge and save data                                            #
#                                                                              #
# Project:   Asymptomatic Influenza Associations                               #
################################################################################

####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(assertr)
library(usdata)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/asymptomatic-flu-proj/')

####################################
# 2. LOAD CLAIMS AND EXTERNAL DATA #
####################################

# Load claims and external data
disease_dat <- read.csv('raw/flu_2022-09_2023-08_p3.csv') # data with covid + rsv
flu_dat <- read.csv('raw/flu_2022-09_2023-08_p2.csv') # flu data with age/gender
race_dat <- read.csv('raw/flu_2022-09_2023-08_race_p2.csv') # flu data with race
all_cause_dat <- read.csv('raw/county_season_ac_v2.csv') # ppl that sought care
weather_dat <- read.csv('raw/tot_2023_crosswalk_county.csv') # humidity
county_dat <- read.csv('raw/analytic_data2024.csv') # county predictors
icu_dat <- read.csv('raw/data-FPBfZ.csv') # icu beds per county

######################################
# 3. FORMAT CLAIMS AND EXTERNAL DATA #
######################################

######################
# IMPUTE CLAIMS DATA #
######################

# Flu data
# Check the percent of rows that will be imputed (80.4%)
nrow(flu_dat[flu_dat$patient_count == "<=5",]) / nrow(flu_dat)
# Impute rows where patient count is <=5 to a random number 1 to 5
flu_dat$patient_count_imp <- flu_dat$patient_count
flu_dat$patient_count_imp[flu_dat$patient_count_imp == '<=5'] <- 
  sample(1:5, 
         length(flu_dat$patient_count_imp[flu_dat$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
flu_dat$patient_count_imp <- as.numeric(flu_dat$patient_count_imp)
# Check the percent of patients that were imputed (7.4%)
sum(flu_dat[flu_dat$patient_count_imp <= 5,]$patient_count_imp) / 
  sum(flu_dat$patient_count_imp)

# Disease data
# Check the percent of rows that will be imputed (77.8%)
nrow(disease_dat[disease_dat$patient_count == "<=5",]) / nrow(disease_dat)
# Impute rows where patient count is <=5 to a random number 1 to 5
disease_dat$patient_count_imp <- disease_dat$patient_count
disease_dat$patient_count_imp[disease_dat$patient_count_imp == '<=5'] <- 
  sample(1:5, 
         length(disease_dat$patient_count_imp[disease_dat$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
disease_dat$patient_count_imp <- as.numeric(disease_dat$patient_count_imp)
# Check the percent of patients that were imputed (2.8%)
sum(disease_dat[disease_dat$patient_count_imp <= 5,]$patient_count_imp) / 
  sum(disease_dat$patient_count_imp)

# Race data
# Check the percent of rows that will be imputed (77.9%)
nrow(race_dat[disease_dat$patient_count == "<=5",]) / nrow(race_dat)
# Impute rows where patient count is <=5 to a random number 1 to 5
race_dat$patient_count_imp <- race_dat$patient_count
race_dat$patient_count_imp[race_dat$patient_count_imp == '<=5'] <- 
  sample(1:5, 
         length(race_dat$patient_count_imp[race_dat$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
race_dat$patient_count_imp <- as.numeric(race_dat$patient_count_imp)
# Check the percent of patients that were imputed (4.1%)
sum(race_dat[race_dat$patient_count_imp <= 5,]$patient_count_imp) / 
  sum(race_dat$patient_count_imp)

###############
# FILTER DATA #
###############

# Flu data
flu_filt <- flu_dat |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  filter(!(patient_gender_code %in% c("U","")))
# Check the number of individuals kept (98.5%)
sum(flu_filt$patient_count_imp) / sum(flu_dat$patient_count_imp)
flu_filt |> assert(not_na, county_fips:state_fips)

# Disease data
disease_filt <- disease_dat |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99","")))
# Check the number of individuals kept (99.9%)
sum(disease_filt$patient_count_imp) / sum(disease_dat$patient_count_imp)
disease_filt |> assert(not_na, county_fips:state_fips)

# Race data
race_filt <- race_dat |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing races
  filter(!(race_code %in% c("U","")))
# Check the number of individuals kept (54.5%)
sum(race_filt$patient_count_imp) / sum(race_dat$patient_count_imp)
race_filt |> assert(not_na, county_fips:state_fips)

# Weather data
weather_filt <- weather_dat |>
  # Add leading 0's to fips code
  mutate(county_fips = ifelse(nchar(county_fips) == 4, 
                              paste0("0", county_fips), county_fips)) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99","78", "60", "99999", ""))) |>
  # Missing humidity or temp
  filter(!is.na(temp)) |>
  filter(!is.na(humidity))
# Check the number of rows kept (97.3%)
length(weather_filt$county_fips) / length(weather_dat$county_fips)
# Some days are missing

# County data
county_filt <- county_dat[-1, ]
county_filt <- county_filt |>
  rename('county_fips' = 'X5.digit.FIPS.Code') |>
  filter(!(State.FIPS.Code %in% c("72","99","78", "60", "99999", "", "00")))
# Check the number of rows kept (99.9%)
length(county_filt$county_fips) / length(county_dat$X5.digit.FIPS.Code)

# ICU data
# No county_fips so merge it on
fips_merge <- county_filt |> dplyr::select(c(county_fips, State.Abbreviation, Name)) |>
  mutate(state_name = abbr2state(State.Abbreviation),
         county_name = gsub('County','', Name),
         county_name = gsub('Parish','', county_name),
         county_name = gsub('Census Area','', county_name),
         county_name = trimws(county_name),
         county_name = gsub( " ", "", tolower(county_name)))
icu_dat$County <- gsub( " ", "", tolower(icu_dat$County))
icu_filt <- left_join(icu_dat, fips_merge[, c(1, 4, 5)],
                      by = c('State' = 'state_name',
                             'County' = 'county_name'))
# Manually fill in some missing fips
icu_filt[icu_filt$State == 'District of Columbia',]$county_fips <- '11001'
icu_filt[icu_filt$County == 'anchorage',]$county_fips <- '02020'
icu_filt <- icu_filt |> 
  dplyr::select(c(county_fips, ICU.Beds)) |>
  dplyr::rename('icu_beds' = 'ICU.Beds')

#################################
# COLLAPSE DATA TO COUNTY LEVEL #
#################################

# Flu data
flu_vars <- flu_filt |> mutate(symp_count = fever + myalgia + cough + 
                     sore_throat + short_breath + hypoxemia + 
                     chest_pain + bronchitis + nausea_vom + 
                     diarrhea + fatigue + headache + congestion + sneezing) |>
  mutate(asymp_flu = ifelse((symp_count == 0 & flu == 1), 1, 0)) |>
  group_by(county_fips) |>
  mutate(total_flu = sum(flu),
         total_asymp_flu = sum(asymp_flu),
         percent_asymp_flu = total_asymp_flu/total_flu) |>
  distinct(county_fips, total_flu, total_asymp_flu, percent_asymp_flu, .keep_all = FALSE)

# Flu age
age_vars <- flu_filt |> filter(flu == 1) |>
  mutate(age_grp_count = 1) |>
  group_by(county_fips, age_grp) |>
  mutate(age_grp_total = sum(age_grp_count)) |>
  distinct(county_fips, age_grp, age_grp_total, .keep_all = FALSE) |>
  ungroup() |> group_by(county_fips) |>
  mutate(total = sum(age_grp_total),
         percent = age_grp_total/total) |>
  ungroup() |>
  dplyr::select(c(county_fips, age_grp, percent)) |>
  pivot_wider(names_from = age_grp, values_from = percent)
age_vars[c("0", "1", "2", "3", "4", "5")][is.na(age_vars[c("0", "1", "2", "3", "4", "5")])] <- 0
age_vars <- age_vars |>
  rename('percent_0_4_flu' = '0',
         'percent_5_12_flu' = '1',
         'percent_13_17_flu' = '2',
         'percent_18_49_flu' = '3',
         'percent_50_64_flu' = '4',
         'percent_65_plus_flu' = '5')

# Flu gender
gender_vars <- flu_filt |> filter(flu == 1) |>
  mutate(gender_count = 1) |>
  group_by(county_fips, patient_gender_code) |>
  mutate(gender_total = sum(gender_count)) |>
  distinct(county_fips, patient_gender_code, gender_total, .keep_all = FALSE) |>
  ungroup() |> group_by(county_fips) |>
  mutate(total = sum(gender_total),
         percent_female_flu = gender_total/total) |>
  filter(patient_gender_code == 'F') |>
  ungroup() |>
  dplyr::select(c(county_fips, percent_female_flu))

# Flu race
race_vars <- race_filt |> filter(flu == 1) |>
  mutate(race_count = 1) |>
  group_by(county_fips, race_code) |>
  mutate(race_total = sum(race_count)) |>
  distinct(county_fips, race_code, race_total, .keep_all = FALSE) |>
  ungroup() |> group_by(county_fips) |>
  mutate(total = sum(race_total),
         percent = race_total/total) |>
  ungroup() |>
  dplyr::select(c(county_fips, race_code, percent)) |>
  pivot_wider(names_from = race_code, values_from = percent)
race_vars[c("W", "B", "H", "A")][is.na(race_vars[c("W", "B", "H", "A")])] <- 0
race_vars <- race_vars |>
  rename('percent_white_flu' = 'W',
         'percent_black_flu' = 'B',
         'percent_hispanic_flu' = 'H',
         'percent_asian_flu' = 'A')

# COVID & RSV
disease_vars <- disease_filt |> group_by(county_fips) |>
  mutate(total_covid = sum(covid),
         total_rsv = sum(rsv)) |>
  distinct(county_fips, total_covid, total_rsv, .keep_all = FALSE)

# All cause
all_cause_vars <- all_cause_dat |> 
  filter(season == '2022-2023') |> 
  rename('total_all_cause' = 'all_cause') |>
  dplyr::select(c(county_fips, total_all_cause)) 

# Weather
weather_vars <- weather_filt |> group_by(county_fips) |>
  mutate(avg_temp = mean(temp, na.rm = TRUE),
         avg_humidity = mean(humidity, na.rm = TRUE)) |>
  distinct(county_fips, avg_temp, avg_humidity, .keep_all = FALSE)

# County (already at the county level)
county_vars <- county_filt |> 
  dplyr::select(c(county_fips, Poor.or.Fair.Health.raw.value,
                  Poor.Mental.Health.Days.raw.value, Adult.Obesity.raw.value,
                  Primary.Care.Physicians.raw.value, Unemployment.raw.value,
                  Life.Expectancy.raw.value, Frequent.Mental.Distress.raw.value,
                  Median.Household.Income.raw.value, Population.raw.value,
                  X..Rural.raw.value, Flu.Vaccinations.raw.value)) |>
  dplyr::rename('percent_poor_health' = 'Poor.or.Fair.Health.raw.value',
                'days_poor_mental_health' = 'Poor.Mental.Health.Days.raw.value',
                'percent_obese' = 'Adult.Obesity.raw.value',
                'rate_doctors' = 'Primary.Care.Physicians.raw.value',
                'percent_unemployed' = 'Unemployment.raw.value',
                'life_expectancy' = 'Life.Expectancy.raw.value',
                'percent_mental_distress' = 'Frequent.Mental.Distress.raw.value',
                'median_hh_income' = 'Median.Household.Income.raw.value',
                'population' = 'Population.raw.value',
                'percent_rural' = 'X..Rural.raw.value',
                'percent_vaccinated' = 'Flu.Vaccinations.raw.value') |>
  filter(!(percent_poor_health %in% c(""))) |>
  filter(!(days_poor_mental_health %in% c(""))) |>
  filter(!(percent_obese %in% c(""))) |>
  filter(!(rate_doctors %in% c(""))) |>
  filter(!(percent_unemployed %in% c(""))) |>
  filter(!(life_expectancy %in% c(""))) |>
  filter(!(percent_mental_distress %in% c(""))) |>
  filter(!(median_hh_income %in% c(""))) |>
  filter(!(population %in% c(""))) |>
  filter(!(percent_rural %in% c(""))) |>
  filter(!(percent_vaccinated %in% c(""))) |>
  mutate(percent_poor_health = as.numeric(percent_poor_health),
         days_poor_mental_health = as.numeric(days_poor_mental_health),
         percent_obese = as.numeric(percent_obese),
         rate_doctors = as.numeric(rate_doctors),
         percent_unemployed = as.numeric(percent_unemployed),
         life_expectancy = as.numeric(life_expectancy),
         percent_mental_distress = as.numeric(percent_mental_distress),
         median_hh_income = as.numeric(median_hh_income),
         population = as.numeric(population),
         percent_rural = as.numeric(percent_rural),
         percent_vaccinated = as.numeric(percent_vaccinated))

# ICU (done in previous section)
icu_vars <- icu_filt

##########################################
# 4. MERGE/SAVE CLAIMS AND EXTERNAL DATA #
##########################################

################
# COMBINE DATA #
################

# Merge on age
analysis_dat <- left_join(flu_vars, age_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on gender
analysis_dat <- left_join(analysis_dat, gender_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on race
analysis_dat <- left_join(analysis_dat, race_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on disease
analysis_dat <- left_join(analysis_dat, disease_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on all cause
analysis_dat <- left_join(analysis_dat, all_cause_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on weather
analysis_dat <- left_join(analysis_dat, weather_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on county
analysis_dat <- left_join(analysis_dat, county_vars, 
                          by = c('county_fips' = 'county_fips'))
# Merge on icu
analysis_dat <- left_join(analysis_dat, icu_vars, 
                          by = c('county_fips' = 'county_fips'))
#############
# SAVE DATA #
#############

write.csv(analysis_dat, "tmp/analysis_dat.csv")

################################################################################
################################################################################
