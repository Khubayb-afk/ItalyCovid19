####################################################################################
##                   Monday, 24 February 2020 - Tuesday 5 May 2020                ##
##  1) Data Collation and Preprocessing - An Overview                             ##
##                                                                                ##
##   This section loads in data on the demographic structure of the Italian       ##
##   population, as well as daily case data for Lombardy (Region) and rest of     ##
##   Italy from information from the Italian Health Deparment. This is integrated ##
##   with information on the Age distribution of cases from reports from the      ##
##   Italian Govt to calculate the age distribution of cases for Lombardy 	      ##
##   and the rest of Italy.                                                       ##
##                                                                                ##
##   Assuming an invariant age-distribution of cases over time, this is then      ##
##   used to convert the aggregate daily incidence of cases (by symptom onset)    ##
##   into age-specific daily case incidences.                                     ##
##                                                                                ##
####################################################################################

# Loading Libraries
library(tidyverse)

# Loading in Demographic and Case Data
population <- read.csv("data/italy_population_demography.csv") 
Epicentro_report_onset <- read.csv("data/Epicentro_onset_data.csv")

# Setting Static Variables
smoothing_centre <- as.Date("2020-03-26", format = "%Y-%m-%d")
days_to_smooth_either_side <- 3
age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

####################################################################################
##                                                                                ##
## Initial Preprocessing of Data Including 1st May Incidence Spike Smoothing      ##
##                                                                                ##
##    Ensuring data variables are in correct format, and calculating              ##
##    the incidence of cases outside Italy over time.                             ##
##    Also contains code to correct for incidence spike that occurred on          ##
##    May1st. Initial linear interpolation of case numbers between                ##
##    April 30th and 1st May, followed by distribution of cases                   ##
##    equally around 3 days either side.                                          ##
##                                                                                ##
####################################################################################

# Loading in Raw Case Onset Incidence Data
cases_by_onset <- Epicentro_report_onset %>%
  mutate(Date = as.Date(x = Date, format = "%d/%m/%Y")) %>%
  select(Date, Cases_Inside_Lombardy, Cases_Italy) %>%
  mutate(Cases_Outside_Lombardy = Cases_Italy - Cases_Inside_Lombardy)

# Smoothing of Incidence Spike Around 26th March
## Creating vector of dates to smooth over
total_days_smoothing <- 2 * days_to_smooth_either_side + 1
dates_vector <- c()
start_date <- smoothing_centre - days_to_smooth_either_side
for (i in 1:total_days_smoothing) {
  dates_vector[i] <- as.character(start_date)
  start_date <- start_date + 1
}

## Linear interpolation of cases on 26th March from cases on 25th March & 27th March
cases_25th_Mar <- cases_by_onset[cases_by_onset$Date == "2020-03-25", 2:4]
cases_26th_Mar <- cases_by_onset[cases_by_onset$Date == "2020-03-26", 2:4]
cases_27th_Mar <- cases_by_onset[cases_by_onset$Date == "2020-03-27", 2:4]
adj_cases_26th_Mar <- (cases_25th_Mar + cases_27th_Mar)/2 # linear intepolation

## Redstributing extra cases to surrounding 3 days either side & 1st May equally 
cases_by_onset[cases_by_onset$Date == "2020-03-26", 2:4] <- adj_cases_26th_Mar
cases_to_distribute <- (cases_26th_Mar - adj_cases_26th_Mar)/total_days_smoothing
rows_for_adjustment <- cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4]
for (i in 1:nrow(rows_for_adjustment)) {
  rows_for_adjustment[i, ] <- rows_for_adjustment[i, ] + cases_to_distribute
}
cases_by_onset[as.character(cases_by_onset$Date) %in% dates_vector, 2:4] <- rows_for_adjustment

####################################################################################
## Calculating Age Distribution of Cases Inside and Outside Lombardy              ##
##                                                                                ##
####################################################################################

# Age Distribution of Italy's Population
Italy_all_age <- population$proportion[population$location == "Outside"]
Italy_all_age <- Italy_all_age/100

# Age Distributions of Cases for Lombardy and Nationally (Including Lombardy)
Lombardy_case_age<- c(0.9, 1.6, 5.7, 7.9, 13.0, 17.9, 13.3, 14.2, 25.3) 
Lombardy_case_age<- Lombardy_case_age/100 # convert percentage to proportion
Italy_case_age <- c(0.9, 1.6, 5.7, 7.9, 13.0, 17.9, 13.3, 14.2, 25.3) 
Italy_case_age <- Italy_case_age/100 # convert percentage to proportion

# Calculating the Age Distribution of Cases Outside Lombardy
Lombardy_Cases_Onset <- cases_by_onset$Cases_Inside_Lombardy # Lombardy only
Italy_Cases_Onset <- cases_by_onset$Cases_Italy # All of Italy, including Lombardy
Outside_Lombardy_Cases <- sum(Italy_Cases_Onset) - sum(Lombardy_Cases_Onset)
Lombardy_cases_by_age <- sum(Lombardy_Cases_Onset) * Lombardy_case_age# cases by age for Lombardy
Total_cases_by_age <- sum(Italy_Cases_Onset) * Italy_case_age # cases by age for Italy (inc. Lombardy)
Outside_Lombardy_Cases_by_age <- Total_cases_by_age - Lombardy_cases_by_age # cases by age outside Lombardy
Outside_Lombardy_Case_Age_Dist <- Outside_Lombardy_Cases_by_age/Outside_Lombardy_Cases # proportion cases by age outside Lombardy


####################################################################################
## Calculating Age-Disaggregated Onset Incidence Over Time, for Inside &          ##
##  Outside Lombardy								                                              ##
##                                                                                ##
##    Using the daily onset information for Lombardy and outside Lombardy in      ##
##    conjunction with the age distribution of cases for each of these settings,  ##
##    we are able to calculate onset incidence over time disaggregated by         ##
##    age-group, for each location.                                               ##
##                                                                                ## 
####################################################################################

# Calculating the Incidence of Cases In Different Age Groups Over Time for Inside Lombardy
##  Create dataframe of Lombardy case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_Lombardy <- data.frame(age_groups, Lombardy_case_age)
Lombardy_Age_Cases_Time <- expand.grid(Lombardy_case_age, cases_by_onset$Cases_Inside_Lombardy)
colnames(Lombardy_Age_Cases_Time) <- c("Age_Prop", "Cases")
Lombardy_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(Lombardy_cases_by_age))
Lombardy_Age_Cases_Time <- Lombardy_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_Lombardy, by = c("Age_Prop" = "Lombardy_case_age")) %>%
  mutate(Location = "Lombardy")

# Calculating the Incidence of Cases In Different Age Groups Over Time for Outside Lombardy
##  Create dataframe of outside Lombardy case distribution by age
##  Create grid with all combinations of age group and daily case incidence
##  Add date to the dataframe
##  Multiply cases on each day for by proportion occurring in each age group
Age_Dist_df_Outside <- data.frame(age_groups, Outside_Lombardy_Case_Age_Dist)
Outside_Age_Cases_Time <- expand.grid(Outside_Lombardy_Case_Age_Dist, cases_by_onset$Cases_Outside_Lombardy)
colnames(Outside_Age_Cases_Time) <- c("Age_Prop", "Cases")
Outside_Age_Cases_Time$Date <- rep(cases_by_onset$Date, each = length(Outside_Lombardy_Cases_by_age))
Outside_Age_Cases_Time <- Outside_Age_Cases_Time %>%
  mutate(Cases = Cases * Age_Prop) %>%
  left_join(Age_Dist_df_Outside, by = c("Age_Prop" = "Outside_Lombardy_Case_Age_Dist")) %>%
  mutate(Location = "Outside")

# Combining the Dataframes Together
age_disaggregated_counts_df <- rbind(Lombardy_Age_Cases_Time, Outside_Age_Cases_Time)
age_disaggregated_counts_df <- age_disaggregated_counts_df[, c("Date", "Location", "age_groups", "Cases")]
colnames(age_disaggregated_counts_df) <- c("date", "location", "age_groups", "cases")

# Integrating Age-Disaggregated Cases With Population to Calculate Adjustment Factor
#   Used in downstream analyses
raw_adjustment_factor_df <- age_disaggregated_counts_df %>%
  group_by(age_groups, location) %>%
  summarise(cases = sum(cases)) %>%
  left_join(population, by = c("age_groups" = "age_groups", "location" = "location")) %>%
  mutate(nici = population/cases) %>%
  select(age_groups, location, nici)

# Joining the Case Incidence and Adjustment Factor Datasets Together
age_disaggregated_counts_df <- age_disaggregated_counts_df %>%
  left_join(raw_adjustment_factor_df, by = c("age_groups" = "age_groups", "location" = "location"))

# Saving Created Dataset
saveRDS(age_disaggregated_counts_df, file = "data/age_disaggregated_onset_incidence_data.rds")



