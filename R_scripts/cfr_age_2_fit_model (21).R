####################################################################################
##                                                                                ##
##  2) Model Fitting and Parameter Inference                                      ##
##                                                                                ##
##    The age-disaggregated incidence data for Lombardy and outside is then    	  ##
##    integrated with age-disaggregated deaths from Lombardy and outside over 	  ##
##    the same time period. We then fit an age-specific CFR model to these data,  ##
##    which explicitly includes considerations of under-reporting related to      ##
##    age-specific patterns of disease severity as well as the capacity of local  ##
##    health systems.                                                             ##
##                                                                                ##
####################################################################################

# Loading Libraries
set.seed(101092)
source("source/mcmc_and_likelihood_functions21.R")

library(devtools)
#devtools::install_github("mrc-ide/drjacoby", ref = "version1.0")  # uncomment and run this line the first time through
library(drjacoby)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Loading in Age and Location Disaggregated Case Data
age_disaggregated_case_onset <- readRDS("data/age_disaggregated_onset_incidence_data21.rds")
age_cases <- age_disaggregated_case_onset %>%
  group_by(age_groups) %>%
  summarise(counts = sum(cases))
data <- age_disaggregated_case_onset %>%
  arrange(location, age_groups, date)

# Setting Static Variables
age_groups <- c("0_9", "10_19", "20_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80+")
n_age_bands <- length(unique(data$age_groups))
max_date <- as.numeric(max(data$date) - min(data$date)) # most recent date in the date, relative to date start
observed_deaths <- c(10, 11 , 51 , 183, 810 , 3139, 9163, 23524, 60086) # https://www.epicentro.iss.it/coronavirus/bollettino/Bollettino-sorveglianza-integrata-COVID-19_10-febbraio-2021.pdf
prop_deaths_lombardy <- 28738 / 96977
deaths_lombardy <- round(prop_deaths_lombardy * sum(observed_deaths))
deaths_outside <- (1 - prop_deaths_lombardy) * sum(observed_deaths)
death_observation_censoring <- as.numeric(as.Date("2019-12-01") - min(data$date)) # Assume no deaths detecte before 21st Jan in Italy CDC numbers
deaths_x <- round(prop_deaths_lombardy * sum(observed_deaths))
deaths_n <- sum(observed_deaths)
# iss-covid19-bollettino-20210303:
ISS_observed_deaths <- 96977
ISS_observed_cases <- 2953120

# Input data
x <- c(0, n_age_bands,
       death_observation_censoring,
       min(data[data$location == "Outside", "nici"]),
       deaths_x, deaths_n,
       ISS_observed_deaths, ISS_observed_cases,
       observed_deaths,
       data$cases, as.numeric(as.factor(data$age_groups)),
       data$date - min(data$date), as.numeric(data$location == "Lombardy"),
       data$nici)

# Input data for output = prediction
x_output <- x
x_output[1] <- 1
saveRDS(x_output, "data/predicted_x21.rds")

# Parameters
#   z scales Lombardy relative to outside
#   r growth rate (fixed)
#   D detection window (fixed)
df_paramsRR <- data.frame(name = c("m_od", "s_od", "maxday", "z", "r", "D", "cfr_80+", "RR_0_9", "RR_10_19", "RR_20_29", "RR_30_39", "RR_40_49", "RR_50_59", "RR_60_69", "RR_70_79"),
                          min = c(10, 0, max_date, 0, 0, 0, rep(0, 9)),
                          max = c(Inf, Inf, max_date, 1, 0.1, 14, 1, rep(Inf, 8)),
                          init = c(15, 0.35, max_date, 0.005, 0.05, 10, 0.1, rep(1, 8)))
params <- df_paramsRR$init
lL_plane(df_paramsRR$init, x)

### Run MCMC ###################################################################
burnin <- 5000
samples <- 200
mcmc_output <- run_mcmc(data = x,
                        df_params = df_paramsRR,
                        loglike = lL_plane,
                        logprior = lP_plane,
                        burnin = burnin,
                        samples = samples,
                        chains = 1)
# (Run with chains > 1 for convergence checks)
saveRDS(mcmc_output, "data/complete_output21.rds")
saveRDS(mcmc_output$output, "data/MCMC_fitting_output21.rds") # save the MCMC output
saveRDS(list(inputs = x, 
             parameters = df_paramsRR,
             samples = samples,
             burnin = burnin), "data/MCMC_inputs_and_parameters21.rds") # save the data list used to run the MCMC
