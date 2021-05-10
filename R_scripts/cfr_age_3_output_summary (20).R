# Output summary: Age stratified CFR

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Load functions
source("source/mcmc_summary.R")
source("source/mcmc_and_likelihood_functions.R")

# Read input raw data
data <- readRDS("data/age_disaggregated_onset_incidence_data.rds")
observed_deaths <- c(4, 0, 15, 65, 286, 1159, 3367, 8830, 19483)
observed_cases <- c(2099, 3744, 13673, 18755, 31057, 42704, 31777, 33916, 58317)
prop_lombardy <- 16466 / 33209
crude_cfr <- 33209 / 238082

# Read raw MCMC output and inputs used to run the MCMC
m <- readRDS("data/MCMC_fitting_output.rds")
# Data for prediction from likelihood function
data_list <- readRDS("data/predicted_x.rds")

### Processing MCMC output #####################################################
# Estimate age specific CFRs from reference group CFR and relative risks
for (i in grep("RR", names(m))) {
  new_name <- sub("RR", "cfr", names(m)[i])
  m$new <- unname(m$`cfr_80+`) * unlist(unname(m[, i]))
  names(m)[ncol(m)] <- new_name
}

# Remove burn in and non parameter columns
m_sampling <- filter(m, stage == "sampling") %>%
  select(-chain, -rung, -iteration, -stage, -loglikelihood)

# Thin
samples <- 1000
sample_index <- round(seq(1 , nrow(m_sampling), length.out = samples))
m_thin <- m_sampling[sample_index, ]

# Remove log likelihood and logprior columns
predict_par <- m_thin[, c("m_od", "s_od", "maxday", "z", "r", "D", "cfr_80+", "RR_0_9", "RR_10_19",
                          "RR_20_29", "RR_30_39", "RR_40_49", "RR_50_59", "RR_60_69", "RR_70_79")]

# Remove RR columns
m_thin <- select(m_thin, !contains("RR"))

# Summarise
posterior_summary <- summarise_mcmc(m_thin)

################################################################################

### Estimate Expected deaths, proportion of deaths in lombardy #####
fit_prediction <- apply(predict_par, 1, function(row_params){
  lL_plane(row_params, data_list)
})

ag <- unique(data$age_groups)
fit_prediction_deaths <- bind_rows(lapply(fit_prediction, function(x){
  out <- x$age_deaths
  d1 <- data.frame(t(out))
  colnames(d1) <- ag
  d1
  }))
fit_prediction_prop_lombardy <- unlist(lapply(fit_prediction, function(x){x$prop_lombardy}))
fit_prediction_crude_cdr <- unlist(lapply(fit_prediction, function(x){x$crude_cfr}))

fit_prediction_deaths_summary <- summarise_mcmc(fit_prediction_deaths) %>%
  mutate(observed = observed_deaths) %>%
  select(observed, everything())

fit_prediction_prop_lombardy_summary <- summarise_mcmc(fit_prediction_prop_lombardy) %>%
  mutate(observed = prop_lombardy) %>%
  select(observed, everything())

fit_prediction_crude_cfr_summary <- summarise_mcmc(fit_prediction_crude_cdr) %>%
  mutate(observed = crude_cfr) %>%
  select(observed, everything())
################################################################################

### Estimate adjusted cases and infections #####################################
cfr_names <- c("cfr_0_9", "cfr_10_19", "cfr_20_29", "cfr_30_39",
               "cfr_40_49", "cfr_50_59", "cfr_60_69", "cfr_70_79", "cfr_80+")
# Sum over time to reduce dimensionality
data2 <- data %>%
  group_by(location, age_groups) %>%
  summarise(cases = sum(cases),
            nici = unique(nici)) %>%
  ungroup()
# Find the z for "gold standard" cases ascertainment outside of lombardy
best_z <- min(data[data$location == "Outside", "nici"])
# Predict adjusted cases and infections | z, ni/ci
predicted <- bind_rows(lapply(1:nrow(m_thin), function(x, data2, m_thin, best_z){
  par <- m_thin[x,]
  data2 %>% mutate(
    age60 = ifelse(age_groups %in% c("0-9", "10-19", "20-29", "30-39",
                                     "40-49", "50-59"), "U60", "O60"),
    cases_adjusted = cases * ifelse(location == "Lombardy", par$z, 1 /best_z) * nici,
    infections = cases_adjusted,
    all_deaths = cases_adjusted * unlist(par[cfr_names]),
    draw = x)
}, data2 = data2, m_thin = m_thin, best_z = best_z))

################################################################################

#### CFR #######################################################################

# Overall
cfr <- predicted %>%
  group_by(draw) %>%
  summarise(cases_adjusted = sum(cases_adjusted),
            all_deaths = sum(all_deaths),
            CFR = all_deaths / cases_adjusted) %>%
  group_by() %>%
  summarise(CFR = list(summarise_mcmc(CFR))) %>%
  unnest(cols = c(CFR))

# By age
cfr_age <- predicted %>%
  group_by(draw, age_groups) %>%
  summarise(cases_adjusted = sum(cases_adjusted),
            all_deaths = sum(all_deaths),
            CFR = all_deaths / cases_adjusted) %>%
  group_by(age_groups) %>%
  summarise(CFR = list(summarise_mcmc(CFR))) %>%
  unnest(cols = c(CFR))

# By under 60
cfr_60 <- predicted %>%
  group_by(draw, age60) %>%
  summarise(cases_adjusted = sum(cases_adjusted),
            all_deaths = sum(all_deaths),
            CFR = all_deaths / cases_adjusted) %>%
  group_by(age60) %>%
  summarise(CFR = list(summarise_mcmc(CFR))) %>%
  unnest(cols = c(CFR))

################################################################################

### Unadjusted CFR #############################################################
# Overall
unadjusted_cfr <- predicted %>%
  group_by(draw) %>%
  summarise(cases = sum(cases),
            all_deaths = sum(all_deaths),
            unadjustedCFR = all_deaths / cases) %>%
  ungroup() %>%
  summarise(unadjustedCFR = list(summarise_mcmc(unadjustedCFR))) %>%
  unnest(cols = c(unadjustedCFR))

# By age
unadjusted_cfr_age <- predicted %>%
  group_by(draw, age_groups) %>%
  summarise(cases = sum(cases),
            all_deaths = sum(all_deaths),
            unadjustedCFR = all_deaths / cases) %>%
  group_by(age_groups) %>%
  summarise(unadjustedCFR = list(summarise_mcmc(unadjustedCFR))) %>%
  unnest(cols = c(unadjustedCFR))

# By over 60
unadjusted_cfr_60 <- predicted %>%
  mutate(age60 = ifelse(age_groups %in% c("0-9", "10-19", "20-29", "30-39",
                                          "40-49", "50-59"), "U60", "O60")) %>%
  group_by(draw, age60) %>%
  summarise(cases = sum(cases),
            all_deaths = sum(all_deaths),
            unadjustedCFR = all_deaths / cases) %>%
  group_by(age60) %>%
  summarise(unadjustedCFR = list(summarise_mcmc(unadjustedCFR))) %>%
  unnest(cols = c(unadjustedCFR))

# By location
unadjusted_cfr_location <- predicted %>%
  group_by(draw, location) %>%
  summarise(cases = sum(cases),
            all_deaths = sum(all_deaths),
            unadjustedCFR = all_deaths / cases) %>%
  group_by(location) %>%
  summarise(unadjustedCFR = list(summarise_mcmc(unadjustedCFR))) %>%
  unnest(cols = c(unadjustedCFR))

# By location X age
unadjusted_cfr_age_location <- predicted %>%
  group_by(draw, location, age_groups) %>%
  summarise(cases = sum(cases),
            all_deaths = sum(all_deaths),
            unadjustedCFR = all_deaths / cases) %>%
  group_by(location, age_groups) %>%
  summarise(unadjustedCFR = list(summarise_mcmc(unadjustedCFR))) %>%
  unnest(cols = c(unadjustedCFR))
################################################################################

### Plot CFR ###########################################################
m_cfr <- m_thin[,cfr_names]
cfr_age_pd <- m_cfr %>%
  pivot_longer(col = everything(), names_to = "Age", values_to = "CFR", names_prefix = "cfr_")
cfr_age_plot <- ggplot(cfr_age_pd, aes(x = Age, y = CFR)) + 
  geom_boxplot() + 
  theme_bw() +
  ylim(0, 0.3)

################################################################################

### Ascertainment fractions ####################################################
# Overall
af <- predicted %>%
  group_by(draw) %>%
  summarise(af_cases  = sum(cases) / sum(cases_adjusted),
            af_infections = sum(cases) / sum(infections)) %>%
  pivot_longer(cols = c(-draw), names_to = "type", values_to = "af", names_prefix = "af_") %>%
  group_by(type) %>%
  summarise(af = list(summarise_mcmc(af))) %>%
  unnest(col = c(af)) 

# AF by age
af_age <- predicted %>%
  group_by(draw, age_groups) %>%
  summarise(af_cases  = sum(cases) / sum(cases_adjusted),
            af_infections = sum(cases) / sum(infections)) %>%
  pivot_longer(cols = c(-draw, -age_groups), names_to = "type", values_to = "af", names_prefix = "af_") %>%
  group_by(type, age_groups) %>%
  summarise(af = list(summarise_mcmc(af))) %>%
  unnest(col = c(af)) 

# AF by location
af_location <- predicted %>%
  group_by(draw, location) %>%
  summarise(af_cases  = sum(cases) / sum(cases_adjusted),
            af_infections = sum(cases) / sum(infections)) %>%
  pivot_longer(cols = c(-draw, -location), names_to = "type", values_to = "af", names_prefix = "af_") %>%
  group_by(type, location) %>%
  summarise(af = list(summarise_mcmc(af))) %>%
  unnest(col = c(af)) 

# AF ~ Location X Age
af_location_age <- predicted %>%
  group_by(draw, location, age_groups) %>%
  summarise(af_cases  = sum(cases) / sum(cases_adjusted),
            af_infections = sum(cases) / sum(infections)) %>%
  pivot_longer(cols = c(-draw, -location, -age_groups), names_to = "type", values_to = "af", names_prefix = "af_") %>%
  group_by(type, location, age_groups) %>%
  summarise(af = list(summarise_mcmc(af))) %>%
  unnest(col = c(af)) 

af_plot <- ggplot(af_location_age, aes(x = age_groups, y = median, ymin = x_0.025,
                                       ymax = x_0.975, col = type)) + 
  geom_errorbar(width = 0) +
  geom_point() +
  ylab("Ascertainment fraction") + 
  xlab("Age") + 
  facet_grid(location ~ type) + 
  theme_bw()

af_case_location_age_pd <- af_location_age %>%
  ungroup() %>%
  filter(type == "cases") %>%
  mutate(location = ifelse(location == "Lombardy", location, "Rest of Italy"))

af_cases_plot <- ggplot(af_case_location_age_pd,
                        aes(x = age_groups, y = median)) + 
  geom_point() +
  ylab("Case ascertainment fraction") + 
  xlab("Age") + 
  facet_grid(~ location) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white"))
################################################################################

### Hospitalisation ############################################################
# From ISS
proportion_observed_cases_severe <- c(0.063348416, 0.052036199, 0.03280543, 0.082579186,0.090497738,0.097285068,0.136877828, 0.222850679,0.262443439)

proportion_infections_severe <- predicted %>%
  group_by(draw, location) %>%
  mutate(h1 = cases * proportion_observed_cases_severe) %>%
  group_by(draw, age_groups) %>%
  summarise(hosp = sum(h1) / sum(infections)) %>%
  group_by(age_groups) %>%
  summarise(hosp = list(summarise_mcmc(hosp))) %>%
  unnest(col = c(hosp))
proportion_infections_severe$mode[proportion_infections_severe$mode < 0 ] <- 0
# Comparative method
proportion_observed_cases_severe * filter(af_age, type == "infections")$mean
################################################################################

### Observed case age distribution #############################################
case_age_dist <- data %>%
  mutate(location = ifelse(location == "Lombardy", location, "Rest of\nItaly")) %>%
  group_by(location, age_groups) %>%
  summarise(cases = sum(cases)) %>%
  group_by(location) %>%
  mutate(proportion = cases / sum(cases)) %>%
  ungroup() %>%
  mutate(proportion = ifelse(location == "Lombardy", proportion, -proportion))
sum(case_age_dist$proportion)

################################################################################

### Onset to death dist sampled over m_od & s_od pairs #########################
# Uncertainty
o2dpred <- bind_rows(lapply(1:nrow(m_thin), function(x, m_thin){
  m_od <- unlist(m_thin[x, "m_od"])
  s_od <- unlist(m_thin[x,"s_od"])
  alpha <- 1 / (s_od * s_od)
  beta <- 1 / (m_od * s_od * s_od)
  d <- seq(0, 100, 0.1)
  y <- dgamma(d, alpha, beta)
  data.frame(x = d, y = y, draw = x)
}, m_thin = m_thin))
# Posterior mode params
m_od <- posterior_summary$mode[rownames(posterior_summary) == "m_od"]
s_od <- posterior_summary$mode[rownames(posterior_summary) == "s_od"]
alpha <- 1 / (s_od * s_od)
beta <- 1 / (m_od * s_od * s_od)
d <- seq(0, 100, 0.1)
y <- dgamma(d, alpha, beta)
o2d_pred_m <- data.frame(x = d, y = y)

o2d_plot <- ggplot() +
  geom_line(data = o2dpred, aes(x = x, y = y, group = draw), alpha = 0.05) + 
  geom_line(data = o2d_pred_m, aes(x = x, y = y), col = "red") + 
  xlim(0, 60) +
  xlab("Days") +
  ylab("P") +
  theme_bw()
################################################################################


### Observed Deaths/Cases CI ###################################################
study_data <- binom::binom.confint(observed_deaths, observed_cases, method = "exact")
################################################################################

# Output_table #################################################################
# Formatting output to go into manuscript tables
options(scipen = 50)

study_data_table <- data.frame(
  Deaths = study_data$x,
  Cases = study_data$n,
  crude_cfr = paste0(signif(100 * study_data$mean, 3), "% (",
                     signif(100 * study_data$lower, 3), ",", 
                     signif(100 * study_data$upper, 3), ")")
)

output_table_overall <- data.frame(
  "CFR adjusting for censoring" = 
    paste0(signif(100 * unadjusted_cfr$mode, 3), "% (",
           signif(100 * unadjusted_cfr$x_0.025, 3), ",", 
           signif(100 * unadjusted_cfr$x_0.975, 3), ")"),
  "CFR additionally adjusted for demography and under-ascertainment " = 
    paste0(signif(100 * cfr$mode, 3), "% (",
           signif(100 * cfr$x_0.025, 3), ",", 
           signif(100 * cfr$x_0.975, 3), ")"))

output_table_age <- data.frame(
  "Age" = unique(data$age_groups),
  "CFR adjusting for censoring" = 
    paste0(signif(100 * unadjusted_cfr_age$mode, 3), "% (",
           signif(100 * unadjusted_cfr_age$x_0.025, 3), ",", 
           signif(100 * unadjusted_cfr_age$x_0.975, 3), ")"),
  "CFR additionally adjusted for demography and under-ascertainment " = 
    paste0(signif(100 * cfr_age$mode, 3), "% (",
           signif(100 * cfr_age$x_0.025, 3), ",", 
           signif(100 * cfr_age$x_0.975, 3), ")"))
output_table_60 <- data.frame(
  "Age" = c("Over 60", "Under 60"),
  "CFR adjusting for censoring" = 
    paste0(signif(100 * unadjusted_cfr_60$mode, 3), "% (",
           signif(100 * unadjusted_cfr_60$x_0.025, 3), ",", 
           signif(100 * unadjusted_cfr_60$x_0.975, 3), ")"),
  "CFR additionally adjusted for demography and under-ascertainment " = 
    paste0(signif(100 * cfr_60$mode, 3), "% (",
           signif(100 * cfr_60$x_0.025, 3), ",", 
           signif(100 * cfr_60$x_0.975, 3), ")"))

################################################################################
cfr_age_plot 
af_cases_plot

