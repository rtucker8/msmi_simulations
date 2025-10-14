#Purpose: Compare methods to compute state occupation probabilities from an illness-death model
#Author: Rachel Gonzalez
#Date: Sep 11 2025

#Libraries
library(tidyverse)
library(survival)
library(mstate)
library(msmi)
library(philentropy) #calculate KL divergence

# Functions ---------------------------------------------------------------
setwd("/Users/rachelgonzalez/Documents/Dissertation/msmi_simulations")
source("simulation_helper_V2.R")
load("Truth/sim_setting_characteristics.RData") #true state occupation probabilities for each setting

# Run Simulation ----------------------------------------------------------

###############################
## Set Simulation Parameters ##
###############################

#simulations setting corresponding to row in settings dataframe

if (!(Sys.getenv('SLURM_ARRAY_TASK_ID') == "")) {
  j = Sys.getenv('SLURM_ARRAY_TASK_ID') #.slurm script is an array job that runs each simulation setting
  j = as.numeric(j)
} else if (Sys.getenv('SLURM_ARRAY_TASK_ID') == ""){
  print("Simulations being run locally. Please make sure the right value of j is being used.")
  j = 1 #change this value to run different settings when running locally
}


#sample size
sample_size=100

#shape parameter for each transition
s01 <- settings[which(settings$setting == j), "s01"] #healthy to ill
s02 = 1.5 #healthy to dead
s12 = 2 #ill to dead

#median parameter for each transition
m01 <- settings[which(settings$setting == j), "m01"] #healthy to ill
m02 <- settings[which(settings$setting == j), "m02"] #healthy to dead
m12 <- settings[which(settings$setting == j), "m12"] #ill to dead

#log HR for extended model
b <- settings[which(settings$setting == j), "b"]

#upper bound for censoring
theta = 12

#evaluate methods at evenly spaced times between 0 and theta
eval_times = seq(0, theta, length.out = theta+1)
eval_times = eval_times[eval_times != 0] #remove time = 0

#true state occupation probabilities from sim_setting_characteristics.R
truth_all <- empirical_truth[[j]] %>% pivot_wider(names_from = state, values_from = probability)
truth <- truth_all %>% filter(time %in% eval_times)

#number of imputations for MI methods
number.nests <- 10
number.imps <- 20


################################
## Simulate Multiple Datasets ##
################################

#Simulate 500 datasets with reproducible random seeds
random.seeds <- read_csv("randomSeeds.csv")$simulationSeeds
d.sim <- map(random.seeds , function(s) {
  simulate_illness_death(n=sample_size, beta=b,
                         shape01 = s01, median01 = m01,
                         shape02 = s02,  median02 = m02,
                         shape12 = s12, median12 = m12,
                         return_latent_data = FALSE, theta=theta, seed = s)
})

# #Check semi-markov assumption on each of the 500 datasets
# #output a dataframe with columns for HR and p-value for t1 covariate
# semi_markov_check <- map_dfr(d.sim, function(data) {
#   ill <- subset(data, event1 == 1)
#   cox_semi_markov_test <- coxph(Surv(t2-t1, event2) ~ t1, data = ill)
#   return(summary(cox_semi_markov_test)$coefficients[,c("exp(coef)", "Pr(>|z|)")])
# })
#
# names(semi_markov_check) <- c("HR", "p_value")
#
# ggplot(semi_markov_check) + geom_histogram(aes(x=HR), bins=30) +
#   labs(title="Distribution of Hazard Ratios for t1 Covariate in Semi-Markov Datasets",
#        x="Hazard Ratio", y="Count") +
#   theme_minimal()
#
# #should be around 5% significant if semi-markov assumption holds
# semi_markov_check %>% summarise(significant_proportion = mean(p_value < 0.05),
#                                 HR_mean = mean(HR))

################################
##   Method: Aalen-Johansen   ##
################################

#create transition matrix for illness death model
tmat <- transMat(x = list(c(2,3), c(3), c()),
                 names= c("Healthy", "Ill", "Death"))

#function to compute AJ state occupation probabilities for a given dataset
aj_estimation <- function(df) {

  #prep data
  df.long <- msprep(data = df, trans = tmat,
                    time = c(NA, "t1", "t2"),
                    status = c(NA, "event1", "event2"))

  #Markov model without covariates- fully nonparametric
  c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=df.long, method="breslow")
  msf0 <- msfit(object = c0, vartype="aalen", trans = tmat)

  #state occupation probabilities
  pt0 <- probtrans(msf0, predt = 0, variance=FALSE)[[1]]

  return(pt0)
}

aj_probabilities <- map(d.sim, aj_estimation)

#AJ estimates for each state at eval_times
aj_probabilities <- map(aj_probabilities, function(df) {
  tibble(
    time = eval_times,
    pstate1 = map_dbl(eval_times, ~get_step_value(df$time, df$pstate1, .x)),
    pstate2 = map_dbl(eval_times, ~get_step_value(df$time, df$pstate2, .x)),
    pstate3 = map_dbl(eval_times, ~get_step_value(df$time, df$pstate3, .x))
  )
})

aj_results <- bind_rows(aj_probabilities, .id = "simulation") %>%
  select(simulation, time, pstate1, pstate2, pstate3)

#distance between AJ and truth over time
aj_distance_results <- map(aj_probabilities, ~compute_distance(truth, .x, cols_reference = c("Healthy", "Ill", "Dead"), times = eval_times))
aj_distance_results <- bind_rows(aj_distance_results, .id = "simulation") %>%
  select(simulation, time, KL, hellinger) %>% mutate(time = factor(time), method = "Aalen-Johansen")

#bias of AJ estimation over time
aj_results_plot <- aj_results %>% left_join(truth, by = "time") %>%
  mutate(pstate1_bias = pstate1 - Healthy,
         pstate2_bias = pstate2 - Ill,
         pstate3_bias = pstate3 - Dead,
         time = factor(round(time, 2))) %>%
  select(time, simulation, pstate1_bias, pstate2_bias, pstate3_bias, pstate1, pstate2, pstate3) %>%
  pivot_longer(
    cols = c(pstate1_bias, pstate2_bias, pstate3_bias,
             pstate1, pstate2, pstate3),
    names_to = c("state", "type"),
    names_pattern = "pstate(\\d+)(_bias)?",
    values_to = "value"
  ) %>%
  mutate(type = ifelse(type == "_bias","bias", "estimate")) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  select(time, simulation, state, bias, estimate) %>%
  filter(time != 0)


################################
##   Method: Marginal MI      ##
################################

imps_marginal <- map2(d.sim, random.seeds, function(d, s) {
  msmi.nested.impute(dat = d, M = number.nests, R=number.imps,  method="marginal", seed=s) %>% list_flatten()
})

#Obtain empirical state occupation probabilities for each imputed dataset (MxR of them)
empirical_probs_list <- map(imps_marginal, function(imp_list) {
  map(imp_list, function(df) {
    get_empirical_probs(df, eval_times)})
})

#Average empirical probabilities across imputations for each simulated dataset (Rubin's Rules Estimate)
mi_estimate <- map(empirical_probs_list, function(imp_list) {
  bind_rows(imp_list, .id = "imputation") %>%
    group_by(time) %>%
    summarise(across(starts_with("p"), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
})

mi_results <- bind_rows(mi_estimate, .id = "simulation") %>%
  select(simulation, time, pHealthy, pIll, pDead)

#Distance between marginal MI and truth over time
mi_distance_results <- map(mi_estimate, ~compute_distance(truth, .x, cols_estimate =c("pHealthy", "pIll", "pDead"), cols_reference = c("Healthy", "Ill", "Dead"), times = eval_times))
mi_distance_results <- bind_rows(mi_distance_results, .id = "simulation") %>%
  select(simulation, time, KL, hellinger)  %>% mutate(time = factor(time), method = "Marginal MI")

#Bias of marginal MI estimates over time
mi_results_plot <- mi_results %>% rename(pstate1 = pHealthy, pstate2 = pIll, pstate3=pDead) %>%
  left_join(truth, by = "time") %>%
  mutate(pstate1_bias = pstate1 - Healthy,
         pstate2_bias = pstate2- Ill,
         pstate3_bias = pstate3 - Dead,
         time = factor(round(time, 2))) %>%
  select(time, simulation, pstate1_bias, pstate2_bias, pstate3_bias, pstate1, pstate2, pstate3) %>%
  pivot_longer(
    cols = c(pstate1_bias, pstate2_bias, pstate3_bias,
             pstate1, pstate2, pstate3),
    names_to = c("state", "type"),
    names_pattern = "pstate(\\d+)(_bias)?",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(type == "_bias","bias", "estimate")
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>%
  select(time, simulation, state, bias, estimate) %>%
  filter(time != 0)

################################
##      Method: Cox MI        ##
################################

imps_cox <- map2(d.sim, random.seeds, function(d, s) {
 msmi.nested.impute(dat = d, M = number.nests, R=number.imps,  method="cox", seed=s) %>% list_flatten()
})

#Obtain empirical state occupation probabilities for each imputed dataset (MxR of them)
cox_empirical_probs_list <- map(imps_cox, function(imp_list) {
  map(imp_list, function(df) {get_empirical_probs(df, eval_times)})
})

#Combine empirical probabilities across imputations for each simulated dataset (sum over both M and R)
cox_mi_estimate <- map(cox_empirical_probs_list, function(imp_list) {
  bind_rows(imp_list, .id = "imputation") %>%
    group_by(time) %>%
    summarise(across(starts_with("p"), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
})

cox_mi_results <- bind_rows(cox_mi_estimate, .id = "simulation") %>%
  select(simulation, time, pHealthy, pIll, pDead)

#Distance between Cox MI and truth over time
cox_mi_distance_results <- map(cox_mi_estimate, ~compute_distance(truth, .x, cols_reference = c("Healthy", "Ill", "Dead"), cols_estimate =c("pHealthy", "pIll", "pDead"), times = eval_times))
cox_mi_distance_results <- bind_rows(cox_mi_distance_results, .id = "simulation") %>%
  select(simulation, time, KL, hellinger)  %>% mutate(time = factor(time), method = "Cox MI")

#Bias of Cox MI estimates over time
cox_mi_results_plot <- cox_mi_results %>% rename(pstate1 = pHealthy, pstate2 = pIll, pstate3=pDead) %>%
  left_join(truth, by = "time") %>%
  mutate(pstate1_bias = pstate1 - Healthy,
         pstate2_bias = pstate2- Ill,
         pstate3_bias = pstate3 - Dead,
         time = factor(round(time, 2))) %>%
  select(time, simulation, pstate1_bias, pstate2_bias, pstate3_bias, pstate1, pstate2, pstate3) %>%
  pivot_longer(
    cols = c(pstate1_bias, pstate2_bias, pstate3_bias,
             pstate1, pstate2, pstate3),
    names_to = c("state", "type"),
    names_pattern = "pstate(\\d+)(_bias)?",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(type == "_bias","bias", "estimate")
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>%
  select(time, simulation, state, bias, estimate) %>%
  filter(time != 0)

################################
##       Compare Methods      ##
################################

#Distance Figure
distance_results <- bind_rows(aj_distance_results, mi_distance_results, cox_mi_distance_results) %>%
  pivot_longer(cols = c(KL, hellinger), names_to = "distance", values_to = "value")

ggplot(data=distance_results) + geom_boxplot(aes(x=time, y=value, fill=method ), outliers= F) +
  labs(title="Distance of State Occupation Probability Estimates from Empirical Probabilities",
       x="Time Point", y="Distance Measure") +
  ylim(0, 0.3) +
  scale_fill_brewer(type="qual", palette = 6) +
  facet_wrap(~distance, nrow=2, ncol=1)
ggsave(paste0("Output/Figures/Distance/distance_results_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".pdf"),
       width=10, height=8)

#Create and save bias dataset
bias_results_plot <- bind_rows(aj_results_plot %>% mutate(method = "Aalen-Johansen"),
                               mi_results_plot %>% mutate(method = "Marginal MI"),
                               cox_mi_results_plot %>% mutate(method = "Cox MI"))
truth <- truth %>% pivot_longer(cols = c("Healthy", "Ill", "Dead"), names_to = "state", values_to = "truth") %>%
  mutate(state = case_when(state == "Healthy" ~ "1",
                           state == "Ill" ~ "2",
                           state == "Dead" ~ "3"),
         time = factor(round(time, 2)))
bias_results_plot <- bias_results_plot %>% left_join(truth, by = c("time", "state"))
write_csv(bias_results_plot, paste0("Output/bias_data_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".csv"))

#Bias Figure 1 (Faceted Boxplots)
ggplot(bias_results_plot) + geom_boxplot(aes(x=time, y=bias, fill=method), outliers = F) +
  geom_hline(yintercept = 0, color="grey70", linetype="dashed") + facet_wrap(~state, nrow=3, ncol=1) +
  labs(title="Bias of State Occupation Probability Estimates",
       x="Time", y="Bias") +
  scale_fill_brewer(type = "qual", palette = 6)

ggsave(paste0("Output/Figures/Bias/bias_results_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".pdf"),
       height = 10, width = 8)

#Bias Figure 2 (Smooth Curves with interquartile bands)
bias_summary <- bias_results_plot %>% group_by(method, state, time) %>%
  summarise(mean_estimate = mean(estimate),
            q1_estimate = quantile(estimate, .25),
            q3_estimate = quantile(estimate, .75),
            .groups = "drop") %>%
  mutate(state = case_when(state == "1" ~ "Healthy",
                           state == "2" ~ "Ill",
                           state == "3" ~ "Dead"))
truth <- truth %>% mutate(state = case_when(state == "1" ~ "Healthy",
                                            state == "2" ~ "Ill",
                                            state == "3" ~ "Dead"))
ggplot() +
  geom_line(data = truth, aes(x=time, y=truth), linetype = "dashed") +
  geom_line(data = bias_summary, aes(x=as.numeric(time), y=mean_estimate, color=method, group=method), linewidth=1.2) +
  geom_ribbon(data = bias_summary, aes(x=as.numeric(time), ymin=q1_estimate, ymax=q3_estimate, fill=method), alpha=0.2) +
  facet_wrap(~factor(state, levels = c("Healthy", "Ill", "Dead"), ordered=TRUE))
ggsave(paste0("Output/Figures/Bias/bias_ribbon_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".pdf"),
       width = 6.5, height = 4)

#Variation Figure
#compute monte carlo standard deviation and compare between methods
se_empiric <- bias_results_plot %>% group_by(method, state, time) %>%
  summarise(se = sd(estimate), .groups = "drop")

ggplot(se_empiric) + geom_point(aes(x = time, y = se, color=method, shape=method)) +
  geom_line(aes(x=time, y=se, color=method, group=method)) + facet_wrap(~state) +
  labs(title="Empirical Variability of State Occupation Probability Estimates",
       x="Time", y="Emprical Standard Error") +
  scale_color_brewer(type = "qual", palette = 6) +
  theme(text = element_text(size=10),
        axis.title = element_text(size=11),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size=13, hjust = 0.5))
ggsave(paste0("Output/Figures/Variability/variability_results_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".pdf"),
       height = 4, width = 6.5)


#Correlation Figure
#scatterplot of AJ vs Marginal MI estimates (only makes sense in settings where HR = 1)
corr_plot <- bias_results_plot %>%
  select(time, simulation, state, method, estimate) %>%
  pivot_wider(names_from = method, values_from = estimate)

ggplot(corr_plot) + geom_point(aes(x=`Aalen-Johansen`, y=`Marginal MI`, color=time), alpha=.8) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(~state) +
  labs(title="Association of Estimates Between Methods",
       x="Aalen-Johansen Estimate", y="Marginal MI Estimate") +
  theme(text = element_text(size=9),
        axis.title = element_text(size=11),
        plot.title = element_text(size=12, hjust = 0.5))
ggsave(paste0("Output/Figures/Correlation/correlation_results_k", as.character(s01), "_HR_",as.character(round(exp(b),2)), "_m01_" , as.character(m01), "_m02_", as.character(m02), "_m12_", as.character(m12),  ".pdf"),
       height = 4, width = 7)

#Notes:
#don't want to enforce proportional hazards for treatment vs control via simulation design
