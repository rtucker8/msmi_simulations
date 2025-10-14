#Purpose: Characterize different simulation settings
#Author: Rachel Gonzalez
#Date: September 5 2025

#Libraries
library(tidyverse)
library(survival)
library(survminer)
library(kableExtra)
library(mstate)
library(mici)
library(philentropy) #calculate KL divergence
library(patchwork)

setwd("/Users/rachelgonzalez/Documents/Dissertation/msmi_simulations")
# Functions ---------------------------------------------------------------
source("simulation_helper_V2.R")

# Simulation Settings -----------------------------------------------------
b<- c(rep(0,4), rep(log(0.5), 4))

m01 <- c(rep(c(3,3,8,8), 2))
m02 <- c(rep(c(8,8,3,3), 2))
m12 <- c(rep(c(3,6), 4))

settings <- data.frame(b, m01, m02, m12, s01 =1.3, s02 = 1.5, s12 = 2) %>%
  bind_rows(data.frame(b, m01, m02, m12, s01 =0.7, s02 = 1.5, s12 = 2)) %>%
  mutate(setting = 1:nrow(.), theta = 12) %>%
  select(setting, everything())

settings <- bind_rows(settings, c(setting = 17, b=0, m01=8, m02=1, m12=3, s01=1.3, s02=1.5, s12=2, theta=12))

# Operating Characteristics ------------------------------------------------

#Simulate a large dataset under each setting
set.seed(2027)

dsets <- pmap(settings, function(b, m01, m02, m12, s01, s02, s12, theta, ...) {
 simulate_illness_death(
    n= 8000000,
    beta = b,
    median01 = m01,
    median02 = m02,
    median12 = m12,
    shape01 = s01,
    shape02 = s02,
    shape12 = s12,
    theta = theta,
    return_latent_data = TRUE,
    return_censored_data=FALSE
  )
})

#Calculate proportion of individuals that experience illness
settings$ever_ill <- map_dbl(dsets, ~ mean(.x[, "event1"] == 1))

# #Plot the density of sojourn times with violin plots
# violin <- function(d) {
#   d <- d %>%
#     select(sojourn01, sojourn02, sojourn12) %>%
#     pivot_longer(cols = everything(), names_to = "transition", names_prefix = "sojourn", values_to = "time") %>%
#     filter(!is.na(time))
#
#   ggplot(data=d, mapping = aes(x=transition, y=time, fill=transition)) +
#     geom_violin()
#
# }



#Calculate the true state occupation probabilities and plot over time
empirical_truth <- map(dsets, function(x) {
  get_empirical_probs(x, times = seq(0, 12, by = 0.25))  %>%
    select(time, pHealthy, pIll, pDead) %>%
    pivot_longer(cols = c("pHealthy", "pIll", "pDead"), names_to = "state", names_prefix = "p", values_to = "probability")
})

empirical_truth_long <- bind_rows(empirical_truth, .id = "setting") %>% mutate(setting = as.numeric(setting),
                                                                 state = factor(state, levels = c("Healthy", "Ill", "Dead"), ordered=TRUE))

ggplot(empirical_truth_long) +
  geom_line(aes(x=time, y=probability, color=state, linetype = state), linewidth = 1.2) +
  ggforce::facet_wrap_paginate(~setting, ncol = 4, nrow = 4, page = 1) +
  labs(title="", x="Time", y="State Occupation Probability") +
  scale_color_brewer(type="qual", palette = 1) +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,12,2), limits=c(0,12)) +
  theme(text = element_text(size=14),
    axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="grey90"))
ggsave("model_dynamics.pdf", width=10, height=8)


save(empirical_truth, settings, file = "Truth/sim_setting_characteristics.RData")


#Extra Code:

# #verify assumption
# # Filter datset to those who became ill
# ill <- subset(d_large_n$d.observed, event1 == 1)
# # Cox model: include time of illness onset as covariate
# markov_test <- coxph(Surv(t2-t1, event2) ~ t1, data = ill)
# summary(markov_test)

# Compare Empirical Estimates with Truth -------------------------------------------------------------------
# #(semi-markov case only)
# truth = get_truth(shape01 = s01, median01 = m01,
#           shape02 = s02, median02 = m02,
#           shape12 = s12, median12 = m12, t=eval_times)
# check_truth = empirical_truth %>% left_join(truth) %>%
#   select(time, pHealthy, p0.star, pIll, p1.star, pDead, p2.star, pDeadMinusIll, p2.0.star, pDeadWithIll, p2.1.star)
