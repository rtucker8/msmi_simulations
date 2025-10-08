#Purpose: Build helper functions for multi-state model simulations (no redundant functions from msmi package)
#Author: Rachel Gonzalez
#Date: October 6 2025


# General Purpose Functions -----------------------------------------------

#Resample function that fixes behavior or sample function when n=1
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Return value of a step function at a given time
get_step_value <- function(step_times, step_values, t) {
  
  # Find the largest time in step_times that is less than or equal to t
  idx <- max(which(step_times <= t), na.rm = TRUE)
  
  # If no such time exists, return NA
  if (is.na(idx)) {
    return(NA)
  }
  
  # Return the corresponding value from step_values
  return(step_values[idx])
}

#Compute the KL Divergence and Hellinger distance between two discrete probability distributions
# rows of df_reference and df_estimate must correspond to the same time points as eval_times
compute_distance <- function(df_reference, df_estimate, 
                             cols_reference = c("pHealthy", "pIll", "pDead"),
                             cols_estimate = c("pstate1", "pstate2", "pstate3"), 
                             times, epsilon = 1e-12) {
  
  distance_df <- lapply(times, function(t){
    
    p_ref <- as.numeric(df_reference[t, cols_reference]) + epsilon
    p_est <- as.numeric(df_estimate[t, cols_estimate]) + epsilon
    
    probs <- rbind(p_est, p_ref)
    
    # Compute KL divergence: philentropy expects rows = distributions
    # and computes KL(row1 || row2)
    KL_result <- philentropy::KL(probs, unit = "log2")
    hellinger_result <- (1/sqrt(2))*sqrt(sum((sqrt(p_est) - sqrt(p_ref))^2))
    
    distance_df = tibble(
      time = t,
      KL= KL_result,
      hellinger = hellinger_result)
  })
  
  return(bind_rows(distance_df))
}
# Data Generation Process -------------------------------------------------

#Calculate scale parameters from Weibull distribution given median and shape
get_weibull_scale <- function(median, shape) {
  scale = log(2)/(median^shape);
  return(scale)
}

#Generate Weibull distributed survival times (no covariates)
generate_weibull <- function(n, shape, median, scale=NULL, seed=NULL) {
  
  set.seed(seed)
  
  # Generate uniform random variables
  u <- runif(n, 0, 1)
  
  if (is.null(scale)) {
    #Use Weibull shape and desired median survival time to get scale
    scale = get_weibull_scale(median, shape)
  }
  # Generate survival times from Weibull distribution using inverse CDF
  t <- (-log(u)/scale)^(1/shape)
  
  return(t)
}

#Generate Cox-Weibull distributed survival times 
generate_cox_weibull <- function(n, shape, median, X, beta) {
  # Generate uniform random variables
  u <- runif(n, 0, 1)
  #Use Weibull shape and desired median survival time to get scale
  scale = get_weibull_scale(median, shape)
  # Generate survival times from a Cox PH model with Weibull baseline hazard using inverse CDF
  t <- (-log(u)/(scale*exp(beta*X)))^(1/shape)
  return(t)
}

#Simulate from a multi-state model with three states: healthy (0), ill (1), dead (2)
#n- sample size (numeric)
#beta- effect (logHR) of T1 on lambda_12. When beta = 0, the MSM satisfies the Markov assumption. Default = 0. (numeric)
#shapeij- Weibull distribution shape for the i -> j transition (numeric > 0)
#medianij - median for the sojourn time between i -> j, used to calculate the Weibull scale for the i -> j transition (numeric > 0)
#return_latent_data- indicates if the function should return the underlying survival times in addition to the censored data that is observed. (logical)
#theta- upper bound for the uniform distribution for censoring (numeric)
simulate_illness_death <- function(n, beta = 0,
                                   shape01, median01,
                                   shape02, median02, 
                                   shape12, median12,
                                   return_latent_data = FALSE, 
                                   return_censored_data = TRUE,
                                   theta = NULL,
                                   seed=sample(1:.Machine$integer.max, 1)) {
  #Set up random seed for data generation
  if (is.null(seed)) {
    warning("Please provide a seed for reproducibility. A random seed has been generated.")
  }
  
  set.seed(seed)
  
  #Determine first event (illness or death)
  sojourn01 <- generate_weibull(n, shape01, median01) 
  sojourn02 <- generate_weibull(n, shape02, median02) 
  
  #Generate time from illness to death only for those with illness
  sojourn12 <- rep(NA, n)
  ill_ids <- which(sojourn01 < sojourn02)
  adjusted_scale_12 <- get_weibull_scale(median12, shape12) * exp(beta * sojourn01[ill_ids])
  sojourn12[ill_ids] <- generate_weibull(length(ill_ids), shape12, median = NULL, scale = adjusted_scale_12)
  
  temp = data.frame(sojourn01 = sojourn01, sojourn02 = sojourn02, sojourn12 = sojourn12)
  
  #Time of entry into each state and event indicators
  temp = temp %>% mutate(id = 1:n,
                         t1 = if_else(sojourn01 < sojourn02, sojourn01, sojourn02),
                         event1 = if_else(sojourn01 < sojourn02, 1, 0),
                         t2 = if_else(sojourn01 < sojourn02, sojourn01 + sojourn12, sojourn02),
                         event2 = 1) %>%
    select(id, t1, event1, t2, event2, sojourn01, sojourn02, sojourn12)
  
  #Add in uniform censoring times
  temp_censoring <- temp %>% mutate(C = runif(n, 0, theta),
                                    event1 = if_else(C < t1, 0, event1),
                                    t1 = if_else(C < t1, C, t1),
                                    event2 = if_else(C < t2, 0, event2),
                                    t2 = if_else(C < t2, C, t2))
  
  temp_censoring <- temp_censoring %>% select(id, t1, event1, t2, event2, sojourn12, sojourn01, sojourn02)
  
  if (return_latent_data == FALSE & return_censored_data == TRUE) {
    return(d.observed = temp_censoring)
  } else if (return_censored_data == FALSE & return_latent_data == TRUE) {
    return(d.latent = temp)
  } else {
    return(list(d.observed = data.frame(temp_censoring), d.latent = data.frame(temp)))
  }
  
}


# Empirical Probability in State ------------------------------------------

get_empirical_probs <- function(df, times) {
  
  # Initialize an empty list to store results for each time
  results <- lapply(times, function(t) {
    # logical conditions for each state
    is_healthy <- df$t1 > t & df$t2 > t
    is_ill <- df$t1 <= t & df$t2 > t
    is_dead_no_illness <- df$t2 <= t & df$event1 == 0
    is_dead_with_illness <- df$t2 <= t & df$event1 == 1
    
    # Calculate proportions in each state
    n <- nrow(df)
    pHealthy <- sum(is_healthy) / n
    pIll <- sum(is_ill) / n
    pDeadMinusIll <- sum(is_dead_no_illness) / n
    pDeadWithIll <- sum(is_dead_with_illness) / n
    pDead <- pDeadMinusIll + pDeadWithIll
    
    tibble(
      time = t,
      pHealthy = pHealthy,
      pIll = pIll,
      pDeadMinusIll = pDeadMinusIll,
      pDeadWithIll = pDeadWithIll,
      pDead = pDead
    )
  })
  
  bind_rows(results)
}

# True Probability in State -----------------------------------------------

#probability in healthy state
p0.truth <- function(shape01, scale01, shape02, scale02, t) {
  return(exp(-scale01*(t^shape01) - scale02*(t^shape02)))
}

#probability in ill state
p1.truth <- function(shape01, scale01, shape02, scale02, shape12, scale12, t){
  
  lambda01 <- function(u) scale01 * shape01 * (u^(shape01 - 1))
  
  compute_p1 <- function(t) {
    #inner integral: ∫ lambda23(u) du from r to t
    S2 <- function(r0) {
      exp(-scale12*((t-r0)^shape12))
    }
    #outer integral: ∫ lambda12(r) * p0(r) * S2(r) ds
    integrand <- function(r) {
      lambda01(r) * exp(-scale01*(r^shape01) - scale02*(r^shape02)) * S2(r)
    }
    
    p1 <- integrate(integrand, lower=0, upper=t)$value
  }
  
  p1.star <- sapply(t, compute_p1)
  
  return(p1.star)
}

#Probability of being dead at time t not having previously been ill
p2_0.truth<- function(shape01, scale01, shape02, scale02, shape12, scale12, t){
  
  compute_p2_0 <- function(t) {
    #probability of bing healthy untill time r: S_0(t) = exp(-scale01*t^shape01 - scale02*t^shape02))
    integrand <- function(r) {
      return(exp(-scale01*(r^shape01) - scale02*(r^shape02))*scale02*shape02*(r^(shape02-1)))
    }
    
    p2_0 <- integrate(integrand, lower=0, upper=t)$value
  }
  
  p2_0.star <- sapply(t, compute_p2_0)
  
  
  return(p2_0.star)
}

get_truth <- function(shape01, median01, shape02, median02, shape12, median12, t){
  
  #transform median to scale parameters
  scale01 = get_weibull_scale(median01, shape01)
  scale02 = get_weibull_scale(median02, shape02)
  scale12 = get_weibull_scale(median12, shape12)
  
  #times to evaluate true probabilities
  times = t
  
  #probability in healthy state
  p0.star <- p0.truth(shape01=shape01, scale01=scale01, shape02=shape02, scale02=scale02, t=times)
  
  #probability in ill state
  p1.star <- p1.truth(shape01=shape01, scale01=scale01, shape02=shape02, scale02=scale02, shape12 = shape12, scale12=scale12, t=times)
  
  #probability in dead state
  p2.0.star <- p2_0.truth(shape01=shape01, scale01=scale01, shape02=shape02, scale02=scale02, shape12 = shape12, scale12=scale12, t=times)
  p2.1.star <- 1- p0.star - p1.star - p2.0.star
  p2.star <- p2.0.star + p2.1.star
  
  probs.star <- data.frame(time = times,
                           p0.star=p0.star,
                           p1.star=p1.star, 
                           p2.0.star=p2.0.star,
                           p2.1.star=p2.1.star,
                           p2.star=p2.star)
  
  return(probs.star)
}

