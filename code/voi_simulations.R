# VOI under risk aversion
library(tidyverse)
library(mvtnorm)

# vNM Utility functions
CRRA <- function(c, gamma = 1) {
  if (gamma == 1) {
    return(log(c))
  } else if (gamma > 0) {
    return(c^(1-gamma) / (1-gamma))
  } else {
    return(NA)
  }
}

CRRA_inv <- function(c, gamma = 1) {
  if (gamma == 1) {
    return(exp(c))
  } else if (gamma > 0) {
    return((1-gamma)^(1/(1-gamma))*(c^(1/(1-gamma))))
  } else {
    return(NA)
  }
}

CARA <- function(c, alpha = 0) {
  if (alpha == 0 | is.na(alpha)) {
    return(c)
  } else {
    return((1-exp(-alpha*c))/alpha)
  }
}

CARA_inv <- function(c, alpha = 0) {
  if (alpha == 0 | is.na(alpha)) {
    return(c)
  } else {
    return(-(log(1-alpha*c))/alpha)
  }
}

# Expected Value of perfect information (compared in certainty equivalents)
fcn_EVPI <- function(gamma, action_state, p) {
  U <- function(c) CARA(c, gamma)
  U_inv <- function(c) CARA_inv(c, gamma)
  utility_table <- U(action_state)
  max_a <- apply(utility_table, 2, which.max)  # Actions that maximise utility in each state (indices)
  max_EU <- which.max(utility_table %*% p) # Action that maximises expected utility (index)
  VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  VPI <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p) - U_inv(utility_table[max_EU,] %*% p)
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU)
}

# Expected Value of Partial Perfect Information (compared in certainty equivalents)
fcn_EVPXI <- function(gamma, action_state, p, Y) {
  U <- function(c) CARA(c, gamma)
  U_inv <- function(c) CARA_inv(c, gamma)
  utility_table <- U(action_state)
  n_y <- max(unique(Y))
  max_ay <- lapply(1:n_y, function(y) which.max(utility_table[,Y == y]%*%(p[Y == y]/sum(p[Y == y])))) %>% unlist() # Actions that maximise utility in each experimental outcome
  max_EU <- which.max(utility_table %*% p) # Action that maximises expected utility (index) across all states
  max_a <- max_ay[Y]
  VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  VPI <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p) - U_inv(utility_table[max_EU,] %*% p)
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU)
}

fcn_VOI_simulation  <- function(voi_problem, gamma_seq) {
  if(is.null(voi_problem$evpxi)) {
    voi_problem$evpxi <- FALSE
  }
  if (!voi_problem$evpxi) {
    voi_result  <- lapply(gamma_seq, fcn_EVPI, action_state = voi_problem$action_state, p = voi_problem$p) 
  } else {
    voi_result  <- lapply(gamma_seq, fcn_EVPXI, action_state = voi_problem$action_state, 
                       p = voi_problem$p, Y = voi_problem$Y) 
  }
  evpi_seq <- voi_result %>%
    lapply(function(x) c(VPI = x$VPI, VPI_value = x$VPI_value, max_EU = x$max_EU)) %>%
    bind_rows(.id = 'gamma')
  evpi_seq$gamma <- gamma_seq
  return(evpi_seq)
}
