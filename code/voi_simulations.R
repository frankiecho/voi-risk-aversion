# VOI under risk aversion
library(tidyverse)
library(mvtnorm)
library(ggpubr)
library(future)
library(future.apply)
library(Rmpfr)

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
    res <- -(log(1-alpha*c))/alpha
    if (is.na(res)) {
      return(res)
    }
    return(res)
  }
}

# Expected Value of perfect information (compared in certainty equivalents)
fcn_EVPI <- function(gamma, action_state, p) {
  U <- function(c) CARA(c, gamma)
  U_inv <- function(c) CARA_inv(c, gamma)
  utility_table <- U(action_state)
  max_a <- apply(utility_table, 2, which.max) %>% unlist() # Actions that maximise utility in each state (indices)
  max_EU <- which.max(utility_table %*% p) # Action that maximises expected utility (index)
  VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  VPI <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p) - U_inv(utility_table[max_EU,] %*% p)
  if (any(is.na(VPI))) {
    warning('Some EVPI values are NA')
  }
  # Calculate the performance of the optimal action in the presence of uncertainty
  action_worst_outcomes <- apply(action_state, 1, min)
  max_EU_worst <- 1 + n_actions - match(action_worst_outcomes[max_EU], sort(action_worst_outcomes))
  
  action_best_outcomes <- apply(action_state, 1, max)
  max_EU_best <- 1 + n_actions - match(action_best_outcomes[max_EU], sort(action_best_outcomes))
  
  action_mean_outcomes <- apply(action_state, 1, mean)
  max_EU_mean <- 1 + n_actions - match(action_mean_outcomes[max_EU], sort(action_mean_outcomes))
  
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU, 
       action_worst_outcomes = max_EU_worst,
       action_mean_outcomes = max_EU_mean,
       action_best_outcomes = max_EU_best)
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
  VPI <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p) - U_inv( utility_table[max_EU,] %*% p)
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
    lapply(function(x) c(VPI = x$VPI, VPI_value = x$VPI_value, max_EU = x$max_EU, 
                         action_best_outcomes = x$action_best_outcomes,
                         action_mean_outcomes = x$action_mean_outcomes,
                         action_worst_outcomes = x$action_worst_outcomes)) %>%
    bind_rows(.id = 'gamma')

  evpi_seq$gamma <- gamma_seq
  return(evpi_seq)
}

# Simulation function for any distribution
# action_state: matrix of payoffs for each action in each state
fcn_voi_simulation_distribution <- function(action_state_sim_func, n_y = 100, gamma_seq = seq(-5,5,0.1)) {
  fcn_sim <- function() {
    prob <- list()
    prob$action_state <- action_state_sim_func()
    n_states <- ncol(prob$action_state)
    n_actions <- nrow(prob$action_state)
    colnames(prob$action_state) <- paste0('s', 1:n_states)
    rownames(prob$action_state) <- paste0('a', 1:n_actions)
    p <- runif(n_states)
    prob$p <- p/sum(p)
    prob$Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])
    prob$evpxi <- n_y < n_states
    prob
  }
  sim <- replicate(nsims, fcn_sim(), simplify = F)
  sim_results <- future_lapply(sim, fcn_VOI_simulation, gamma_seq = gamma_seq)
  sim_table <- sim_results %>%
    #lapply(function(x) mutate(x, VPI = VPI - x$VPI[x$gamma == 0])) %>%
    #lapply(function(x) x$VPI) %>%
    bind_rows(.id = 'run_index')
  return(sim_table)
}

fcn_summarise_table <- function(sim_table,variable = "VPI") {
  out <- sim_table %>%
    mutate(var := !!as.name(variable)) %>%
    filter(is.infinite(var) == F) %>%
    group_by(gamma) %>%
    summarise(median := median(var, na.rm = T), lb := quantile(var, 0.05, na.rm = T), 
              ub := quantile(var, 0.95, na.rm = T), llb := quantile(var, 0.01, na.rm = T), 
              uub := quantile(var, 0.99, na.rm = T)) 
}

fcn_plt_sim_table <- function(table, ribbon = TRUE, col = NULL, additional_gg = NULL) {
  plt <- ggplot(table, aes(x = gamma, y = median))
  
  if (!is.null(additional_gg)) {
    plt <- plt + additional_gg
  }
  
  if (ribbon) {
    plt <- plt +
      geom_ribbon(aes(ymin = llb, ymax = uub), fill = 'gray80') +
      geom_ribbon(aes(ymin = lb, ymax = ub), fill = 'gray70')
  }
  
  if (is.null(col)) {
    plt +
      geom_line(color = 'gray20', size = 1) +
      theme_pubr()
  } else {
    plt +
      geom_line(aes_string(color = col), size = 1) +
      theme_pubr()
  }
}

fcn_plot_simulations <- function(action_state) {
  order_var <- c('action_worst_outcomes','action_mean_outcomes','action_best_outcomes')
  action_sim <- fcn_voi_simulation_distribution(action_state)
  action_sim_table <- fcn_summarise_table(action_sim)
  action_sim_order <- lapply(order_var, \(x) fcn_summarise_table(action_sim, variable = x))
  names(action_sim_order) <- order_var
  action_sim_order <- action_sim_order %>% bind_rows(.id = 'name')
  plt <- fcn_plt_sim_table(action_sim_table, additional_gg = list(geom_vline(xintercept = 0, color = 'gray50'))) +
    scale_y_continuous("Certainty-equivalent\n Value of Perfect Information") +
    scale_x_continuous("Risk Aversion Coefficient")
  order_plt <- action_sim_order %>%
    mutate(name = factor(name, levels = order_var, labels = c('Minimum', 'Mean', 'Maximum'))) %>%
    fcn_plt_sim_table(F, 'name', additional_gg = list(geom_hline(yintercept = 1, color = 'gray50'),
                                                      geom_vline(xintercept = 0, color = 'gray50'))) +
    scale_y_reverse("Rank") +
    guides(color = guide_legend(title="Order payoffs by")) +
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    annotate("text", x = -2.5, y = 0, label = "Risk-loving", vjust = 0, hjust = 0.5)+
    annotate("text", x = 2.5, y = 0, label = "Risk-averse", vjust = 0, hjust = 0.5)
  list(order_plt = order_plt, plt = plt)
}
