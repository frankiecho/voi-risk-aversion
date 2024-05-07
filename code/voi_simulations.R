# VOI under risk aversion
library(tidyverse)
library(mvtnorm)
library(ggpubr)
library(future)
library(future.apply)
library(Rmpfr)
library(matrixStats)
library(DescTools)
library(modi)

okabe_ito_colors = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2",   "#999999")

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
    res <- c
  } else {
    res <- (1-exp(-alpha*c))/alpha
  }
  return(res)
}

CARA_inv <- function(c, alpha = 0) {
  if (alpha == 0 | is.na(alpha)) {
    res <- c
  } else {
    res <- -(log(1-alpha*(c)))/alpha
  }
  return(res)
}

exp_u <- function(c, z = 0) {
  return(c^(1/z))
}

inv_exp_u <- function(u, z = 0) {
  return(u^(z))
}

fcn_is_real <- function(x) {
  !is.infinite(x) & !is.null(x) & !is.na(x)
}

# Function to calculate certainty equivalents of CARA function
CARA_CE_pref <- function(x, p=rep(1,length(x))/length(x), lambda = 0, mpfrPrec = NULL) {
  gamma = lambda
  if (is.null(mpfrPrec)) {
    ce <- CARA_inv(weighted.mean(CARA(x, gamma), p), gamma)
    if (!fcn_is_real(ce)) {
      # Retry with higher numerical precision with recursion until the computation returns a non-infinite result
      mpfrPrec <- 2000
      max_mpfrPrec <- 50000
      while (!fcn_is_real(ce) & mpfrPrec <= max_mpfrPrec) {
        ce <- CARA_CE_pref(x, p, lambda, mpfrPrec = mpfrPrec)
        mpfrPrec <- mpfrPrec + 2000
      }
      if (!fcn_is_real(ce)) {
        stop("CE is NA")
      }
    }
    return(ce)
  } else {
    x_mpfr <- mpfrArray(x, mpfrPrec)
    ce_mpfr <- CARA_inv(weighted.mean(CARA(x_mpfr, gamma), p), gamma)
    ce <- as.numeric(ce_mpfr)
    return(ce)
  }
}

CRRA_CE_pref <- function(x, p=rep(1,length(x))/length(x), lambda = 0, mpfrPrec = NULL) {
  gamma = lambda
  if (is.null(mpfrPrec)) {
    ce <- inv_exp_u(weighted.mean(exp_u(x, gamma), p), gamma)
    if (!fcn_is_real(ce)) {
      # Retry with higher numerical precision with recursion until the computation returns a non-infinite result
      mpfrPrec <- 2000
      max_mpfrPrec <- 50000
      while (!fcn_is_real(ce) & mpfrPrec <= max_mpfrPrec) {
        ce <- CARA_CE_pref(x, p, lambda, mpfrPrec = mpfrPrec)
        mpfrPrec <- mpfrPrec + 2000
      }
      if (!fcn_is_real(ce)) {
        stop("CE is NA")
      }
    }
    return(ce)
  } else {
    x_mpfr <- mpfrArray(x, mpfrPrec)
    ce_mpfr <- inv_exp_u(weighted.mean(exp_u(x_mpfr, gamma), p), gamma)
    ce <- as.numeric(ce_mpfr)
    return(ce)
  }
}

pref_define <- function(pref = 'CE') {
  if (pref == 'CE') {
    pref_func <- CRRA_CE_pref
  } else if (pref == 'MV'){
    pref_func <- MV
  } else if (pref == 'MCVaR') {
    pref_func <- MCVaR
  } else {
    stop(paste("Preference function not defined for", pref))
  }
  pref_func
}

# Mean-variance utility function given x, probability p, and lambda weighting parameter 
MV <- function(x, p=rep(1,length(x))/length(x), lambda=0) {
  if (lambda > 1 | lambda < -1) {
    stop("Lambda is only defined for -1 to 1")
  }
  return((1-abs(lambda))*weighted.mean(x, p) - lambda*weighted.var(x, p))
}

MCVaR <- function(x, p=rep(1,length(x))/length(x), lambda=0, beta = 0.1) {
  if (lambda > 1 | lambda < -1) {
    stop("Lambda is only defined for -1 to 1")
  }
  if(lambda >= 0) {
    VaR <- Quantile(x, weights = p*length(p), probs = beta)
    CVaR <- weighted.mean(x[x <= VaR], p[x <= VaR])
  } else {
    VaR <- Quantile(x, weights = p*length(p), probs = 1-beta)
    CVaR <- weighted.mean(x[x >= VaR], p[x >= VaR])
  }
  return((1-abs(lambda))*weighted.mean(x, p*length(p)) + abs(lambda)*CVaR)
}

# Expected Value of perfect information (compared in certainty equivalents)
fcn_EVPI <- function(gamma, action_state, p, pref = 'MCVaR') {
  # U <- function(c) CARA(c, gamma)
  # U_inv <- function(c) CARA_inv(c, gamma)
  # utility_table <- U(action_state)
  # max_a <- apply(utility_table, 2, which.max) %>% unlist() # Actions that maximise utility in each state (indices)
  # max_EU <- which.max(utility_table %*% p) # Action that maximises expected utility (index)
  # VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  # V_certainty <- U_inv(utility_table[cbind(max_a, seq_along(max_a))] %*% p)
  # V_uncertainty <- U_inv(utility_table[max_EU,] %*% p)
  
  
  pref_func <- pref_define(pref)
  utility_table <- action_state
  max_a <- apply(utility_table, 2, which.max) %>% unlist() # Actions that maximise utility in each state (indices)
  max_EU <- which.max(apply(utility_table, 1, pref_func, p = p, lambda = gamma))[[1]] # Action that maximises mean-variance utility
  VPI_value <- action_state[cbind(max_a, seq_along(max_a))] %*% p - action_state[max_EU,] %*% p
  V_certainty <- pref_func(action_state[cbind(max_a, seq_along(max_a))], p, gamma)
  V_uncertainty <- pref_func(action_state[max_EU,], p, gamma)

  VPI <- V_certainty - V_uncertainty
  if (any(is.na(VPI))) {
    warning('Some EVPI values are NA')
  }
  # Calculate the performance of the optimal action in the presence of uncertainty
  action_worst_outcomes <- apply(action_state, 1, min)
  max_EU_worst <- 1 + n_actions - match(action_worst_outcomes[max_EU], sort(action_worst_outcomes))
  
  action_best_outcomes <- apply(action_state, 1, max)
  max_EU_best <- 1 + n_actions - match(action_best_outcomes[max_EU], sort(action_best_outcomes))

  action_mean_outcomes <- action_state %*% p
  max_EU_mean <- 1 + n_actions - match(action_mean_outcomes[max_EU], sort(action_mean_outcomes))
  
  # Proportion where the optimal action under certainty and under uncertainty is the same
  p_same_action <- (max_EU == max_a) %*% p
  
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU, 
       V_certainty = V_certainty,
       V_uncertainty = V_uncertainty,
       action_worst_outcomes = max_EU_worst,
       action_mean_outcomes = max_EU_mean,
       action_best_outcomes = max_EU_best,
       p_same_action = p_same_action)
}


# Expected Value of Partial Perfect Information (compared in certainty equivalents)
# gamma: Risk aversion coefficient
# action_state: action state "payoff" matrix
# p: probability of each state (before resolving uncertainty)
# Y: a list of probable outcomes of the experiment, in a list of lists. In each list, the index of the states that are possible conditional on the experiment is given. 
#   For example, a three-state system with the first list of Y being c(1,2) and second being c(3) means that the experiment can have two outcomes - (1) where state 3 is ruled out
#   but it is still uncertain over state 1 or 2, or (2) where state 3 is confirmed
# pY: probability of the outcome(s) of the experiments in Y
fcn_EVSI <- function(gamma, action_state, p=rep(1,ncol(action_state))/ncol(action_state), Y, pY=length(Y), pref = 'CE') {
  pref_func <- pref_define(pref)
  utility_table <- apply(action_state, 1, pref_func, p = p, lambda = gamma) # utility of each action under complete uncertainty
  max_EU <- which.max(utility_table)[[1]] # Action that maximises utility under uncertainty
  
  n_y <- length(Y) # Number of possible outcomes of the experiment
  utility_table_Y <- lapply(Y, function(y) {
    if (length(y)==1) return(sapply(action_state[,y], \(x) pref_func(x, p = p[y]/sum(p[y]), lambda = gamma)))
    apply(action_state[,y], 1, pref_func, p = p[y]/sum(p[y]), lambda = gamma) 
  }) %>% bind_rows() # The certainty equivalent of each outcome in Y, after resolving uncertainty in experiments
  max_a <- apply(utility_table_Y, 1, which.max) # Which action maximises utility for each experimental outcome
  
  VPI_value <- sum(sapply(1:n_y, \(y) mean(action_state[cbind(max_a[y], Y[[y]])])) * pY)
  outcomes_certainty <- do.call(c, lapply(1:n_y, \(y) action_state[cbind(max_a[y], Y[[y]])]))
  p_certainty <- do.call(c, lapply(1:n_y, \(y) pY[y] * p[Y[[y]]]/sum(Y[[y]])))
  
  V_certainty <- pref_func(outcomes_certainty, p_certainty, gamma)
  V_uncertainty <- pref_func(action_state[max_EU,], p, gamma)
  
  VPI <- V_certainty - V_uncertainty
  if (any(is.na(VPI))) {
    warning('Some EVPI values are NA')
  }
  # Calculate the performance of the optimal action in the presence of uncertainty
  action_worst_outcomes <- apply(action_state, 1, min)
  max_EU_worst <- 1 + n_actions - match(action_worst_outcomes[max_EU], sort(action_worst_outcomes))
  
  action_best_outcomes <- apply(action_state, 1, max)
  max_EU_best <- 1 + n_actions - match(action_best_outcomes[max_EU], sort(action_best_outcomes))
  
  action_mean_outcomes <- action_state %*% p
  max_EU_mean <- 1 + n_actions - match(action_mean_outcomes[max_EU], sort(action_mean_outcomes))
  
  # Proportion where the optimal action under certainty and under uncertainty is the same
  p_same_action <- (max_EU == max_a) %*% pY
  
  list(VPI = VPI, VPI_value = VPI_value, max_EU = max_EU, 
       V_certainty = V_certainty,
       V_uncertainty = V_uncertainty,
       action_worst_outcomes = max_EU_worst,
       action_mean_outcomes = max_EU_mean,
       action_best_outcomes = max_EU_best,
       p_same_action = p_same_action)
}

fcn_VOI_simulation  <- function(voi_problem, gamma_seq = NULL, pref = "CE") {
  
  if (is.null(gamma_seq)) {
    if (pref=='CE') {
      gamma_seq <- seq(-5,5,.05)
    } else {
      gamma_seq <- seq(-1,1,.01)
    }
  }
  
  if(is.null(voi_problem$evsi)) {
    voi_problem$evsi <- FALSE
  }
  if (!voi_problem$evsi) {
    voi_result  <- lapply(gamma_seq, fcn_EVPI, action_state = voi_problem$action_state, p = voi_problem$p, pref = pref) 
  } else {
    voi_result  <- lapply(gamma_seq, fcn_EVSI, action_state = voi_problem$action_state, 
                       p = voi_problem$p, Y = voi_problem$Y, pY = voi_problem$pY) 
  }
  evpi_seq <- voi_result %>%
    lapply(function(x) c(VPI = x$VPI, VPI_value = x$VPI_value, max_EU = x$max_EU, 
                         action_best_outcomes = x$action_best_outcomes,
                         action_mean_outcomes = x$action_mean_outcomes,
                         action_worst_outcomes = x$action_worst_outcomes,
                         V_certainty = x$V_certainty,
                         V_uncertainty = x$V_uncertainty,
                         p_same = x$p_same)) %>%
    bind_rows(.id = 'gamma')

  evpi_seq$gamma <- gamma_seq
  return(evpi_seq)
}

# Simulation function for any distribution
# action_state: matrix of payoffs for each action in each state
fcn_voi_simulation_distribution <- function(action_state_sim_func, n_y = 100, pref = 'CE', nsims = 500) {
  fcn_sim <- function() {
    prob <- list()
    prob$action_state <- action_state_sim_func()
    n_states <- ncol(prob$action_state)
    n_actions <- nrow(prob$action_state)
    colnames(prob$action_state) <- paste0('s', 1:n_states)
    rownames(prob$action_state) <- paste0('a', 1:n_actions)
    prob$p <- rep(1/n_states, n_states)
    prob$Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])
    prob$evpxi <- n_y < n_states
    prob
  }
  sim <- future_replicate(nsims, fcn_sim(), simplify = F)
  sim_results <- future_lapply(sim, fcn_VOI_simulation, pref = pref)
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
      geom_ribbon(aes(ymin = lb, ymax = ub), fill = 'gray60')
  }
  
  if (is.null(col)) {
    plt +
      geom_line(color = 'gray20', linewidth = 1) +
      theme_pubr()
  } else {
    plt +
      geom_line(aes_string(color = col), size = 1) +
      theme_pubr()
  }
}

fcn_plot_simulations <- function(action_state, pref = 'CE') {
  order_var <- c('action_worst_outcomes','action_mean_outcomes','action_best_outcomes')
  values_cert_uncert <- c("V_certainty", "V_uncertainty")
  n_actions <- nrow(action_state())
  action_sim <- fcn_voi_simulation_distribution(action_state, pref = pref)
  action_sim_table <- fcn_summarise_table(action_sim)
  action_sim_order <- lapply(order_var, \(x) fcn_summarise_table(action_sim, variable = x))
  names(action_sim_order) <- order_var
  action_sim_order <- action_sim_order %>% bind_rows(.id = 'name')
  plt <- fcn_plt_sim_table(action_sim_table, additional_gg = list(geom_vline(xintercept = 0, color = 'gray50'))) +
    scale_y_continuous("Risk-adjusted\n Value of Perfect Information") +
    scale_x_continuous("Risk Aversion Coefficient") +
    theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank())
  
  v_certainty <- lapply(values_cert_uncert, \(x) fcn_summarise_table(action_sim, variable = x))
  names(v_certainty) <- values_cert_uncert
  v_certainty <- v_certainty %>% bind_rows(.id = 'name')
  v_plt <- v_certainty %>%
    mutate(name = factor(name, levels = values_cert_uncert, labels = c('Value under certainty', 'Value under uncertainty'))) %>%
    fcn_plt_sim_table(F, 'name', additional_gg = list(geom_vline(xintercept = 0, color = 'gray50'))) +
    guides(color = guide_legend(title="")) +
    scale_color_manual(values = okabe_ito_colors[4:5])+
    scale_y_continuous("Risk-adjusted Value") +
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank())

  order_plt <- action_sim_order %>%
    mutate(name = factor(name, levels = order_var, labels = c('Minimum', 'Mean', 'Maximum'))) %>%
    fcn_plt_sim_table(F, 'name', additional_gg = list(geom_hline(yintercept = 1, color = 'gray50'),
                                                      geom_vline(xintercept = 0, color = 'gray50'))) +
    scale_y_reverse("Rank of action \nchosen under uncertainty") +
    guides(color = guide_legend(title="Order payoffs by")) +
    scale_color_manual(values = okabe_ito_colors[1:3])+
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank()) +
    coord_cartesian(ylim = c(n_actions, -n_actions*0.1)) +
    annotate("text", x = ifelse(pref=='CE',-2.5, -0.5), y = -n_actions*0.05, label = "Risk-loving", vjust = 0.5, hjust = 0.5)+
    annotate("text", x = ifelse(pref=='CE', 2.5, 0.5), y = -n_actions*0.05, label = "Risk-averse", vjust = 0.5, hjust = 0.5)
  
  p_same_sim_table <- fcn_summarise_table(action_sim, variable = 'p_same')
  p_same_plt <- fcn_plt_sim_table(p_same_sim_table, additional_gg = list(geom_vline(xintercept = 0, color = 'gray50'))) +
    scale_y_continuous("Probability information \ndoes not change action") +
    scale_x_continuous("Risk Aversion Coefficient") +
    theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank())
    
  list(order_plt = order_plt, v_plt = v_plt, plt = plt, p_same_plt = p_same_plt)
}

fcn_plot_corr_mat <- function(sigma_list, limits = c(NA, NA)) {
  names(sigma_list) <- 1:length(sigma_list)
  sigma_df <- lapply(sigma_list, function(sigma) {
    colnames(sigma) <- rownames(sigma) <- 1:nrow(sigma)
    melted_cormat <- reshape2::melt(sigma)
  }) %>%
    bind_rows(.id = 'name')
  ggplot(sigma_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_y_reverse() +
    scale_fill_gradient2(midpoint = 0, limits = limits) +
    facet_wrap(~name) +
    theme_minimal()
}

fcn_plt_ce <- function(df, n_actions = 2) {
  df <- df %>%
    pivot_longer(paste0('a', 1:n_actions), names_to = 'actions', values_to = 'ce') %>%
    mutate(actions = factor(actions, paste0('a', 1:n_actions), as.character(1:n_actions))) 
  df %>%
    ggplot(aes(y = ce, x = gamma, color = actions))+
    geom_hline(yintercept = 0, color = 'gray50') +
    geom_vline(xintercept = 0, color = 'gray50') +
    geom_line(linewidth = 1)+
    geom_point(data = filter(df, gamma==0), size = 2) +
    scale_x_continuous("Risk Aversion Coefficient") +
    scale_y_continuous("Value of system under uncertainty") +
    scale_color_manual("Action", values=okabe_ito_colors[c(1,3,5)[1:n_actions]], drop=F) +
    theme_pubr() +
    theme(legend.position = 'right', 
          panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank())
}

fcn_plt_vpi <- function(df, n_actions = 2) {
  df <- df %>%
    mutate(actions = factor(max_EU, 1:n_actions, as.character(1:n_actions))) %>%
    select(gamma, VPI, max_EU, actions) 
  
  df %>%
    ggplot(aes(x = gamma, y = VPI)) +
    geom_hline(yintercept = 0, color = 'gray50') +
    geom_vline(xintercept = 0, color = 'gray50') +
    geom_line(color = 'gray50', linewidth = 1) +
    geom_line(aes(color = actions), linewidth = 1) +
    geom_point(data = filter(df, gamma==0), aes(color = actions), size = 2) +
    scale_x_log10("Risk Aversion Coefficient") +
    scale_y_continuous("Value of Perfect Information") +
    scale_color_manual("Optimal \nAction \nunder \nUncertainty", values=okabe_ito_colors[c(1,3,5)[1:n_actions]], drop=F) +
    theme_pubr() +
    theme(legend.position = 'right', 
          panel.border = element_rect(linewidth = 0.5, fill = NA),
          axis.line = element_blank())
}

