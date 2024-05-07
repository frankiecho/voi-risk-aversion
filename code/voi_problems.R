library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)
library(reshape2)
library(patchwork)
library(future)

source("code/voi_simulations.R")

pref_list_long <- c("CE"="Certainty Equivalent", "MV"="Mean-Variance", "MCVaR"="Mean-Expected Shortfall")
pref_list <- c("CE", "MV", "MCVaR")

e1_ce_plt_list <- list()
e2_ce_plt_list <- list()
e1_vpi_plt_list <- list()
e2_vpi_plt_list <- list()

for (pref in pref_list[1]) {
  
  ## Example 1: 2 actions and two states
  if (pref=='CE') {
    gamma_seq <- seq(0.05, 5,.05)
  } else {
    gamma_seq <- seq(-1,1,.01)
  }
  e1 <- list()
  n_states <- 2
  n_actions <- 2
  a1 <- c(1.35, .55)
  a2 <- c(1.00, 1.00)
  e1$p <- c(0.5, 0.5)
  action_state <- matrix(c(a1, a2), byrow = T, ncol = length(a1))
  colnames(action_state) <- paste0('s', 1:n_states)
  rownames(action_state) <- paste0('a', 1:n_actions)
  e1$action_state <- action_state
  e1_voi <- fcn_VOI_simulation(e1, pref = pref, gamma_seq = gamma_seq) 
  
  e1_vpi_plt <- fcn_plt_vpi(e1_voi, n_actions)
  
  e1_ce <- lapply(gamma_seq, function(gamma) apply(e1$action_state, 1, pref_define(pref), lambda=gamma, p = e1$p)) %>%
    bind_rows()
  e1_ce$gamma <- gamma_seq
  e1_ce_plt <- fcn_plt_ce(e1_ce, n_actions)
  ggsave(e1_ce_plt / e1_vpi_plt, filename = paste0("plots/e1_obj_value_", pref, ".png"), height = 8, width = 6)
  
  ## Example 1: 3 actions and 3 states
  e2 <- list()
  n_states <- 3
  n_actions <- 3
  a1 <- c(0.689, 0.582, 0.547)
  a2 <- c(0.729, 0.674, 0.484)
  a3 <- c(0.745, 0.710, 0.332)
  e2$p <- c(0.4, 0.2, 0.4)
  action_state <- matrix(c(a1, a2, a3), byrow = T, ncol = length(a1))
  colnames(action_state) <- paste0('s', 1:n_states)
  rownames(action_state) <- paste0('a', 1:n_actions)
  e2$action_state <- action_state
  e2_voi <- fcn_VOI_simulation(e2, pref = pref, gamma_seq = gamma_seq)
  e2_vpi_plt <- fcn_plt_vpi(e2_voi, n_actions)
  
  e2_ce <- lapply(gamma_seq, function(g) apply(e2$action_state, 1, pref_define(pref), lambda=g, p = e2$p)) %>%
    bind_rows()
  e2_ce$gamma <- gamma_seq
  e2_ce_plt <- fcn_plt_ce(e2_ce, n_actions)
  
  ggsave(e2_ce_plt / e2_vpi_plt, filename = paste0("plots/e2_obj_value_", pref, ".png"), height = 6, width = 5)
  
  e1_ce_plt_list[[pref]] <- e1_ce_plt + ggtitle(pref_list_long[pref]) + theme(axis.line.x = element_blank(), axis.title.x = element_blank(),
                                                                axis.text.x = element_blank(), axis.ticks.x = element_blank())
  e1_vpi_plt_list[[pref]] <- e1_vpi_plt + ggtitle(pref_list_long[pref])
  e2_ce_plt_list[[pref]] <- e2_ce_plt + ggtitle(pref_list_long[pref]) + theme(axis.line.x = element_blank(), axis.title.x = element_blank(),
                                                              axis.text.x = element_blank(), axis.ticks.x = element_blank())
  e2_vpi_plt_list[[pref]] <- e2_vpi_plt + ggtitle(pref_list_long[pref])
}

e1_plt <- wrap_plots(e1_ce_plt_list) / wrap_plots(e1_vpi_plt_list)+ plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'a')
e2_plt <- wrap_plots(e2_ce_plt_list) / wrap_plots(e2_vpi_plt_list)+ plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'a')
ggsave(e1_plt, filename='plots/e1_plt.png', width = 12, height = 7)
ggsave(e2_plt, filename='plots/e2_plt.png', width = 12, height = 7)

## Several states
set.seed(1111)
n_states <- 100
n_actions <- 100

for (pref in pref_list) {
  
  ## Draw from uniform distribution like Holden et al (2024) ------
  plan(multisession)
  unif_action_state_sim <- function() matrix(runif(n_states*n_actions), nrow = n_actions)
  unif_plt <- fcn_plot_simulations(unif_action_state_sim, pref = pref)
  print('Completed unif_plt')
  
  ## Draw from an exponential distribution like Holden et al (2024) ------
  exp_action_state_sim <- function() matrix(rexp(n_states*n_actions, 1), nrow = n_actions)
  exp_plt <- fcn_plot_simulations(exp_action_state_sim, pref = pref)
  print('Completed exp_plt')
  
  ## Draw from an poisson distribution with heterogeneous rate parameters ------
  #pois_action_state_sim <- function() rpois(n_states*n_actions, rep(rbinom(n_actions, 30, 0.1), each=n_states)) %>%
  #  matrix(nrow = n_actions, byrow = T)
  #pois_plt <- fcn_plot_simulations(pois_action_state_sim, pref = pref)
  #print('Completed pois_plt')
  
  ## Draw from an lognormal distribution like Holden et al (2024) ------
  lnorm_action_state_sim <- function() rlnorm(n_states*n_actions, rep(0.1, each=n_states), 0.8) %>%
    matrix(nrow = n_actions, byrow = T)
  lnorm_plt <- fcn_plot_simulations(lnorm_action_state_sim, pref = pref)
  print('Completed lnorm_plt')
  
  ## Draw from a negative lognormal distribution ------
  neg_lnorm_action_state_sim <- function() -rlnorm(n_states*n_actions, rep(0.1, each=n_states), 0.8) %>%
    matrix(nrow = n_actions, byrow = T)
  neg_lnorm_plt <- fcn_plot_simulations(neg_lnorm_action_state_sim, pref = pref)
  print('Completed neg_lnorm')
  
  ## Draw from an kurtotic distribution ------
  t_action_state_sim <- function() rt(n_states*n_actions, 3) %>%
    matrix(nrow = n_actions, byrow = T)
  t_plt <- fcn_plot_simulations(t_action_state_sim, pref = pref)
  print('Completed t_plt')
  
  ## No trade-off in mean and variance
  mu <- seq(0.5, 1.5, length.out = n_actions)
  sd <- seq(0.5, 0.5, length.out = n_actions)
  mv_no_tradeoff_action_state_sim <- function() rnorm(n_states*n_actions, 
                                          rep(mu, each=n_states),
                                          rep(sd, each=n_states)) %>%
    matrix(nrow = n_actions, byrow = T)
  mv_no_tradeoff_plt <- fcn_plot_simulations(mv_no_tradeoff_action_state_sim, pref = pref)
  
  ## Explicit trade-off in mean and variance
  mu <- seq(0.5, 1.5, length.out = n_actions)
  sd <- seq(0.1, 3, length.out = n_actions)
  mv_action_state_sim <- function() rnorm(n_states*n_actions, 
                                            rep(mu, each=n_states),
                                            rep(sd, each=n_states)) %>%
    matrix(nrow = n_actions, byrow = T)
  mv_plt <- fcn_plot_simulations(mv_action_state_sim, pref = pref)
  
  ## Write plots -----
  fig1 <- unif_plt$order_plt + ggtitle("Uniform") + unif_plt$v_plt + unif_plt$plt + 
    exp_plt$order_plt + ggtitle("Exponential") + exp_plt$v_plt + exp_plt$plt + 
    plot_layout(guides='collect', byrow = F, nrow = 3) & theme(legend.position = "bottom") &
    plot_annotation(tag_levels = 'a')
  fig1
  ggsave(fig1, filename = paste0("plots/homogeneous_dist_", pref, ".png"), width = 10, height = 11, dpi = 300)
  
  fig2 <- lnorm_plt$order_plt + ggtitle("Lognormal") + lnorm_plt$v_plt + lnorm_plt$plt + 
    neg_lnorm_plt$order_plt + ggtitle("Lognormal (negative values)") + neg_lnorm_plt$v_plt + neg_lnorm_plt$plt + 
    t_plt$order_plt + ggtitle("T-distribution") + t_plt$v_plt + t_plt$plt +
    plot_layout(guides='collect',byrow = F, nrow = 3) & theme(legend.position = "bottom") &
    plot_annotation(tag_levels = 'a')
  fig2
  ggsave(fig2, filename = paste0("plots/lnorm_dist_", pref, ".png"), width = 12, height = 11, dpi = 300)
  
  fig3 <- #pois_plt$order_plt + ggtitle("Poisson") + pois_plt$v_plt + pois_plt$plt +
    mv_no_tradeoff_plt$order_plt + ggtitle("Normal distribution") + mv_no_tradeoff_plt$v_plt + mv_no_tradeoff_plt$plt +
    mv_plt$order_plt + ggtitle("Mean-Variance tradeoff") + mv_plt$v_plt + mv_plt$plt + 
    plot_layout(guides='collect',byrow = F, nrow = 3) & theme(legend.position = "bottom") &
    plot_annotation(tag_levels = 'a')
  fig3
  ggsave(fig3, filename = paste0("plots/mv_dist_", pref, ".png"), width = 10, height = 11, dpi = 300)
}
## Cleanup --------
plan(sequential)
