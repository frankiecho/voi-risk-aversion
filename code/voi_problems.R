library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)
library(reshape2)
library(patchwork)
library(future)

source("~/Documents/GitHub/voi-risk-aversion/code/voi_simulations.R")

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

gamma_seq <- seq(-10,10,.1)
## Example 1: 2 actions and two states
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
fcn_VOI_simulation(e1, gamma_seq) %>%
  ggplot(aes(x = gamma, y = VPI, color = factor(max_EU))) +
  geom_line() +
  theme_minimal()

e1_ce <- lapply(gamma_seq, function(gamma) apply(CARA(e1$action_state, gamma) %*% e1$p, 1, 
                                        function(x) CARA_inv(x, gamma))) %>%
  bind_rows()
e1_ce$gamma <- gamma_seq
e1_ce %>%
  pivot_longer(paste0('a', 1:n_actions), names_to = 'actions', values_to = 'ce') %>%
  ggplot(aes(y = ce, x = gamma, color = actions))+
  geom_line()+
  theme_minimal()

## Example 1: 3 actions and 3 states
e1 <- list()
n_states <- 3
n_actions <- 3
a1 <- c(0.689, 0.582, 0.547)
a2 <- c(0.729, 0.674, 0.484)
a3 <- c(0.745, 0.710, 0.332)
e1$p <- c(0.4, 0.2, 0.4)
action_state <- matrix(c(a1, a2, a3), byrow = T, ncol = length(a1))
colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)
e1$action_state <- action_state
fcn_VOI_simulation(e1, gamma_seq) %>%
  ggplot(aes(x = gamma, y = VPI, color = as.factor(max_EU))) +
  geom_line() +
  theme_minimal()

e1_ce <- lapply(gamma_seq, function(gamma) apply(CARA(e1$action_state, gamma) %*% e1$p, 1, 
                                                 function(x) CARA_inv(x, gamma))) %>%
  bind_rows()
e1_ce$gamma <- gamma_seq
e1_ce %>%
  pivot_longer(paste0('a', 1:n_actions), names_to = 'actions', values_to = 'ce') %>%
  ggplot(aes(y = ce, x = gamma, color = actions))+
  geom_line()+
  theme_minimal()

## Several states
set.seed(1111)
n_states <- 50
n_actions <- 50
A <- lhs::randomLHS(n = n_actions, k = n_states)
action_state <- A
p <- runif(n_states)
p <- rep(1, n_states)
p <- p/sum(p)
n_y <- 5 # Partial information experiment that resolves uncertainty to n_y groups of possible outcomes
Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])

## Parameters -------
set.seed(11111)
n_states <- 20
n_actions <- 20



## Draw from uniform distribution like Holden et al (2024) ------
plan(multisession)
unif_action_state_sim <- function() matrix(runif(n_states*n_actions), nrow = n_actions)
unif_plt <- fcn_plot_simulations(unif_action_state_sim)

## Draw from an exponential distribution like Holden et al (2024) ------
exp_action_state_sim <- function() matrix(rexp(n_states*n_actions, 1), nrow = n_actions)
exp_plt <- fcn_plot_simulations(exp_action_state_sim)

## Draw from an poisson distribution with heterogeneous rate parameters ------
pois_action_state_sim <- function() rpois(n_states*n_actions, rep(sample(0:5, n_actions, replace = T), each=n_states)) %>%
  matrix(nrow = n_actions, byrow = T)
pois_plt <- fcn_plot_simulations(pois_action_state_sim)

## Draw from an lognormal distribution like Holden et al (2024) ------
lnorm_action_state_sim <- function() rlnorm(n_states*n_actions, rep(0.1, each=n_states)) %>%
  matrix(nrow = n_actions, byrow = T)
lnorm_plt <- fcn_plot_simulations(lnorm_action_state_sim)

## Draw from a negative lognormal distribution ------
neg_lnorm_action_state_sim <- function() -rlnorm(n_states*n_actions, rep(0.1, each=n_states)) %>%
  matrix(nrow = n_actions, byrow = T)
neg_lnorm_plt <- fcn_plot_simulations(neg_lnorm_action_state_sim)

## Draw from an kurtotic distribution ------
t_action_state_sim <- function() rt(n_states*n_actions, 3) %>%
  matrix(nrow = n_actions, byrow = T)
t_plt <- fcn_plot_simulations(t_action_state_sim)

## No trade-off in mean and variance
mu <- seq(1, 1.5, length.out = n_actions)
sd <- seq(0.5, 0.5, length.out = n_actions)
mv_no_tradeoff_action_state_sim <- function() rnorm(n_states*n_actions, 
                                        rep(mu, each=n_states),
                                        rep(sd, each=n_states)) %>%
  matrix(nrow = n_actions, byrow = T)
mv_no_tradeoff_plt <- fcn_plot_simulations(mv_no_tradeoff_action_state_sim)

## Explicit trade-off in mean and variance
mu <- seq(1, 1.5, length.out = n_actions)
sd <- seq(0.1, 3, length.out = n_actions)
mv_action_state_sim <- function() rnorm(n_states*n_actions, 
                                          rep(mu, each=n_states),
                                          rep(sd, each=n_states)) %>%
  matrix(nrow = n_actions, byrow = T)
mv_plt <- fcn_plot_simulations(mv_action_state_sim)

## Write plots -----
fig1 <- unif_plt$order_plt + ggtitle("Uniform") + unif_plt$v_plt + unif_plt$plt + 
  exp_plt$order_plt + ggtitle("Exponential") + exp_plt$v_plt + exp_plt$plt +
  pois_plt$order_plt + ggtitle("Poisson") + pois_plt$v_plt + pois_plt$plt + plot_layout(guides='collect',byrow = F) & theme(legend.position = "bottom")
fig1
ggsave(fig1, filename = "plots/homogeneous_dist.png", width = 12, height = 10, dpi = 300)

fig2 <- lnorm_plt$order_plt + ggtitle("Lognormal") + lnorm_plt$v_plt + lnorm_plt$plt + 
  neg_lnorm_plt$order_plt + ggtitle("Lognormal (negative values)") + neg_lnorm_plt$v_plt + neg_lnorm_plt$plt + 
  t_plt$order_plt + ggtitle("T-distribution") + t_plt$v_plt + t_plt$plt +
  plot_layout(guides='collect',byrow = F, ncol = 3) & theme(legend.position = "bottom")
fig2
ggsave(fig2, filename = "plots/lnorm_dist.png", width = 12, height = 10, dpi = 300)


fig3 <- mv_no_tradeoff_plt$order_plt + ggtitle("Normal distribution") + mv_no_tradeoff_plt$v_plt + mv_no_tradeoff_plt$plt +
  mv_plt$order_plt + ggtitle("Mean-Variance tradeoff") + mv_plt$v_plt + mv_plt$plt + plot_layout(guides='collect',byrow = F, ncol=2) & theme(legend.position = "bottom")
fig3
ggsave(fig3, filename = "plots/mv_dist.png", width = 10, height = 10, dpi = 300)
