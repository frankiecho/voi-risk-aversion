library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)
library(reshape2)

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

gamma_seq <- seq(-20,20,.1)
## Example 1: 2 actions and two states
e1 <- list()
n_states <- 2
n_actions <- 2
a1 <- c(1, 1)
a2 <- c(1.5, 0.9)
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

## Draw from an inverse Wishart distribution
set.seed(11111)
n_states <- 100
n_actions <- 20

nsims <- 100
mu <- runif(n_actions)
set.seed(1224598)
fcn_wish_sim <- function(df = 1, n_y = 5) {
  Sigma <- LaplacesDemon::rinvwishart(n_actions+df, diag(n_actions))
  prob <- list(mu = mu, sigma = Sigma)
  action_state <- t(rmvnorm(n_states, mean = mu, sigma = Sigma))
  colnames(action_state) <- paste0('s', 1:n_states)
  rownames(action_state) <- paste0('a', 1:n_actions)
  prob$action_state <- action_state
  p <- runif(n_states)
  prob$p <- p/sum(p)
  prob$Y <- sort(rep(1:n_y, 1+n_states/n_y)[1:n_states])
  prob$evpxi <- T
  prob
}
wish_sim <- replicate(nsims, fcn_wish_sim(), simplify = F)
wish_sim_results <- lapply(wish_sim, fcn_VOI_simulation, gamma_seq = gamma_seq)
names(wish_sim_results) <- 1:nsims
wish_sim_table <- wish_sim_results %>%
  #lapply(function(x) mutate(x, VPI = VPI - x$VPI[x$gamma == 0])) %>%
  #lapply(function(x) x$VPI) %>%
  bind_rows(.id = 'run_index')
wish_sim_mean <- wish_sim_table %>%
  group_by(gamma) %>%
  summarise(VPI = mean(VPI, na.rm = T))
wish_sim_table %>%
  #ggplot(aes(x = gamma, y = VPI, group = name)) +
  ggplot(aes(x = gamma, y = VPI)) +
  geom_line(aes(group = run_index), size = 1, alpha = 0.1) +
  geom_line(data= wish_sim_mean, color = 'blue', size = 1) +
  coord_cartesian(ylim = c(0,3), xlim = c(-15, 15)) +
  theme_minimal()

# Fineness of partition ----------------
fineness <- rep(c(5,10,25,50), 50) %>% sort()
wish_sim_evpxi <- lapply(fineness, function(ny) fcn_wish_sim(n_y = ny))
wish_sim_evpxi_results <- lapply(wish_sim_evpxi, fcn_VOI_simulation, gamma_seq = gamma_seq)
names(wish_sim_evpxi_results) <- paste0(1:length(fineness), '_', fineness)
wish_sim_evpxi_results %>%
  bind_rows(.id='name') %>%
  separate(name, c('run_index', 'fineness')) %>%
  group_by(fineness, gamma) %>%
  summarise(VPI = mean(VPI)) %>%
  ggplot(aes(y = VPI, x = gamma, color = as.numeric(fineness), group = as.numeric(fineness))) +
  geom_line(size=1)+
  theme_minimal()


# Degrees of freedom -----------
gamma_seq <- seq(-20,20,.5)
df_seq <- rep(c(1,5,10), 100) %>% sort()
wish_sim_df <- lapply(df_seq, fcn_wish_sim)
names(wish_sim_df) <- paste(1:length(df_seq), df_seq, sep = '_')
wish_sim_df_results <- lapply(wish_sim_df, fcn_VOI_simulation, gamma_seq = gamma_seq) %>%
  bind_rows(.id='name') %>%
  separate(name, c('run_id', 'df'), sep = "_") #%>%
wish_sim_df_mean <- wish_sim_df_results %>%
  group_by(df, gamma) %>%
  summarise(VPI = mean(VPI) )
wish_sim_df_results %>%
  #ggplot(aes(x = gamma, y = VPI, group = name)) +
  ggplot(aes(x = gamma, y = VPI)) +
  geom_line(aes(group = as.numeric(run_id)), size = 1, alpha = 0.1) +
  geom_line(data = wish_sim_df_mean, size = 1, color = 'red') +
  facet_wrap(~df) +
  theme_minimal()

# Find the gamma value where the value of information is the highest
wish_sim_df_max_gamma <- wish_sim_df_results %>%
  group_by(run_id, df) %>%
  summarise(max_gamma = gamma[which.max(VPI)]) %>%
  ungroup()

ggplot(wish_sim_df_max_gamma, aes(x = max_gamma)) +
  geom_histogram() +
  facet_wrap(~df) +
  theme_minimal()


