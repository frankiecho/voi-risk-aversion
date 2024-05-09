library(patchwork)

source("code/voi_simulations.R")

set.seed(103905)

S <- 1000 # Number of draws from binomial distribution

gamma_seq = seq(0.01, 10, 0.01)
#gamma_seq[1] <- gamma_seq[1] + 1e-10
pref <- 'CE'

# Variation of Canessa (2015) Example 1 of Chytrid, with number of turtles modelled as binomial distribution with the same mean as Canessa
n <- 20

n_states <- S*2
n_actions <- 2

a11 <- rbinom(S, n, 0.5)
a21 <- rbinom(S, n, 0.5)
a12 <- rbinom(S, n, 0.675)
a22 <- rbinom(S, n, 0.275)

action_state <- matrix(c(a11, a21, a12, a22), byrow = T, ncol = S*2)
colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)

Y <- list(c(1:S), c((S+1):(S*2)))
pY <- c(0.5, 0.5)

e1 <- list()
e1$action_state <- action_state
e1$p <- rep(1/(S*2), S*2)
e1$Y <- Y
e1$pY <- c(0.5, 0.5)
e1$evsi <- TRUE

e1_voi <- fcn_VOI_simulation(e1, pref = 'CE', gamma_seq = gamma_seq)

e1_vpi_plt <- fcn_plt_vpi(e1_voi, n_actions)+ geom_vline(xintercept = 1)
e1_vpi_plt

e1_ce <- lapply(gamma_seq, function(gamma) apply(e1$action_state, 1, pref_define(pref), lambda=gamma, p = e1$p)) %>%
  bind_rows()
e1_ce$V_certainty <- e1_voi$V_certainty
e1_ce$gamma <- gamma_seq
e1_ce_plt <- fcn_plt_ce(e1_ce, n_actions)+ 
  geom_line(data = e1_voi, aes(y = V_uncertainty, x = gamma), color = 'white', linetype = 2, linewidth = 1) +
  geom_vline(xintercept = 1)
e1_ce_plt


# Variation of Canessa (2015) Example 2 of turtle reintroduction, where p is the prob success parameter of respective binomial distributions
e2 <- list()
n_states <- S*3
n_actions <- 3
pa1 <- c(0.689, 0.582, 0.547)
pa2 <- c(0.729, 0.674, 0.484)
pa3 <- c(0.745, 0.710, 0.332)
e2$pY <- c(0.4, 0.2, 0.4)
e2$Y <- list(c(1:S), S+c(1:S), S*2+c(1:S))
e2$p <- rep(1/(S*3), S*3)
e2$evsi <- TRUE

n_individuals <- 10

a1 <- do.call(c, lapply(pa1, \(p) rbinom(S, n_individuals, p)))
a2 <- do.call(c, lapply(pa2, \(p) rbinom(S, n_individuals, p)))
a3 <- do.call(c, lapply(pa3, \(p) rbinom(S, n_individuals, p)))

action_state <- matrix(c(a1, a2, a3), byrow = T, ncol = n_states)
colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)
e2$p <- rep(1/n_states, n_states)
e2$action_state <- action_state
e2_voi <- fcn_VOI_simulation(e2, pref = pref, gamma_seq = gamma_seq)
e2_vpi_plt <- fcn_plt_vpi(e2_voi, n_actions) + geom_vline(xintercept = 1)

e2_ce <- lapply(gamma_seq, function(g) apply(e2$action_state, 1, pref_define(pref), lambda=g, p = e2$p)) %>%
  bind_rows()
e2_ce$gamma <- gamma_seq
e2_ce$V_certainty <- e2_voi$V_certainty
e2_ce_plt <- fcn_plt_ce(e2_ce, n_actions) + 
  geom_line(data = e2_voi, aes(y = V_certainty, x = gamma), color = 'gray30', linewidth = 1) +
  geom_line(data = e2_voi, aes(y = V_uncertainty, x = gamma), color = 'white', linetype = 2, linewidth = 1) +
  geom_vline(xintercept = 1)
e2_ce_plt

(e1_ce_plt / e1_vpi_plt) | (e2_ce_plt / e2_vpi_plt)

