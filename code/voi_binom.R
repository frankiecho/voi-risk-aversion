
S <- 200
n <- 200
p <- 0.5

# Give mu and np and get n and p
gamma_seq = seq(0.01, 2, 0.01)
pref <- 'CE'

n_states <- S*2
n_actions <- 2

a11 <- rbinom(S, 100/p, p)
a21 <- rbinom(S, 100/p, p)
a12 <- rbinom(S, 135/p, p)
a22 <- rbinom(S, 55/p, p)

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

e1_voi <- fcn_VOI_simulation(e1, pref = pref, gamma_seq = gamma_seq)

e1_vpi_plt <- fcn_plt_vpi(e1_voi, n_actions)
e1_vpi_plt

e1_ce <- lapply(gamma_seq, function(gamma) apply(e1$action_state, 1, pref_define(pref), lambda=gamma, p = e1$p)) %>%
  bind_rows()
e1_ce$gamma <- gamma_seq
e1_ce_plt <- fcn_plt_ce(e1_ce, n_actions)
e1_ce_plt


## Example 1: 3 actions and 3 states
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

n_individuals <- 2

a1 <- do.call(c, lapply(pa1, \(p) rbinom(S, n_individuals, p)))
a2 <- do.call(c, lapply(pa2, \(p) rbinom(S, n_individuals, p)))
a3 <- do.call(c, lapply(pa3, \(p) rbinom(S, n_individuals, p)))

action_state <- matrix(c(a1, a2, a3), byrow = T, ncol = n_states)
colnames(action_state) <- paste0('s', 1:n_states)
rownames(action_state) <- paste0('a', 1:n_actions)
e2$p <- rep(1/n_states, n_states)
e2$action_state <- action_state
e2_voi <- fcn_VOI_simulation(e2, pref = pref, gamma_seq = gamma_seq)
e2_vpi_plt <- fcn_plt_vpi(e2_voi, n_actions)

e2_ce <- lapply(gamma_seq, function(g) apply(e2$action_state, 1, pref_define(pref), lambda=g, p = e2$p)) %>%
  bind_rows()
e2_ce$gamma <- gamma_seq
e2_ce_plt <- fcn_plt_ce(e2_ce, n_actions)

e2_vpi_plt
e2_ce_plt
