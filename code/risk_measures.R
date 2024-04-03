## Plot preference functions
library(ggplot2)
library(ggpubr)
library(patchwork)

pref_list_long <- c("CE"="Certainty Equivalent", "MV"="Mean-Variance", "MCVaR"="Mean-Expected Shortfall")
pref_list <- names(pref_list_long)

plt_list <- list()
for (pref in pref_list) {
  pref_func <- pref_define(pref)
  N <- 100000
  p <- rep(1, N)/N
  sigma_vec <- 0:3
  if (pref == 'CE') {
    l_vec <- seq(-6,6,0.1)
  } else {
    l_vec <- seq(-1,1,0.1)
  }
  norm_dist <- lapply(sigma_vec, \(x) rnorm(N, 1, x))
  outcome <- sapply(l_vec, \(l) sapply(norm_dist, \(x) pref_func(x, lambda=l))) %>%
    t() %>%
    as.data.frame()
  names(outcome) <- as.character(sigma_vec)
  outcome$lambda <- l_vec
  plt_list[[pref]] <- outcome %>%
    pivot_longer(cols = all_of(as.character(sigma_vec))) %>%
    ggplot(aes(x = lambda, y = value, color = name)) +
    geom_hline(yintercept = 0, color = 'gray50') +
    geom_vline(xintercept = 0, color = 'gray50') +
    annotate("text", x = ifelse(pref=='CE',-3, -0.5), y = 14, label = "Risk-loving", vjust = 0, hjust = 0.5)+
    annotate("text", x = ifelse(pref=='CE', 3, 0.5), y = 14, label = "Risk-averse", vjust = 0, hjust = 0.5) +
    scale_color_manual("Standard \nDeviation", values=okabe_ito_colors[c(1,3,5,7)]) +
    geom_line(linewidth = 1) +
    theme_pubr() +
    scale_x_continuous("Risk Aversion Coefficient") +
    scale_y_continuous("Risk-adjusted Value", limits = c(-15, 15)) +
    theme(legend.position = 'right', panel.border = element_rect(linewidth = 1, fill = NA)) +
    ggtitle(pref_list_long[pref])
  if (pref != pref_list[1]) {
    plt_list[[pref]] <- plt_list[[pref]] + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

plt_out <- wrap_plots(plt_list) + plot_layout(guides='collect')
ggsave("plots/sd_plt.png", plt_out, width = 9, height = 4)
