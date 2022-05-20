#!/usr/bin/env Rscript
source("sim.R")
library(pbapply)
set.seed(42)
conv_rate <- function(sim_type, n_sim, n_obs, p_c) {
  length(unique((simulate_ivpim(sim_type, n_sim, n_obs = n_obs, p_c = p_c))$index)) / n_sim
}
# conv_rate("obs", 10, 100, 0.5)

n_sim <- 2000
n_obs <- (1:10) * 40
p_c <- (1:10 / 10)
sim_combs <- expand.grid(n_obs, p_c)
names(sim_combs) <- c("n_obs", "p_c")
rct_conv_rates <- apply(sim_combs, 1, function(x) {
  conv_rate("rct", n_sim, as.numeric(x[1]), as.numeric(x[2]))
})
obs_conv_rates <- apply(sim_combs, 1, function(x) {
  conv_rate("obs", n_sim, as.numeric(x[1]), as.numeric(x[2]))
})
ivpim_conv_rates <- data.frame(sim_combs, rct_conv_rates, obs_conv_rates)
write.csv(ivpim_conv_rates, "./out/ivpim_conv_rates.csv", row.names = F)

q("no")
