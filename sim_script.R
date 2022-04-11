#!/usr/bin/env Rscript
source("sim.R")

args = commandArgs(trailingOnly = TRUE)
# args = c(sim_type, seed, n_sim, p_c, alpha, beta_1, beta_2, model)
# e.g.
# when you run with
# % Rscript --vanilla sim_script.R rct 42 2000 0.6 1 0.5 -0.7 logit
# then
# args = c("rct", "42", "2000", "0.6", "1", "0.5", "-0.7", "logit")

# progress bar setup
library(pbapply)
pboptions(type = "timer")

# fix n_obs = c(200, 400, 800, 1600)
# if you think you need more simulation, i.e., increase n_sim,
# change seed and re-simulate, then append the data frame to previous one
n_obs <- c(200, 400, 800, 1600)
sim_script <- function(sim_type, seed, n_sim, p_c, alpha, beta_1, beta_2, model) {
  set.seed(as.numeric(seed))
  sim_results <- do.call(rbind, lapply(n_obs, function(n) {
    print(Sys.time())
    simulate_ivpim(sim_one = sim_type, n_sim = n_sim, n_obs = n, p_c = p_c,
                   alpha = alpha, beta_1 = beta_1, beta_2 = beta_2, model = model)
  }))
  file_name <- paste0(c("sim_type", "seed", "n_sim", "p_c", "alpha", "beta_1", "beta_2", "model"),
                      "=", args, collapse = ";", ".csv")
  write.csv(sim_results, paste0("./out/", file_name), row.names = F)
}

# call the sim_script() function above by passing args
sim_script(
  args[1], # sim_type
  as.numeric(args[2]), # seed
  as.numeric(args[3]), # n_sim
  as.numeric(args[4]), # p_c
  as.numeric(args[5]), # alpha
  as.numeric(args[6]), # beta_1
  as.numeric(args[7]), # beta_2
  args[8] # model
)

q("no")
