library(upim)
library(complyr)
sim_one_rct <- function(n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "logit") {
  if (model == "probit") {
    error_dist <- "normal"
  } else if (model == "logit") {
    error_dist <- "gumbel"
  } else {
    stop("Model must be either probit of logit!")
  }
  df <- dgp_rct(n = n_obs, p_c = p_c, alpha = alpha, beta_1 = beta_1, beta_2 = beta_2, beta_3 = 0, error_dist = error_dist)
  # We use three methods to analyze this data set: ITT, PP, IV
  # ITT: intention-to-treat
  itt_fit <- pim_fit(y = df$y, X = df[, c("z", "x1", "x2")], link = model)
  # PP: per-protocol
  df_pp <- df[df$z == df$a, ]
  pp_fit <- pim_fit(y = df_pp$y, X = df_pp[, c("a", "x1", "x2")], link = model)
  # IV: instrumental variable
  ps_model <- glm(df$z ~ 1, family = binomial(link = "logit"), x = TRUE)
  iv_fit <- complyr::ivpim(y = df$y, z = df$z, a = df$a, X = df[, c("x1", "x2")],
                           ps_model =  ps_model, link = "logit")
  sum_stat <- function(model_fit) {
    # get summary stats for downstream analysis
    ss <- c(model_fit$termcd, model_fit$coef,
            sqrt(diag(model_fit$vcov)))
    names(ss) <- c("termcd", "alpha_hat", "beta_1_hat", "beta_2_hat",
                   "alpha_se", "beta_1_se", "beta_2_se")
    ss
  }
  results <- as.data.frame(rbind(sum_stat(itt_fit),
                                 sum_stat(pp_fit),
                                 sum_stat(iv_fit)))
  results$methods <- c("itt", "pp", "iv")
  results$n_obs <- n_obs
  results$p_c <- p_c
  results$alpha <- alpha
  results$beta_1 <- beta_1
  results$beta_2 <- beta_2
  results$model <- model
  results
}

sim_one_obs <- function(n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "logit") {
  if (model == "probit") {
    error_dist <- "normal"
  } else if (model == "logit") {
    error_dist <- "gumbel"
  } else {
    stop("Model must be either probit of logit!")
  }
  df <- dgp_obs(n = n_obs, p_c = p_c, alpha = alpha, beta_1 = beta_1, beta_2 = beta_2, beta_3 = 0, error_dist = error_dist)
  # fit propensity score models
  ps_model_correct <- glm(z ~ x1_t + x2_t, family = binomial(link = "logit"), data = df, x = T)
  ps_model_misspecified <- glm(z ~ x1 + x2, family = binomial(link = "logit"), data = df, x = T)
  ps_model_marginal <- glm(z ~ 1, family = binomial(link = "logit"), data = df, x = T)
  # We benchmark IVPIMs with different PS models to PP analysis
  # PP: per-protocol
  df_pp <- df[df$z == df$a, ]
  pp_fit <- pim_fit(y = df_pp$y, X = df_pp[, c("a", "x1", "x2")], link = model)
  # IVPIM
  iv_correct_fit <- ivpim(y = df$y, z = df$z, a = df$a, X = model.matrix(~ x1 + x2 - 1, df),
                              ps_model = ps_model_correct, link = model)
  iv_misspecified_fit <- ivpim(y = df$y, z = df$z, a = df$a, X = model.matrix(~ x1 + x2 - 1, df),
                                   ps_model = ps_model_misspecified, link = model)
  iv_marginal_fit <- ivpim(y = df$y, z = df$z, a = df$a, X = model.matrix(~ x1 + x2 - 1, df),
                               ps_model = ps_model_marginal, link = model)
  sum_stat <- function(model_fit) {
    # get summary stats for downstream analysis
    ss <- c(model_fit$termcd, model_fit$coef,
            sqrt(diag(model_fit$vcov)))
    names(ss) <- c("termcd", "alpha_hat", "beta_1_hat", "beta_2_hat", "alpha_se", "beta_1_se", "beta_2_se")
    ss
  }
  results <- as.data.frame(rbind(sum_stat(pp_fit),
                                 sum_stat(iv_correct_fit),
                                 sum_stat(iv_misspecified_fit),
                                 sum_stat(iv_marginal_fit)))
  results$methods <- c("pp", "iv-correct", "iv-misspecified", "iv-marginal")
  results$n_obs <- n_obs
  results$p_c <- p_c
  results$alpha <- alpha
  results$beta_1 <- beta_1
  results$beta_2 <- beta_2
  results$model <- model
  results
}

simulate_ivpim <- function(sim_one = "rct", n_sim = 5, ...) {
  # use pbreplicate (replicate with progress bar)
  # to repeatedly simulate n_sim times
  sim_args <- list(...) # ... stores arguments to sim_one_*() functions
  print(c(sim_type = sim_one, n_sim = n_sim, unlist(sim_args)))
  sim_result_list <- pbapply::pbreplicate(n_sim, { # replicate n_sim times
    regular_iter <- FALSE # flag to mark regular iteration
    while(!regular_iter) {
      try({
        if (sim_one == "rct") {
          one_sim <- suppressWarnings(do.call(sim_one_rct, sim_args))
        } else if (sim_one == "obs") {
          one_sim <- suppressWarnings(do.call(sim_one_obs, sim_args))
        } else { # default is rct, can add other sim_one types in the future
          one_sim <- suppressWarnings(do.call(sim_one_rct, sim_args))
        }
        if (all(one_sim$termcd == 1)) { # regular if all termcd == 1
          regular_iter <- TRUE
        }
      }, silent = TRUE)
    }
    # return the simulation result of this iteration, strip termcd (1st column)
    sim_result <- one_sim[, -1]
  }, simplify = F)
  # add simulation index (each index corresponds to the same data frame), then return all simulation results
  all_sim_result <- do.call(rbind, sim_result_list)
  all_sim_result$index <- rep(seq(n_sim), each = NROW(all_sim_result) / n_sim)
  all_sim_result
}

# simulate_ivpim(sim_one = "rct", n_sim = 10, n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "logit")
# simulate_ivpim(sim_one = "obs", n_sim = 10, n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "logit")
#
# simulate_ivpim(sim_one = "rct", n_sim = 10, n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "probit")
# simulate_ivpim(sim_one = "obs", n_sim = 10, n_obs = 200, p_c = 0.7, alpha = 1, beta_1 = 0.5, beta_2 = -0.7, model = "probit")


