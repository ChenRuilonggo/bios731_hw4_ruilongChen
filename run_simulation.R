#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript scripts/run_simulation.R <n> <task_id> <batch_size>")
}

n <- as.integer(args[1])
task_id <- as.integer(args[2])
batch_size <- as.integer(args[3])

nsim_total <- 500
mu_true <- c(0, 5, 10, 20)
K <- length(mu_true)

sigma2_prior <- 100
gibbs_iter <- 10000
gibbs_burn <- 2000
vb_max_iter <- 1000
vb_tol <- 1e-8

dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", "raw"), showWarnings = FALSE, recursive = TRUE)

# 关键：直接读取你已有的函数文件
source("/projects/YangLabData/Ruilong/731/R/helpers.R")
source("/projects/YangLabData/Ruilong/731/R/gibbs_mixture.R")
source("/projects/YangLabData/Ruilong/731/R/vb_mixture.R")

start_rep <- (task_id - 1) * batch_size + 1
end_rep <- min(task_id * batch_size, nsim_total)

if (start_rep > nsim_total) {
  stop("This task_id is out of range for nsim_total.")
}

simulate_mixture_data <- function(n, mu_true) {
  K <- length(mu_true)
  c_true <- sample(1:K, size = n, replace = TRUE)
  y <- rnorm(n, mean = mu_true[c_true], sd = 1)

  tibble(
    y = y,
    c_true = c_true
  )
}

run_one_sim <- function(n, sim_id, mu_true) {
  set.seed(100000 + n * 1000 + sim_id)

  dat <- simulate_mixture_data(n = n, mu_true = mu_true)
  y <- dat$y

  gibbs_time <- system.time({
    fit_gibbs <- gibbs_mixture(
      y = y,
      K = length(mu_true),
      sigma2 = sigma2_prior,
      n_iter = gibbs_iter,
      burn = gibbs_burn
    )
  })

  vb_time <- system.time({
    fit_vb <- vb_mixture(
      y = y,
      K = length(mu_true),
      sigma2 = sigma2_prior,
      max_iter = vb_max_iter,
      tol = vb_tol
    )
  })

  # Gibbs 95% credible interval
  gibbs_ci_lower <- apply(
    fit_gibbs$mu_draws[fit_gibbs$keep, , drop = FALSE],
    2,
    quantile,
    probs = 0.025
  )

  gibbs_ci_upper <- apply(
    fit_gibbs$mu_draws[fit_gibbs$keep, , drop = FALSE],
    2,
    quantile,
    probs = 0.975
  )

  # VB approximate 95% interval
  vb_ci_lower <- fit_vb$m - 1.96 * sqrt(fit_vb$s2)
  vb_ci_upper <- fit_vb$m + 1.96 * sqrt(fit_vb$s2)

  gibbs_df <- tibble(
    method = "Gibbs",
    n = n,
    sim_id = sim_id,
    component = 1:length(mu_true),
    mu_true = mu_true,
    mu_hat = fit_gibbs$mu_post_mean,
    ci_lower = gibbs_ci_lower,
    ci_upper = gibbs_ci_upper,
    covered = as.integer(mu_true >= gibbs_ci_lower & mu_true <= gibbs_ci_upper),
    elapsed_time = unname(gibbs_time["elapsed"])
  )

  vb_df <- tibble(
    method = "VB",
    n = n,
    sim_id = sim_id,
    component = 1:length(mu_true),
    mu_true = mu_true,
    mu_hat = fit_vb$m,
    ci_lower = vb_ci_lower,
    ci_upper = vb_ci_upper,
    covered = as.integer(mu_true >= vb_ci_lower & mu_true <= vb_ci_upper),
    elapsed_time = unname(vb_time["elapsed"])
  )

  bind_rows(gibbs_df, vb_df)
}

all_results <- vector("list", length = end_rep - start_rep + 1)

idx <- 1
for (sim_id in start_rep:end_rep) {
  message("Running n = ", n, ", sim_id = ", sim_id)
  all_results[[idx]] <- run_one_sim(n = n, sim_id = sim_id, mu_true = mu_true)
  idx <- idx + 1
}

results_batch <- bind_rows(all_results)

outfile <- file.path(
  "results",
  "raw",
  paste0("sim_n", n, "_task", task_id, ".csv")
)

write_csv(results_batch, outfile)

message("Saved results to: ", outfile)