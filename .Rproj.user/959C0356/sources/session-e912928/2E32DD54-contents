gibbs_mixture <- function(y,
                          K = 2,
                          sigma2 = 100,
                          n_iter = 10000,
                          burn = 2000,
                          init_c = NULL) {
  
  n <- length(y)
  
  if (is.null(init_c)) {
    init_c <- sample(1:K, size = n, replace = TRUE)
  }
  
  c_curr <- init_c
  
  mu_curr <- rep(mean(y), K)
  
  mu_draws <- matrix(NA_real_, nrow = n_iter, ncol = K)
  c_draws  <- matrix(NA_integer_, nrow = n_iter, ncol = n)
  
  for (iter in seq_len(n_iter)) {
    
    # update cluster assignments
    for (i in seq_len(n)) {
      
      log_prob <- sapply(seq_len(K), function(k) {
        -0.5 * (y[i] - mu_curr[k])^2
      })
      
      log_prob <- log_prob - max(log_prob)
      
      prob <- exp(log_prob)
      prob <- prob / sum(prob)
      
      c_curr[i] <- sample_categorical(prob)
    }
    
    # update means
    for (k in seq_len(K)) {
      
      idx <- which(c_curr == k)
      
      n_k <- length(idx)
      S_k <- sum(y[idx])
      
      post_var <- 1 / (n_k + 1 / sigma2)
      post_mean <- post_var * S_k
      
      mu_curr[k] <- rnorm(
        1,
        mean = post_mean,
        sd = sqrt(post_var)
      )
    }
    
    relabeled <- relabel_by_mu(mu_curr, c_curr)
    
    mu_curr <- relabeled$mu
    c_curr  <- relabeled$c
    
    mu_draws[iter, ] <- mu_curr
    c_draws[iter, ]  <- c_curr
  }
  
  keep <- (burn + 1):n_iter
  
  mu_post_mean <- colMeans(mu_draws[keep, , drop = FALSE])
  
  c_post_mode <- apply(
    c_draws[keep, , drop = FALSE],
    2,
    function(x) {
      as.integer(names(sort(table(x), decreasing = TRUE))[1])
    }
  )
  
  list(
    mu_draws = mu_draws,
    c_draws = c_draws,
    mu_post_mean = mu_post_mean,
    c_post_mode = c_post_mode,
    keep = keep
  )
}