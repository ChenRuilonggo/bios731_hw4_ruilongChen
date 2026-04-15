compute_elbo <- function(y, m, s2, phi, sigma2) {
  
  n <- length(y)
  K <- length(m)
  
  term1 <- sum(
    -0.5 * log(2 * pi * sigma2) -
      (m^2 + s2) / (2 * sigma2)
  )
  
  term2 <- -n * log(K)
  
  term3 <- 0
  
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      term3 <- term3 + phi[i, k] * (
        -0.5 * log(2 * pi) -
          0.5 * ((y[i] - m[k])^2 + s2[k])
      )
    }
  }
  
  term4 <- sum(
    0.5 * log(2 * pi * exp(1) * s2)
  )
  
  eps <- 1e-12
  
  term5 <- -sum(phi * log(phi + eps))
  
  term1 + term2 + term3 + term4 + term5
}


vb_mixture <- function(y,
                       K = 2,
                       sigma2 = 100,
                       max_iter = 1000,
                       tol = 1e-8) {
  
  n <- length(y)
  
  q <- quantile(
    y,
    probs = seq(0.25, 0.75, length.out = K)
  )
  
  m <- as.numeric(q)
  s2 <- rep(1, K)
  
  phi <- matrix(
    1 / K,
    nrow = n,
    ncol = K
  )
  
  elbo <- numeric(max_iter)
  
  for (iter in seq_len(max_iter)) {
    
    # update phi
    for (i in seq_len(n)) {
      
      log_phi <- sapply(seq_len(K), function(k) {
        y[i] * m[k] -
          0.5 * (m[k]^2 + s2[k])
      })
      
      log_phi <- log_phi - max(log_phi)
      
      phi[i, ] <- exp(log_phi)
      phi[i, ] <- phi[i, ] / sum(phi[i, ])
    }
    
    # update mu
    for (k in seq_len(K)) {
      
      s2[k] <- 1 / (
        1 / sigma2 + sum(phi[, k])
      )
      
      m[k] <- s2[k] * sum(phi[, k] * y)
    }
    
    ord <- order(m)
    
    m <- m[ord]
    s2 <- s2[ord]
    phi <- phi[, ord, drop = FALSE]
    
    elbo[iter] <- compute_elbo(
      y = y,
      m = m,
      s2 = s2,
      phi = phi,
      sigma2 = sigma2
    )
    
    if (iter > 1) {
      if (abs(elbo[iter] - elbo[iter - 1]) < tol) {
        elbo <- elbo[1:iter]
        break
      }
    }
  }
  
  c_hat <- max.col(phi)
  
  list(
    m = m,
    s2 = s2,
    phi = phi,
    c_hat = c_hat,
    elbo = elbo
  )
}