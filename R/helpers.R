sample_categorical <- function(prob) {
  sample.int(length(prob), size = 1, prob = prob)
}

relabel_by_mu <- function(mu, c_vec) {
  ord <- order(mu)
  mu_new <- mu[ord]
  
  old_to_new <- integer(length(mu))
  old_to_new[ord] <- seq_along(ord)
  
  c_new <- old_to_new[c_vec]
  
  list(mu = mu_new, c = c_new)
}