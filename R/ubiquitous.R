
equal_distribution_model = function(alphabet, length) {
  ps = rep(1 / length(alphabet), length(alphabet))
  res = sample(alphabet, length, replace = TRUE, prob = ps)
  return(res)
}

multinomial_distribution_model = function(alphabet, length, probabilities) {
  res = sample(alphabet, length, replace = TRUE, prob = probabilities)
  return(res)
}
