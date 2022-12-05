
equal_distribution_model = function(alphabet, length) {
  ps = rep(1 / length(alphabet), length(albphabet))
  sample(alphabet, length, replace = TRUE, prob = ps)
}
