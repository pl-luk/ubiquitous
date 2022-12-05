#' @title equal_distribution_model
#'
#' @description Generate a random sequence in the equal distribution model with a specific length over a specified alphabet.
#' The probability \eqn{\mathbb{P}(\textbf{s})} of a sequence \eqn{\textbf{s}} over an alphabet \eqn{\Sigma} with length \eqn{n} is calculated by
#' \eqn{\mathbb{P}(\textbf{s}) = \frac{1}{|\Sigma|^n}}
#' @param alphabet The alphabet the sequence should be based on as a vector
#' @param length The length of the generated sequence as an integer
#' @return A vector that contains the letters of the generated sequence
equal_distribution_model = function(alphabet, length) {
  ps = rep(1 / length(alphabet), length(alphabet))
  res = sample(alphabet, length, replace = TRUE, prob = ps)
  return(res)
}

#' @title multinomial_distribution_model
#'
#' @description  Generate a random sequence in the multinomial distribution model with a specific length over a specified alphabet
#' with specified probabilities. Given an alphabet \eqn{\Sigma} and nonnegative Numbers \eqn{(p_i)_{i \in \Sigma}} with sum \eqn{1}
#' (\eqn{\forall i \in \Sigma: p_i \in [0, 1]}) the probability \eqn{\mathbb{P}(\textbf{s})} of a sequence \eqn{\textbf{s}} over \eqn{\Sigma} with length \eqn{n}
#' is calculated by \eqn{\mathbb{P}(\textbf{s}) = \prod_{l = 1}^{n} p_{s_l}}
#' @param alphabet The alphabet the sequence should be based on as a vector
#' @param length The length of the generated sequence as an integer
#' @param probabilities A vector of probabilities so that the \eqn{i}-th element is the probability of the \eqn{i}-th letter of the alphabet.
#' @return A vector that contains the letters of the generated sequence
multinomial_distribution_model = function(alphabet, length, probabilities) {
  res = sample(alphabet, length, replace = TRUE, prob = probabilities)
  return(res)
}
