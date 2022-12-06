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

#' @title markow_model
#'
#' @description Generate a random sequence in the markow model with a specific lenght over a specified alphabet
#' with a specified transition matrix and initial probabilities. Given an index set \eqn{I} with \eqn{r} elements then a \eqn{r\times r} matrix \eqn{(p_{ij})_{i,j \in I}}
#' of nonnegative numbers is a transition matrix if \eqn{\forall i \in I: \sum_{j \in I} p_{ij} = 1}. \eqn{p_{ij}} is the probability that a state \eqn{i} is followed
#' by state \eqn{j}. Given an Alphabet \eqn{\Sigma}, nonnegative numbers \eqn{(p_i)_{i \in \Sigma}} with sum \eqn{1} and a transition matrix \eqn{P = (p_{ij})_{i,j \in \Sigma}}
#' the probability \eqn{\mathbb{P}(\textbf{s})} of a sequence \eqn{(\textbf{s})} over \eqn{\Sigma} with length \eqn{n} is calculated by
#' \eqn{\mathbb{P}(\textbf{s}) = p_{s_1}\prod_{l = 2}^n p_{s_{n-1}s_n}}
#' @param alphabet The alphabet the sequence should be based on as a vector
#' @param length The length of the generated sequence as an integer
#' @param matrix The transition matrix so that \code{matrix[i, j]} describes the probability that a state \eqn{i} is followed by state \eqn{j}
#' @param initial_probabilities A vector of probabilities that a sequence is starting with the \eqn{i}-th letter of the alphabet. If this argument is empty
#' an equal distribution is assumed (which does not matter for large sequences)
#' @return A vector that contains the letters of the generated sequence. If either the sum of the rows of the transition matrix or the sum of the initial probabilities
#' are not equal to 1 a warning is thrown
markow_model = function(alphabet, length, matrix, initial_probabilities = NULL ) {

  #update rownames and colnames for correct indexing
  rownames(matrix) = alphabet
  colnames(matrix) = alphabet

  if(is.null(initial_probabilities)) {
    initial_probabilities = rep(1 / length(alphabet), length(alphabet))
  }

  #alert if sum of probabilities are not equal to 1
  if(sum(initial_probabilities) != 1) {
    warning("Sum of initial_probabilities != 1")
  }

  if(prod(rowSums(matrix) == rep(1, length(a))) != 1) {
    warning("Sum of one ore more row in the transition matrix != 1")
  }

  res = character(length)
  res[1] = sample(alphabet, 1, prob = initial_probabilities)
  for(i in 2:length) {
    res[i] = sample(alphabet, 1, prob = matrix[res[i - 1],])
  }

  return(res)
}
