library(seqinr)
library(Biostrings)
library(stringr)

AMINO_ALPHABET = c(names(AMINO_ACID_CODE), "*")
DNA_ALPHABET = c("A", "C", "G", "T")

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
#' @return A vector that contains the letters of the generated sequence. If sum of probabilities is not equal to 1 a warning is thrown
multinomial_distribution_model = function(alphabet, length, probabilities) {

  if(sum(probabilities) != 1) {
    warning("Sum of probabilities != 1")
  }

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

#' @title rel_freq
#'
#' @description Get the relative frequencies of words (pairs) in a given sequence with a specified alphabet.
#' This can also used as a approximation for the required probabilities for the markow model or the multinomial distribution model. The relative frequency
#' \eqn{\hat{p_i}} for a given state \eqn{i} in a sequence \eqn{\textbf{s}_{[1, L]}} with length \eqn{L}
#' is given by \eqn{\hat{p_i} = \frac{\text{\# } i \text{ in }\textbf{s}_{[1, L]}}{L}} whilst the relative frequency of \eqn{\hat{p_{ij}}} for
#' a given state \eqn{i} that is followed by state \eqn{j} in \eqn{\textbf{s}} is determined by
#' \eqn{\hat{p_{ij}} = \frac{\text{\# } ij \text{ in } \textbf{s}_{[1, L]}}{\text{\# } i \text{ in }\textbf{s}_{[1, L - 1]}}}. If the single \eqn{\hat{p_i}} are
#' clearly differing than an equal distribution can be ruled out.
#' If the \eqn{\hat{p_{ij}}} are clearly differing from the single \eqn{\hat{p_j}} a multinomial distribution can be ruled out.
#'
#' @param sequence The sequence to be analysed
#' @param alphabet The underlying alphabet of the input sequence
#' @param wordsize The size of the words to be analysed
#' @param start The start of the sequence
#' @param count Switch if absolute instead of relative frequencies are desired
#'
#' @returns A vector containing the absolute or relative frequencies of the possible words with length \code{wordsize} found in the input sequence
rel_freq = function(sequence, alphabet, wordsize = 1, start = 0, count = FALSE) {
  res = count(sequence, wordsize, start, freq = !count, alphabet = alphabet)
  return(res)
}

#' @title initialize_markow
#'
#' @description Calculates the transition matrix assuming a markow model specified by a given sequence with the respective alphabet. The exact calculation
#' is described in \link{rel_freq}
#'
#' @param sequence The sequence the markow model should be based on
#' @param alphabet The alphabet that is used by the sequence
#' @param round The number of digits the contents of the transition matrix should be rounded to. If this value is zero or lower the values will stay as calculated.
#'
#' @return The calculated transition matrix
initialize_markow = function(sequence, alphabet, round = 2) {
  Fdi = rel_freq(sequence, alphabet, wordsize = 2, count = TRUE)
  Fmo = rel_freq(sequence, alphabet, count = TRUE)

  matrix_di = matrix(Fdi, length(alphabet), length(alphabet), byrow = TRUE)
  matrix_mo = matrix(Fmo, length(alphabet), length(alphabet), byrow = FALSE)

  res = matrix_di / matrix_mo

  if(round > 0) {
    res = round(res, round)
  }

  colnames(res) = alphabet
  rownames(res) = alphabet

  return(res)
}

#' @title initialize_markow
#'
#' @description Calculates the probability vector assuming a multinomial distribution model specified by a given sequence with the respective alphabet.
#' The exact calculation is described in \link{rel_freq}
#'
#' @param sequence The sequence the multinomial distribution model should be based on
#' @param alphabet The alphabet that is used by the sequence
#' @param round The number of digits the contents of the probability vector should be rounded to. If this value is zero or lower the values will stay as calculated.
#'
#' @return The calculated probability vector
initialize_multinomial = function(sequence, alphabet, round = 2) {
  res = rel_freq(sequence, alphabet)

  if(round > 0) {
    res = round(res, round)
  }

  return(res)
}

#' @title gc_plot
#'
#' @description Calculate the GC content of a given sequence in frames with a given length and plot the GC content as a function of the \eqn{i}-th frame.
#' The way this is calculated is described in \code{\link{rel_freq}}. By viewing the plot it is possible to estimate if GC contents are random or not. Given
#' the relative frequencies \eqn{f} of a variable \eqn{x} in a frame with length \eqn{l} the true and unknown probability of \eqn{x} \eqn{p} can be described as
#' an interval \eqn{p \in \left[f - \frac{2\sigma}{\sqrt{l}}, f + \frac{2\sigma}{\sqrt{l}} \right]} wheras \eqn{\sigma} is the standard deviation of \eqn{x}.
#' The interval describes the \eqn{95\%} confidence interval thus in \eqn{5\%} of all cases the true value of \eqn{p} is element of the interval. Since the true
#' standard diviation \eqn{\sigma = \sqrt(p(1-p))} cannot be calculated it has to be estimated by \eqn{\sigma = \sqrt{f(1 - f)}}. If used on the GC contents the
#' random probability is given by \eqn{p = \frac{1}{4} + \frac{1}{4} = \frac{1}{2}} and thus the interval
#' \eqn{\left[\frac{1}{2} - \frac{1}{\sqrt{l}}, \frac{1}{2} + \frac{1}{\sqrt{l}}\right]} which is also plotted. So if most of the GC contents are not in the
#' plotted interval a non random GC content can be assumed.
#'
#' @param sequence The sequence to be analysed
#' @param frame_length The length of the individual analysing frames
#' @param plot_resolution The number of how much frequencies should be displayed
gc_plot = function(sequence, frame_length, plot_resolution = 500) {
  x = round(seq(frame_length / 2 , length(sequence) - frame_length / 2 - 1, length.out = plot_resolution))
  FL2 = numeric(plot_resolution)
  for (i in 1:plot_resolution) FL2[i] = GC(sequence[(x[i] - frame_length / 2):(x[i] + frame_length / 2 - 1)])
  plot(x, FL2, type = "l", ylim = c(0,1), xlab = "frame", ylab = "GC-content")
  abline(h = (0.5 + 1 / sqrt(frame_length)), col = "red")
  abline(h = (0.5 - 1 / sqrt(frame_length)), col = "red")
}

#' @title gc_content
#'
#' @description This proceedure operates in the exact same way as \code{\link{gc_plot}} however instead of graphing the results of the analysis it returns its
#' results. It is slower than \code{\link{gc_plot}} since it calculates all relative frequencies. If the deviation result is close to zero the GC distribution
#' is probably random. If it's closer to one the distribution is probably not random.
#'
#' @param sequence The sequence to be analysed
#' @param frame_length The length of the individual analysing frames
#' @param interpretation A switch to determine if a vector with the relative frequencies should be returned or a value of how much the GC contents deviate from
#' a random distribution.
#'
#' @returns Either a vector of relative GC contents or a single number indicating the deviation of a random GC distribution based on the value of
#' \code{interpretation}
gc_content = function(sequence, frame_length, interpretation = TRUE) {

    FL = numeric(length(sequence) - frame_length + 1)
    for(i in 1:(length(sequence) - frame_length + 1)) {
      FL[i] = GC(sequence[i:(i + frame_length- 1)])
    }

    if(!interpretation) {
      return(FL)
    } else {
      return(mean(FL > (0.5 + 1 / sqrt(frame_length)) | FL < (0.5 - 1 / sqrt(frame_length))))
    }
}

#' @title dna_to_aa_sequence
#'
#' @description Utility function to resolve \code{translate} duplication from \code{Biostrings} and \code{seqinr} to apply a genetic code. Let
#' \eqn{\mathscr{C} \subseteq \mathscr{N}^3} and two different symbols \eqn{\text{START, STOP}}. The elements of \eqn{\mathscr{C}} are being called codons. A
#' genetic code \eqn{g} is a surjective function \eqn{g: \mathscr{C} \longrightarrow \mathscr{A} \cup \{\text{START, STOP}\}}. \eqn{\mathscr{N} = \{A, C, G, T\}},
#' and the amino acids \eqn{\mathscr{A} = \{A, R, N, D, C, Q, E, G, H, I, L , K, M, F, P, S, T, W, Y, V\}}.
#'
#' @param sequence Input DNA Sequence
#'
#' @return Translated amino acid sequence
dna_to_aa_sequence = function(sequence) {
  return(seqinr::translate(sequence))
}


#' @title match_orfs
#'
#' @description Calculate all reading frames with a given length. Given a DNA Sequence \eqn{\textbf{s} = \textbf{s}_{[1, 3n]}} with length \eqn{3n}. Then a
#' triplet sequence \eqn{(t_k)_{k \in [1, n]}} from \eqn{\mathscr{N}^3} with \eqn{t_k = \textbf{s}_{[3(k - 1) + 1, k]}} is called a reading frame of \eqn{\textbf{s}}.
#' This output can be accomplished by switching \code{filter_orfs} to \code{FALSE}. For \code{filter_orfs = TRUE} all open reading frames are calculated.
#' Given a genetic code \eqn{g} (\code{\link{dna_to_aa_sequence}}) and a DNA Sequence \eqn{\textbf{s}} with a partial sequence \eqn{\textbf{s}'} with length \eqn{3n}.
#' A reading frame \eqn{(t_k)_{k \in [1, n]}} for \eqn{\textbf{s}'} is called an open reading frame regarfing \eqn{q} if \eqn{\forall k \in [1,n]: t_k \in \mathscr{C}},
#' \eqn{g(t_1) = \text{START}}, \eqn{g(t_n) = \text{STOP}} and \eqn{\forall k \in [1, k - 1]: g(t_k) \neq \text{STOP}}.
#'
#' @param aa_sequence The input amino acid sequence based on \code{AMINO_ALPHABET}
#' @param max_rf_length The maximum length of the reading frames
#' @param filter_orfs Only return open reading frames
#'
#' @returns Either a \code{XStringViews} object of all reading frames or all open reading frames depending on \code{filter_orfs} without the STOP codon.
match_orfs = function(aa_sequence, max_rf_length = 300, filter_orfs = TRUE) {
  AAs = AAString(paste(aa_sequence, collapse = ""))
  LRP = matchLRPatterns("M", "*", max_rf_length, AAs)

  #Filter ORFs from RFs
  valid = rep(TRUE, length(LRP))
  if(filter_orfs){
    for(i in 1:length(LRP)) {

      if(str_count(as.character(LRP[i]), "\\*") > 1) {
        valid[i] = FALSE
      }
    }
  }

  return(LRP[valid] - .5) #remove stop codon to get correct widths
}

#' @title rf_plot
#'
#' @description Plot the frequencies of the width of the input reading frames in a histogram. Assuming the original sequence was a multinomially distributed random
#' DNA Sequence a length can be calculated so that an open reading frame greater than that length is probably not described by a random distribution and thus
#' biologically relevant. This length can be approximated by a quantile with a specified level of confidence. The quantile is also plotted. The probability that
#' an open reading frame with length \eqn{\geq k} is found in a random distribution can be calculated as
#' \eqn{\mathbb{P}(k) = \mathbb{P}(\text{STOP}) \cdot (1 - \mathbb{P}(\text{STOP})^{k-1}) = \frac{3}{64} \cdot \left(1 - \frac{3}{64}\right)^{k - 1}} which
#' can also easily be calculated by \code{1 - ecdf(width(orfs))(k - 1)}. This equals to the p-value of how random the open reading frame actually is.
#'
#'
#' @param rfs Input reading frames/open reading frames
#' @param sequence Original DNA Sequence to initialize multinomial distribution
#' @param max_rf_length The maximum length of the reading frames
#' @param quantile The level of confidence needed for quantile calculation
#' @param show_orf_quantile Plot the quantile
#'
#' @return The lowest probably (confidence level) not randomly distributed open reading frame length
rf_plot = function(rfs, sequence, max_rf_length = 300, quantile = .95, show_orf_quantile = TRUE) {
  random_dna = multinomial_distribution_model(DNA_ALPHABET, 200000, initialize_multinomial(sequence, DNA_ALPHABET))
  random_rf = match_orfs(dna_to_aa_sequence(random_dna), max_rf_length = max_rf_length, filter_orfs = FALSE)
  #Save time by only calculating widths of orfs
  random_widths = width(random_rf)[start(random_rf)!=c(0, start(random_rf)[-length(random_rf)])]-1
  q = quantile(random_widths, quantile)

  hist(width(rfs), breaks = 1000)
  abline(v = q, col = "red")
  return(q)
}


#' @title read_dna_sequence
#'
#' @description Read a file of specified filetype containing DNA information and return the DNA sequence as a vector of characters.
#'
#' @param path The path to the file to be read
#' @param filetype The type of the file to be read currently supported are the values \code{"fasta"} for .fasta files and \code{"dat"} for .dat files
#'
#' @returns The DNA sequence as a vector of characters
read_dna_sequence = function(path, filetype){
  res = NULL
  if(filetype == "fasta"){
    res = read.fasta(path)
    res = res[[1]]
    res = as.character(res)
  }

  if(filetype == "dat"){
    res = scan(path, what="character")
  }

  return(toUpper(res))
}
