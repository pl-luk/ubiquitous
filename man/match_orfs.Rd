% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ubiquitous.R
\name{match_orfs}
\alias{match_orfs}
\title{match_orfs}
\usage{
match_orfs(aa_sequence, max_rf_length = 300, filter_orfs = TRUE)
}
\arguments{
\item{aa_sequence}{The input amino acid sequence based on \code{AMINO_ALPHABET}}

\item{max_rf_length}{The maximum length of the reading frames}

\item{filter_orfs}{Only return open reading frames}
}
\value{
Either a \code{XStringViews} object of all reading frames or all open reading frames depending on \code{filter_orfs} without the STOP codon.
}
\description{
Calculate all reading frames with a given length. Given a DNA Sequence \eqn{\textbf{s} = \textbf{s}_{[1, 3n]}} with length \eqn{3n}. Then a
triplet sequence \eqn{(t_k)_{k \in [1, n]}} from \eqn{\mathscr{N}^3} with \eqn{t_k = \textbf{s}_{[3(k - 1) + 1, k]}} is called a reading frame of \eqn{\textbf{s}}.
This output can be accomplished by switching \code{filter_orfs} to \code{FALSE}. For \code{filter_orfs = TRUE} all open reading frames are calculated.
Given a genetic code \eqn{g} (\code{\link{dna_to_aa_sequence}}) and a DNA Sequence \eqn{\textbf{s}} with a partial sequence \eqn{\textbf{s}'} with length \eqn{3n}.
A reading frame \eqn{(t_k)_{k \in [1, n]}} for \eqn{\textbf{s}'} is called an open reading frame regarfing \eqn{q} if \eqn{\forall k \in [1,n]: t_k \in \mathscr{C}},
\eqn{g(t_1) = \text{START}}, \eqn{g(t_n) = \text{STOP}} and \eqn{\forall k \in [1, k - 1]: g(t_k) \neq \text{STOP}}.
}
