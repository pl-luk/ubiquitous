% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ubiquitous.R
\name{dna_to_aa_sequence}
\alias{dna_to_aa_sequence}
\title{dna_to_aa_sequence}
\usage{
dna_to_aa_sequence(sequence)
}
\arguments{
\item{sequence}{Input DNA Sequence}
}
\value{
Translated amino acid sequence
}
\description{
Utility function to resolve \code{translate} duplication from \code{Biostrings} and \code{seqinr} to apply a genetic code. Let
\eqn{\mathscr{C} \subseteq \mathscr{N}^3} and two different symbols \eqn{\text{START, STOP}}. The elements of \eqn{\mathscr{C}} are being called codons. A
genetic code \eqn{g} is a surjective function \eqn{g: \mathscr{C} \longrightarrow \mathscr{A} \cup \{\text{START, STOP}\}}. \eqn{\mathscr{N} = \{A, C, G, T\}},
and the amino acids \eqn{\mathscr{A} = \{A, R, N, D, C, Q, E, G, H, I, L , K, M, F, P, S, T, W, Y, V\}}.
}
