
#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            For gene representation of derivation trees.
#     Package: xegaGpGene
#

#' Decode a derivation tree.
#' 
#' @description \code{xegaGpDecodeGene} decodes a derivation tree.
#'
#' @details The recursive algorithm for the decoder is imported 
#'          from package \code{xegaDerivationTrees}.
#'
#' @param gene      A derivation tree.
#' @param lF        The local configuration of the genetic algorithm.
#' 
#' @return The decoded gene.
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees decodeCDT
#' @export
xegaGpDecodeGene<-function(gene, lF)
{ 
	xegaDerivationTrees::decodeCDT(gene$gene1, lF$Grammar$ST)
}

