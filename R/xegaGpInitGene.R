
#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            For gene representation of derivation trees.
#     Package: xegaGpGene
#

#' Generates a gene as a random derivation tree.
#'
#' @description For a given grammar, \code{xegaGpInitGene()} 
#'              generates a gene as a random derivation tree
#'              with a depth-bound.
#'
#' @details In the derivation tree representation of 
#'          package \code{xegaGp}, \emph{gene} is a list with 
#'          \enumerate{
#'          \item \code{$evaluated}: Boolean: TRUE if the fitness is known.
#'          \item \code{$fit}:       The fitness of the genotype of 
#'                                  \code{$gene1}         
#'          \item \code{$gene1}:     a derivation tree.
#'          }
#'
#'          This representation makes implementation of several 
#'          code optimizations and generalizations easier.
#'
#'          The algorithm for generating a complete derivation tree 
#'          with a depth-bound
#'          is imported from package \code{xegaDerivationTrees}. 
#' 
#' @param lF  Local configuration of the genetic algorithm.
#'
#' @return Derivation tree.
#'
#' @family Initialization
#'
#' @examples
#' gene<-xegaGpInitGene(lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees randomDerivationTree
#' @export
xegaGpInitGene<-function(lF)
{ gene1<-xegaDerivationTrees::randomDerivationTree(lF$Grammar$Start, 
						 lF$Grammar, lF$MaxDepth())
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1))
}

