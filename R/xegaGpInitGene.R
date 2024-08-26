
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
#'          package \code{xegaGpGene}, a \emph{gene} is a list with 
#'          \enumerate{
#'          \item \code{$gene1}:     a derivation tree.
#'          \item \code{$fit}:       The fitness of the genotype of 
#'                                  \code{$gene1}         
#'          \item \code{$evaluated}: Boolean: TRUE if the fitness is known.
#'          \item \code{$evalFail}:   Has the evaluation of the gene failed?
#'    \item \code{$var}:        The cumulative variance of the fitness 
#'                      of all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item \code{$sigma}:      The standard deviation of the fitness of 
#'                      all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item \code{$obs}:        The number of evaluations of a gene.
#'                      (For stochastic functions)
#'          }
#'
#'          The algorithm for generating a complete derivation tree 
#'          with a depth-bound
#'          is imported from the package \code{xegaDerivationTrees}. 
#' 
#' @param lF  Local configuration of the genetic algorithm.
#'
#' @return Derivation tree.
#'
#' @family Gene Generation
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

