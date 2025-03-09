
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
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees randomDerivationTree
#' @export
xegaGpInitGene<-function(lF)
{ gene1<-xegaDerivationTrees::randomDerivationTree(lF$Grammar$Start, 
						 lF$Grammar, lF$MaxDepth())
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1))
}


#' Generates a gene as a random derivation tree from a random integer vector.
#'
#' @description For a given grammar, \code{xegaGpInitGene()} 
#'              generates a gene as a random derivation tree
#'              with a depth-bound. This function uses almost the same 
#'              initialization algorithm as for grammatical evolution.
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
#' gene<-xegaGpInitGeneGe(lFxegaGpGene)
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees generateCDT
#' @export
xegaGpInitGeneGe<-function(lF)
{ vec<-sample(1000, max(100,2^lF$MaxDepth()), replace=TRUE)
  gene1<-xegaDerivationTrees::generateCDT(sym=lF$Grammar$Start, 
           kvec=vec, complete=TRUE, G=lF$Grammar, maxdepth=lF$MaxDepth())
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1$tree))
}

#' Configure the initialization function of grammar-based genetic programming.
#'
#' @description \code{xegaGpInitGeneFactory()} implements the 
#'              creation of a complete derivation tree by one 
#'              an algorithm  selected 
#'              by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "InitGene" returns \code{xegaGpInitGene()}.
#'              \item "InitGeneGe" returns \code{xegaGpInitGeneGe()}.
#'              }
#'
#' @param method     String specifying the mutation function.
#'
#' @return Initialization function for genes.
#'
#' @family Configuration
#'
#' @examples
#' InitGene<-xegaGpInitGeneFactory("InitGene")
#' gene1<-InitGene(lFxegaGpGene)
#' InitGene<-xegaGpInitGeneFactory("InitGeneGe")
#' gene2<-InitGene(lFxegaGpGene)
#'
#' @export
xegaGpInitGeneFactory<-function(method="InitGene") {
if (method=="InitGene") {f<- xegaGpInitGene}
if (method=="InitGeneGe") {f<- xegaGpInitGeneGe}
if (!exists("f", inherits=FALSE))
        {stop("sgp InitGene label ", method, " does not exist")}
return(f)
}

