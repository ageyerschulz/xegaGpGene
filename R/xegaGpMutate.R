
#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            For gene representation of derivation trees.
#     Package: xegaGpGene
#

#' Mutate a gene.
#'
#' @description \code{xegaGpMutateGene} replaces a randomly selected subtree by
#'               a random derivation tree with the same root symbol 
#'               with small probability.
#'               Depth-bounds are respected.
#'
#' @details  Mutation is controlled by three local parameters: 
#'           \enumerate{
#'            \item \code{lF$MaxMutDepth()} controls the maximal depth of 
#'                  of the new random generation tree.
#'            \item \code{lF$MinMutInsertionDepth()} and 
#'                  \code{lF$MaxMutInsertionDepth()} control the possible 
#'                  insertion points for the new random derivation tree.
#'                  The depth of the insertion node must be 
#'                  between \code{lF$MinMutInsertionDepth()} and
#'                  \code{lF$MaxMutInsertionDepth()}.
#'           }
#' @param g        A derivation tree.
#' @param lF       Local configuration of the genetic algorithm
#'
#' @return A derivation tree.
#'
#' @family Mutation
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' gene<-xegaGpMutateGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom stats runif
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANL
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees randomDerivationTree 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpMutateGene<-function(g, lF)
{ gene<-g$gene1
  anl<-xegaDerivationTrees::treeANL(gene, ST=lF$Grammar$ST, 
				  maxdepth=lF$MaxDepth())
  anl<-xegaDerivationTrees::filterANL(anl, 
                 minb=lF$MinMutInsertionDepth(),
                 maxb=lF$MaxMutInsertionDepth())
  node<-xegaDerivationTrees::chooseNode(anl$ANL)	
  mutgene<-xegaDerivationTrees::randomDerivationTree(
    node$ID, lF$Grammar, min(lF$MaxMutDepth(),node$Rdepth))
  newgene<-xegaDerivationTrees::treeInsert(gene, mutgene, node)
  return(list(evaluated=FALSE, fit=0, gene1=newgene)) }

#' Configure the mutation function of a genetic algorithm.
#'
#' @description \code{xegaGpMutationFactory} implements the selection
#'              of one of the mutation functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "MutateGene" returns \code{MutateGene}.
#'              }
#'
#' @param method A string specifying the mutation function.
#'
#' @return A mutation function for genes.
#'
#' @family Configuration
#'
#' @family configuration
#'
#' @examples
#' Mutate<-xegaGpMutationFactory("MutateGene")
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene1
#' Mutate(gene1, lFxegaGpGene)
#' @export
xegaGpMutationFactory<-function(method="MutateGene") {
if (method=="MutateGene") {f<- xegaGpMutateGene}
if (!exists("f", inherits=FALSE))
        {stop("sgp Mutation label ", method, " does not exist")}
return(f)
}
