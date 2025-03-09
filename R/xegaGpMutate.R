
#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            For gene representation of derivation trees.
#     Package: xegaGpGene
#

#' Mutate a gene.
#'
#' @description \code{xegaGpMutateAllGene()} 
#'               replaces a randomly selected subtree by
#'               a random derivation tree with the same root symbol 
#'               with a small probability.
#'               All non-terminal nodes are considered as insertion points.
#'               Depth-bounds are respected.
#'
#' @details  Mutation is controlled by one local parameter: 
#'           \enumerate{
#'            \item \code{lF$MaxMutDepth()} controls the maximal depth of 
#'                  the new random generation tree.
#'           }
#'           This version of the genetic operator skips the filter loop.
#'
#' @param g        Derivation tree.
#' @param lF       Local configuration of the genetic algorithm.
#'
#' @return Derivation tree.
#'
#' @family Mutation
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' gene<-xegaGpMutateAllGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom stats runif
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees randomDerivationTree 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpMutateAllGene<-function(g, lF)
{ gene<-g$gene1
  anl<-xegaDerivationTrees::treeANL(gene, ST=lF$Grammar$ST, 
				  maxdepth=lF$MaxDepth())
  node<-xegaDerivationTrees::chooseNode(anl$ANL)	
  mutgene<-xegaDerivationTrees::randomDerivationTree(
    node$ID, lF$Grammar, min(lF$MaxMutDepth(),node$Rdepth))
  newgene<-xegaDerivationTrees::treeInsert(gene, mutgene, node)
  a<-newgene
  return(list(evaluated=FALSE, fit=0, gene1=newgene)) }

#' Mutate a gene (with a node filter)
#'
#' @description \code{xegaGpMutateGeneFilter()} replaces 
#'              a randomly selected subtree by
#'              a random derivation tree with the same root symbol 
#'              with a small probability.
#'              Only non-terminal nodes with a depth
#'              between \code{lF$MinMutInsertionDepth()} and
#'              \code{lF$MaxMutInsertionDepth()} are considered 
#'              as tree insertion points.
#'              Depth-bounds are respected.
#'
#' @details  Mutation is controlled by three local parameters: 
#'           \enumerate{
#'            \item \code{lF$MaxMutDepth()} controls the maximal depth of 
#'                  the new random generation tree.
#'            \item \code{lF$MinMutInsertionDepth()} and 
#'                  \code{lF$MaxMutInsertionDepth()} control the possible 
#'                  insertion points for the new random derivation tree.
#'                  The depth of the insertion node must be 
#'                  between \code{lF$MinMutInsertionDepth()} and
#'                  \code{lF$MaxMutInsertionDepth()}.
#'           }
#'
#' @param g        Derivation tree.
#' @param lF       Local configuration of the genetic algorithm.
#'
#' @return Derivation tree.
#'
#' @family Mutation
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' gene<-xegaGpMutateFilterGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene, lFxegaGpGene)
#'
#' @importFrom stats runif
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANL
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees randomDerivationTree 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpMutateFilterGene<-function(g, lF)
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
  a<-newgene
  return(list(evaluated=FALSE, fit=0, gene1=newgene)) }

#' Configure the mutation function of a genetic algorithm.
#'
#' @description \code{xegaGpMutationFactory()} implements the selection
#'              of one of the mutation functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "MutateGene" returns \code{xegaGpMutateAllGene()}.
#'              \item "MutateAllGene" returns \code{xegaGpMutateAllGene()}.
#'              \item "MutateFilterGene" returns \code{xegaGpMutateFilterGene()}.
#'              }
#'
#' @param method     String specifying the mutation function.
#'
#' @return Mutation function for genes.
#'
#' @family Configuration
#'
#' @examples
#' Mutate<-xegaGpMutationFactory("MutateGene")
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene1
#' Mutate(gene1, lFxegaGpGene)
#' @export
xegaGpMutationFactory<-function(method="MutateGene") {
if (method=="MutateGene") {f<- xegaGpMutateAllGene}
if (method=="MutateAllGene") {f<- xegaGpMutateAllGene}
if (method=="MutateFilterGene") {f<- xegaGpMutateFilterGene}
if (!exists("f", inherits=FALSE))
        {stop("sgp Mutation label ", method, " does not exist")}
return(f)
}
