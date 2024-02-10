
#
# (c) 2021 Andreas Geyer-Schulz
#     Grammar-based Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            Gene operations with derivation trees.
#     Package: xegaGpGene
#

#' Crossover of 2 derivation tree genes with node filter.
#'
#' @description \code{xegaGpFilterCross2Gene()} swaps two randomly extracted 
#'              subtrees between 2 genes. Subtrees must have the same
#'              root in order to be compatible. The current implementation 
#'              performs at most \code{maxtrials} trials to find compatible
#'              subtrees. If this fails, the original genes are returned.
#'              Only nodes with a depth
#'              between \code{lF$MinMutInsertionDepth()} and
#'              \code{lF$MaxMutInsertionDepth()} are considered as
#'              candidate roots of derivation trees to be swapped 
#'              by crossover.
#'
#' @details Crossover is controlled by three local parameters:
#'   \itemize{
#'      \item \code{lF$MinCrossDepth()} and 
#'       \code{lF$MaxCrossDepth()} control the possible exchange points 
#'           for subtrees. The depth of the exchange node must be 
#'                  between \code{lF$MinMutInsertionDepth()} and
#'                  \code{lF$MaxMutInsertionDepth()}.
#'      \item \code{lF$MaxTrials()}: Maximal number of trials to find 
#'                compatible subtrees. If compatible subtrees are not 
#'                found, the gene is returned unchanged.
#'          }
#'
#' @param ng1          Derivation tree.
#' @param ng2          Derivation tree.
#' @param lF           Local configuration of the genetic algorithm.
#'
#' @return List of 2 derivation trees.
#'
#' @family Crossover
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene2<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene2, lFxegaGpGene)
#' newgenes<-xegaGpFilterCross2Gene(gene1, gene2,  lFxegaGpGene)
#' xegaGpDecodeGene(newgenes[[1]], lFxegaGpGene)
#' xegaGpDecodeGene(newgenes[[2]], lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANL
#' @importFrom xegaDerivationTrees filterANLid
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees compatibleSubtrees 
#' @importFrom xegaDerivationTrees treeExtract 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpFilterCross2Gene<-function(ng1, ng2, lF)
{ g1<-ng1$gene1
  g2<-ng2$gene1
  anl1<-xegaDerivationTrees::treeANL(g1, lF$Grammar$ST, lF$MaxDepth())
  anl1<-xegaDerivationTrees::filterANL(anl1,
                 minb=lF$MinCrossDepth(),
                 maxb=lF$MaxCrossDepth())
  anl2<-xegaDerivationTrees::treeANL(g2, lF$Grammar$ST, lF$MaxDepth())
  anl2<-xegaDerivationTrees::filterANL(anl2,
                 minb=lF$MinCrossDepth(),
                 maxb=lF$MaxCrossDepth())
  rg<-list()
# TODO: Replace the maxtrials loop ...
for (i in 1: lF$MaxTrials())
  { n1<-xegaDerivationTrees::chooseNode(anl1$ANL)
    ANL2<-filterANLid(anl2, n1$ID)
    if (length(ANL2$ANL)==0) {next}
    n2<-xegaDerivationTrees::chooseNode(ANL2$ANL)
    if (xegaDerivationTrees::compatibleSubtrees(n1, n2, lF$MaxDepth())) {
    subtree1<-xegaDerivationTrees::treeExtract(g1, n1)
    subtree2<-xegaDerivationTrees::treeExtract(g2, n2)
    newg1<-xegaDerivationTrees::treeInsert(g1, subtree2, n1)
    newg2<-xegaDerivationTrees::treeInsert(g2, subtree1, n2)
# print("cross over SUCESS.")
    rg[[1]]<-list(evaluated=FALSE, fit=0, gene1=newg1)
    rg[[2]]<-list(evaluated=FALSE, fit=0, gene1=newg2)
    return(rg)}}
# print("cross over fails. Return genes.")
rg[[1]]<-ng1
rg[[2]]<-ng2
return(rg)
}

#' Crossover of 2 derivation tree genes
#'
#' @description \code{xegaGpAllCross2Gene()} swaps two randomly extracted 
#'              subtrees between 2 genes. Subtrees must have the same
#'              root in order to be compatible. The current implementation 
#'              performs at most \code{maxtrials} trials to find compatible
#'              subtrees. If this fails, the original genes are returned.
#'
#' @details Crossover is controlled by one local parameter:
#'   \itemize{
#'      \item \code{lF$MaxTrials()}: Maximal number of trials to find 
#'                compatible subtrees. If compatible subtrees are not 
#'                found, the gene is returned unchanged.
#'          }
#'
#' @param ng1          Derivation tree.
#' @param ng2          Derivation tree.
#' @param lF           Local configuration of the genetic algorithm.
#'
#' @return List of 2 derivation trees.
#'
#' @family Crossover
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene2<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene2, lFxegaGpGene)
#' newgenes<-xegaGpAllCross2Gene(gene1, gene2,  lFxegaGpGene)
#' xegaGpDecodeGene(newgenes[[1]], lFxegaGpGene)
#' xegaGpDecodeGene(newgenes[[2]], lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANLid
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees compatibleSubtrees 
#' @importFrom xegaDerivationTrees treeExtract 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpAllCross2Gene<-function(ng1, ng2, lF)
{ g1<-ng1$gene1
  g2<-ng2$gene1
  anl1<-xegaDerivationTrees::treeANL(g1, lF$Grammar$ST, lF$MaxCrossDepth())
  anl2<-xegaDerivationTrees::treeANL(g2, lF$Grammar$ST, lF$MaxCrossDepth())
  rg<-list()
# TODO: Replace the maxtrials loop ...
for (i in 1: lF$MaxTrials())
  { n1<-xegaDerivationTrees::chooseNode(anl1$ANL)
    ANL2<-filterANLid(anl2, n1$ID)
    if (length(ANL2$ANL)==0) {next}
    n2<-xegaDerivationTrees::chooseNode(ANL2$ANL)
    if (xegaDerivationTrees::compatibleSubtrees(n1, n2, lF$MaxDepth())) {
    subtree1<-xegaDerivationTrees::treeExtract(g1, n1)
    subtree2<-xegaDerivationTrees::treeExtract(g2, n2)
    newg1<-xegaDerivationTrees::treeInsert(g1, subtree2, n1)
    newg2<-xegaDerivationTrees::treeInsert(g2, subtree1, n2)
# print("cross over SUCESS.")
    rg[[1]]<-list(evaluated=FALSE, fit=0, gene1=newg1)
    rg[[2]]<-list(evaluated=FALSE, fit=0, gene1=newg2)
    return(rg)}}
# print("cross over fails. Return genes.")
rg[[1]]<-ng1
rg[[2]]<-ng2
return(rg)
}

#' Crossover of 2 derivation tree genes.
#'
#' @description \code{xegaGpAllCrossGene()} swaps two randomly extracted 
#'              subtrees between 2 genes. Subtrees must have the same
#'              root in order to be compatible. The current implementation 
#'              performs at most \code{lF$MaxTrials()} 
#'              attempts to find compatible
#'              subtrees. If this fails, the original gene is returned.
#'
#' @details Crossover is controlled by one local parameter:
#'   \itemize{
#'      \item \code{lF$MaxTrials()}: Maximal number of trials to find 
#'                compatible subtrees. If compatible subtrees are not 
#'                found, the gene is returned unchanged.
#'          }
#'
#' @param ng1          Derivation tree.
#' @param ng2          Derivation tree.
#' @param lF           Local configuration of the genetic algorithm.
#'
#' @return List of 1 derivation tree.
#'
#' @family Crossover
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene2<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene2, lFxegaGpGene)
#' newgene<-xegaGpAllCrossGene(gene1, gene2,  lFxegaGpGene)
#' xegaGpDecodeGene(newgene[[1]], lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANLid
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees compatibleSubtrees 
#' @importFrom xegaDerivationTrees treeExtract 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpAllCrossGene<-function(ng1, ng2, lF)
{ g1<-ng1$gene1
  g2<-ng2$gene1
  anl1<-xegaDerivationTrees::treeANL(g1, lF$Grammar$ST, lF$MaxDepth())
  anl2<-xegaDerivationTrees::treeANL(g2, lF$Grammar$ST, lF$MaxDepth())
  rg<-list()
# TODO: Replace the maxtrials loop ...
for (i in 1: lF$MaxTrials())
  { n1<-xegaDerivationTrees::chooseNode(anl1$ANL)
    ANL2<-filterANLid(anl2, n1$ID)
    if (length(ANL2$ANL)==0) {next}
    n2<-xegaDerivationTrees::chooseNode(ANL2$ANL)
    if (xegaDerivationTrees::compatibleSubtrees(n1, n2, lF$MaxDepth())) {
    subtree2<-xegaDerivationTrees::treeExtract(g2, n2)
    newg1<-xegaDerivationTrees::treeInsert(g1, subtree2, n1)
# print("cross over SUCESS.")
    rg[[1]]<-list(evaluated=FALSE, fit=0, gene1=newg1)
    return(rg)}}
# print("cross over fails. Return genes.")
rg[[1]]<-ng1
return(rg)
}

#' Crossover of 2 derivation tree genes with node filter.
#'
#' @description \code{xegaGpFilterCrossGene()} swaps two randomly extracted 
#'              subtrees between 2 genes. Subtrees must have the same
#'              root in order to be compatible. The current implementation 
#'              performs at most \code{lF$maxtrials()} 
#'              attempts to find compatible
#'              subtrees. If this fails, the original gene is returned.
#'              Only nodes with a depth
#'              between \code{lF$MinMutInsertionDepth()} and
#'              \code{lF$MaxMutInsertionDepth()} are considered as
#'              candidate roots of derivation trees to be swapped 
#'              by crossover.
#'
#' @details Crossover is controlled by three local parameters:
#'   \itemize{
#'      \item \code{lF$MinCrossDepth()} and 
#'       \code{lF$MaxCrossDepth()} control the possible exchange points 
#'           for subtrees. The depth of the exchange node must be 
#'                  between \code{lF$MinMutInsertionDepth()} and
#'                  \code{lF$MaxMutInsertionDepth()}.
#'      \item \code{lF$MaxTrials()}: Maximal number of trials to find 
#'                compatible subtrees. If compatible subtrees are not 
#'                found, the gene is returned unchanged.
#'          }
#'
#' @param ng1          Derivation tree.
#' @param ng2          Derivation tree.
#' @param lF           Local configuration of the genetic algorithm.
#'
#' @return List of 1 derivation tree.
#'
#' @family Crossover
#'
#' @examples
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene2<-xegaGpInitGene(lFxegaGpGene)
#' xegaGpDecodeGene(gene1, lFxegaGpGene)
#' xegaGpDecodeGene(gene2, lFxegaGpGene)
#' newgene<-xegaGpFilterCrossGene(gene1, gene2,  lFxegaGpGene)
#' xegaGpDecodeGene(newgene[[1]], lFxegaGpGene)
#'
#' @importFrom xegaDerivationTrees treeANL
#' @importFrom xegaDerivationTrees filterANL
#' @importFrom xegaDerivationTrees filterANLid
#' @importFrom xegaDerivationTrees chooseNode 
#' @importFrom xegaDerivationTrees compatibleSubtrees 
#' @importFrom xegaDerivationTrees treeExtract 
#' @importFrom xegaDerivationTrees treeInsert 
#' @export
xegaGpFilterCrossGene<-function(ng1, ng2, lF)
{ g1<-ng1$gene1
  g2<-ng2$gene1
  anl1<-xegaDerivationTrees::treeANL(g1, lF$Grammar$ST, lF$MaxDepth())
  anl1<-xegaDerivationTrees::filterANL(anl1,
                 minb=lF$MinCrossDepth(),
                 maxb=lF$MaxCrossDepth())
  anl2<-xegaDerivationTrees::treeANL(g2, lF$Grammar$ST, lF$MaxDepth())
  anl2<-xegaDerivationTrees::filterANL(anl2,
                 minb=lF$MinCrossDepth(),
                 maxb=lF$MaxCrossDepth())
  rg<-list()
# TODO: Replace the maxtrials loop ...
for (i in 1: lF$MaxTrials())
  { n1<-xegaDerivationTrees::chooseNode(anl1$ANL)
    ANL2<-filterANLid(anl2, n1$ID)
    if (length(ANL2$ANL)==0) {next}
    n2<-xegaDerivationTrees::chooseNode(ANL2$ANL)
    if (xegaDerivationTrees::compatibleSubtrees(n1, n2, lF$MaxDepth())) {
    subtree2<-xegaDerivationTrees::treeExtract(g2, n2)
    newg1<-xegaDerivationTrees::treeInsert(g1, subtree2, n1)
# print("cross over SUCESS.")
    rg[[1]]<-list(evaluated=FALSE, fit=0, gene1=newg1)
    return(rg)}}
# print("cross over fails. Return genes.")
rg[[1]]<-ng1
return(rg)
}

#' Configure the crossover function of a grammar-based genetic algorithm.
#'
#' @description \code{xegaGpCrossoverFactory()} implements the selection
#'              of one of the crossover functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item Crossover functions with two kids:
#'              \enumerate{
#'              \item "Cross2Gene"       returns \code{xegaGpAllCross2Gene()}.
#'              \item "AllCross2Gene"    returns \code{xegaGpAllCross2Gene()}.
#'              \item "FilterCross2Gene" returns \code{xegaGpFilterCross2Gene()}.
#'              }
#'              \item Crossover functions with one kid:
#'              \enumerate{
#'              \item "AllCrossGene" returns \code{xegaGpAllCrossGene()}.
#'              \item "FilterCrossGene" returns \code{xegaGpFilterCrossGene()}.
#'              }
#'              }
#'
#' @param method    String specifying the crossover function.
#'
#' @return Crossover function for genes.
#'
#' @family Configuration
#'
#' @examples
#' XGeneTwo<-xegaGpCrossoverFactory("Cross2Gene")
#' XGeneOne<-xegaGpCrossoverFactory("FilterCrossGene")
#' gene1<-xegaGpInitGene(lFxegaGpGene)
#' gene2<-xegaGpInitGene(lFxegaGpGene)
#' XGeneTwo(gene1, gene2, lFxegaGpGene)
#' XGeneOne(gene1, gene2, lFxegaGpGene)
#' @export
xegaGpCrossoverFactory<-function(method="Cross2Gene") {
if (method=="Cross2Gene") {f<- xegaGpAllCross2Gene}
if (method=="AllCross2Gene") {f<- xegaGpAllCross2Gene}
if (method=="AllCrossGene") {f<- xegaGpAllCrossGene}
if (method=="FilterCross2Gene") {f<- xegaGpFilterCross2Gene}
if (method=="FilterCrossGene") {f<- xegaGpFilterCrossGene}
if (!exists("f", inherits=FALSE))
        {stop("sgp Crossover label ", method, " does not exist")}
return(f)
}

