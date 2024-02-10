
#
# (c) 2021 Andreas Geyer-Schulz
#     Simple Genetic Programming in R. V0.1
#     Layer: Gene-Level Functions
#            For gene representation of derivation trees.
#     Package: xegaGpGene
#

#' Generate local functions and objects.
#'
#' @description
#' \code{lFxegaPermGene} is a list of functions 
#' which contains a definition of all local objects 
#' required for the use of genetic operators with the 
#  permutation representation. 
#' We refer to this object as local configuration.
#'
#' @importFrom xegaBNF compileBNF
#' @importFrom xegaBNF booleanGrammar
#' @importFrom xegaSelectGene parm
#' @importFrom xegaSelectGene envXOR
#' @export
lFxegaGpGene<-list(
penv=xegaSelectGene::envXOR,
Grammar=xegaBNF::compileBNF(xegaBNF::booleanGrammar()),
replay=xegaSelectGene::parm(0),
verbose=xegaSelectGene::parm(4),
MaxDepth=xegaSelectGene::parm(5),
MaxMutDepth=xegaSelectGene::parm(3),
MinMutInsertionDepth=xegaSelectGene::parm(1),
MaxMutInsertionDepth=xegaSelectGene::parm(5),
MinCrossDepth=xegaSelectGene::parm(1),
MaxCrossDepth=xegaSelectGene::parm(5),
MaxTrials=xegaSelectGene::parm(5),
MutationRate=xegaSelectGene::parm(0.05),
CrossRate=xegaSelectGene::parm(0.2),
Max=xegaSelectGene::parm(1), 
Offset=xegaSelectGene::parm(1),
Eps=xegaSelectGene::parm(0.01),
Elitist=xegaSelectGene::parm(TRUE),
TournamentSize=xegaSelectGene::parm(2),
SelectGene=xegaSelectGene::SelectGeneFactory(method="PropFitDiff"),
SelectMate=xegaSelectGene::SelectGeneFactory(method="Uniform"),
InitGene=xegaGpInitGene,
MutateGene=xegaGpMutateAllGene,
CrossGene=xegaGpFilterCross2Gene,
DecodeGene=xegaGpDecodeGene,
EvalGene=xegaSelectGene::EvalGeneFactory(method="EvalGeneU")
)

