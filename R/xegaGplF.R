
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
#' We enhance the configurability of our code by introducing 
#'  a function factory. The  function factory contains
#'  all the functions that are needed for defining
#'  local functions in genetic operators. The local function  
#'  list keeps the signatures of functions (e.g. mutation functions)
#'  uniform and small. At the same time, variants of functions
#'  can use different local functions. 
#'
#' @details
#'    We use the local function list for 
#'    \enumerate{
#'    \item
#'       replacing all constants by constant functions.
#'       
#'       Rationale: We need one formal argument (the local function list lF)
#'       and we can dispatch multiple functions. E.g.  \code{lF$verbose()}
#'   \item    
#'       dynamically binding a local function with a definition from a
#'       proper function factory. E.g., the selection methods 
#'       \code{lF$SelectGene()} and \code{lF$SelectMate()}.
#'       
#'  \item gene representations which require special functions to handle them:
#'        For example,
#'        \code{lF$InitGene()}, \code{lF$DecodeGene()}, \code{lF$EvalGene()},
#'        \code{lF$ReplicateGene()}, ...
#'       
#'  } 
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

