
#' Genetic operations for grammar-based genetic algorithms.
#' 
#' For derivation tree genes, the \code{xegaGpGene} package provides
#' \itemize{
#' \item Gene initiatilization.
#' \item Decoding of parameters.
#' \item Mutation functions as well as a function factory for configuration.
#' \item Crossover functions as well as a function factory for configuration.
#'       We provide two families of crossover functions:
#'  \enumerate{
#' \item Crossover functions with two kids:
#'       Crossover preserves the genetic information in the gene pool.
#' \item Crossover functions with one kid:
#'       These functions allow the construction of evaluation pipelines
#'       for genes. One advantage of this is a simple control structure 
#'       at the population level.
#' }
#' }
#'
#' @references  Geyer-Schulz, Andreas (1997):
#'          \emph{Fuzzy Rule-Based Expert Systems and Genetic Machine Learning},
#'                Physica, Heidelberg.
#'           (ISBN:978-3-7908-0830-X)
#'
#' @section Derivation Tree Gene Representation:
#'            
#' A derivation tree gene is a named list:
#'   \itemize{
#'    \item \code{$gene1}:     The gene must be a complete derivation tree.
#'    \item \code{$fit}:       The fitness value of the gene
#'                             (for EvalGeneDet() and EvalGeneU()) or
#'                             the mean fitness (for stochastic functions
#'                             evaluated with EvalGeneStoch()).
#'    \item \code{$evaluated}: Boolean. Has the gene been evaluated?
#'    \item \code{$evalFail}:  Boolean. Has the evaluation of the gene failed?
#'    \item \code{$var}:       The variance of the fitness 
#'                             of all evaluations of a gene is updated
#'                             after each evaluation of a gene.
#'                             (For stochastic functions.)
#'    \item \code{$sigma}:     The standard deviation of the fitness of 
#'                             all evaluations of a gene.
#'                             (For stochastic functions.)
#'    \item \code{$obs:}       The number evaluations of a gene.
#'                             (For stochastic functions.)
#'   }
#'
#' @section Abstract Interface of Problem Environment:
#'
#' A problem environment \code{penv} must provide:
#'   \itemize{
#'     \item \code{$f(word, gene, lF)}: 
#'   Function with a word of a language as first argument
#'   which returns a gene 
#'   with evaluated fitness.
#'   
#'   } 
#'
#' @section Abstract Interface of Mutation Functions:
#'
#' Each mutation function has the following function signature:
#'
#'     \code{newGene<-Mutate(gene, lF)}
#'
#' All local parameters of the mutation function configured are 
#' expected in the local function list lF.
#' 
#' @section Local Constants of Mutation Functions:
#'
#' The local constants of a mutation function determine the 
#' the behavior of the function. 
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' \code{lF$MaxMutDepth()} \tab 3  \tab xegaGpMutateGene \cr 
#' \code{lF$MinMutInsertiontDepth()} \tab 3  \tab xegaGpMutateGene \cr 
#' \code{lF$MaxMutInsertiontDepth()} \tab 4  \tab xegaGpMutateGene \cr 
#' }
#'
#' @section Abstract Interface of Crossover Functions:
#'
#' The signatures of the abstract interface to the 2 families 
#' of crossover functions are:
#'
#'     \code{ListOfTwoGenes<-Crossover2(gene1, gene2, lF)} 
#'
#'     \code{newGene<-Crossover(gene1, gene2, lF)}
#'
#' All local parameters of the crossover function configured are 
#' expected in the local function list lF.
#'
#' @section Local Constants of Crossover Functions:
#'
#' \tabular{rcl}{ 
#' \strong{Constant} \tab \strong{Default} \tab \strong{Used in} \cr 
#' \code{lF$MinCrossDepth()} \tab 1  \tab xegaGpCrossGene,  \cr 
#' \code{lF$MaxCrossDepth()} \tab 7  \tab xegaGpCrossGene,  \cr 
#' \code{lF$MaxTrials()}     \tab 5  \tab xegaGpCross2Gene  \cr 
#' }
#'
#' @section The Architecture of the xegaX-Packages:
#' 
#' The xegaX-packages are a family of R-packages which implement 
#' eXtended Evolutionary and Genetic Algorithms (xega).  
#' The architecture has 3 layers, 
#' namely the user interface layer,
#' the population layer, and the gene layer: 
#' 
#' \itemize{
#' \item
#' The user interface layer (package \code{xega}) 
#' provides a function call interface and configuration support
#' for several algorithms: genetic algorithms (sga), 
#' permutation-based genetic algorithms (sgPerm), 
#' derivation free algorithms as e.g. differential evolution (sgde), 
#' grammar-based genetic programming (sgp) and grammatical evolution
#' (sge). 
#'
#' \item
#' The population layer (package \code{xegaPopulation}) contains
#' population related functionality as well as support for 
#' population statistics dependent adaptive mechanisms and parallelization.
#'
#' \item 
#' The gene layer is split in a representation independent and 
#' a representation dependent part:
#' \enumerate{
#' \item 
#'  The representation indendent part (package \code{xegaSelectGene})
#'  is responsible for variants of selection operators, evaluation 
#'  strategies for genes, as well as profiling and timing capabilities.        
#' \item 
#'  The representation dependent part consists of the following packages: 
#' \itemize{
#' \item \code{xegaGaGene} for binary coded genetic algorithms.
#' \item \code{xegaPermGene} for permutation-based genetic algorithms.
#' \item \code{xegaDfGene} for derivation free algorithms as e.g. 
#'                         differential evolution.
#' \item \code{xegaGpGene} for grammar-based genetic algorithms.
#' \item \code{xegaGeGene} for grammatical evolution algorithms.
#' }
#' The packages \code{xegaDerivationTrees} and \code{xegaBNF} support
#' the last two packages:
#' \code{xegaBNF} essentially provides a grammar compiler and 
#' \code{xegaDerivationTrees} an abstract data type for derivation trees.
#' }} 
#'
#' @family Package Description
#'
#' @name xegaGpGene
#' @aliases xegaGpGene
#' @docType package
#' @title Package xegaGpGene.
#' @author Andreas Geyer-Schulz
#' @section Copyright: (c) 2023 Andreas Geyer-Schulz
#' @section License: MIT
#' @section URL: <https://github.com/ageyerschulz/xegaGpGene>
#' @section Installation: From CRAN by \code{install.packages('xegaGpGene')}
NULL

