
#' Prints a random example of crossover for a crossover method given 
#' a random number seed. 
#'
#' @description The function supports the search 
#'              for examples for unit tests for crossover 
#'              functions whose behavior depends on random numbers.
#'
#' @param  FUN        String. Specification of crossover method.
#' @param  s          Integer. Seed of random number generator. 
#' @param  verbose    Boolean. 
#'                    If \code{TRUE} (default), print the example 
#'                    to the console. 
#'
#' @return No return.
#'
#' @family Testing
#'
#' @examples
#' findCrossoverExample(FUN="AllCross2Gene", s=2)
# findCrossoverExample(FUN="FilterCross2Gene", s=19)
#' @export
findCrossoverExample<-function(FUN, s, verbose=TRUE)
{
set.seed(s)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
CROSSOVER<-xegaGpCrossoverFactory(method=FUN)
gene<-CROSSOVER(gene1, gene2, lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
if (verbose)
     {cat(" g1", a, "\n")
      cat(" g2", c, "\n")}
for (i in (1:length(gene)))
{
b<-xegaGpDecodeGene(gene[[i]], lFxegaGpGene)
if (verbose)  {cat("ng", i, ":", b, "\n")}
}
}

