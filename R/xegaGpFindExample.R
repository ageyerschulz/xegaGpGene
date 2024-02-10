
#' Prints a random example of crossover for a crossover method given 
#' a random number seed. 
#'
#' @param  FUN    String. Specification of crossover method.
#' @param  s      Integer. Seed of random number generator. 
#'
#' @return No return.
#'
#' @family Testing
#'
#' @examples
#' findCrossoverExample(FUN="AllCross2Gene", s=2)
# findCrossoverExample(FUN="FilterCross2Gene", s=19)
#' @export
findCrossoverExample<-function(FUN, s)
{
set.seed(s)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
CROSSOVER<-xegaGpCrossoverFactory(method=FUN)
gene<-CROSSOVER(gene1, gene2, lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
cat(" g1", a, "\n")
cat(" g2", c, "\n")
for (i in (1:length(gene)))
{
b<-xegaGpDecodeGene(gene[[i]], lFxegaGpGene)
cat("ng", i, ":", b, "\n")
}
}

