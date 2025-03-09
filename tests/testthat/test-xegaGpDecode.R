
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpDecodeGene OK", 
{
 set.seed(20)
 gene<-xegaGpInitGene(lFxegaGpGene)
 out<- xegaGpDecodeGene(gene, lFxegaGpGene)
 set.seed(20)
 gene1<-xegaGpInitGene(lFxegaGpGene)
 out1<- xegaGpDecodeGene(gene1, lFxegaGpGene)
 expect_identical(out, out1)
}
)

