
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpDecodeGene OK", 
{
 set.seed(20)
 gene<-xegaGpInitGene(lFxegaGpGene)
 out<- xegaGpDecodeGene(gene, lFxegaGpGene)
 expect_identical(out, "NOT(AND(NOT(NOT(NOT(D1))),OR(OR(D1,NOT(NOT(D1))),D1)))")
}
)

