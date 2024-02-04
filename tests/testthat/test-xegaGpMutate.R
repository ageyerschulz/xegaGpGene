
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpMutateGene OK", 
{
 set.seed(21)
gene1<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpMutateGene(gene1, lFxegaGpGene)
b<-xegaGpDecodeGene(gene, lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpMutationFactory MutateGene OK",
 {
 f<-xegaGpMutationFactory(method="MutateGene")
 expect_identical(body(f), body(xegaGpMutateGene))
}
)


test_that("xegaGpMutationFactory sgunknown OK",
 {
 expect_error(
 xegaGpMutationFactory(method="sgunknown"),
 "sgunknown")
}
)

