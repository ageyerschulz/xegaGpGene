
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpMutateAllGene OK", 
{
 set.seed(21)
gene1<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpMutateAllGene(gene1, lFxegaGpGene)
b<-xegaGpDecodeGene(gene, lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpMutateFilterGene OK", 
{
 set.seed(21)
gene1<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpMutateFilterGene(gene1, lFxegaGpGene)
b<-xegaGpDecodeGene(gene, lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpMutationFactory MutateGene OK",
 {
 f<-xegaGpMutationFactory(method="MutateGene")
 expect_identical(body(f), body(xegaGpMutateAllGene))
}
)

test_that("xegaGpMutationFactory MutateAllGene OK",
 {
 f<-xegaGpMutationFactory(method="MutateAllGene")
 expect_identical(body(f), body(xegaGpMutateAllGene))
}
)

test_that("xegaGpMutationFactory MutateFilterGene OK",
 {
 f<-xegaGpMutationFactory(method="MutateFilterGene")
 expect_identical(body(f), body(xegaGpMutateFilterGene))
}
)

test_that("xegaGpMutationFactory sgunknown OK",
 {
 expect_error(
 xegaGpMutationFactory(method="sgunknown"),
 "sgunknown")
}
)

