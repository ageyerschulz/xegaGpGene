
library(testthat)
library(xegaGpGene)

test_that("xegaGpInitGene OK", 
{
set.seed(21)
gene1<-xegaGpInitGene(lFxegaGpGene)
set.seed(2)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
b<-xegaGpDecodeGene(gene2, lFxegaGpGene)
expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpInitGeneGe OK", 
{
set.seed(20)
agene<-xegaGpInitGene(lFxegaGpGene)
set.seed(2)
bgene<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(agene, lFxegaGpGene)
b<-xegaGpDecodeGene(bgene, lFxegaGpGene)
expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpInitGeneFactory InitGene OK",
 {
 f<-xegaGpInitGeneFactory(method="InitGene")
 expect_identical(body(f), body(xegaGpInitGene))
}
)

test_that("xegaGpInitGeneFactory InitGeneGe OK",
 {
 f<-xegaGpInitGeneFactory(method="InitGeneGe")
 expect_identical(body(f), body(xegaGpInitGeneGe))
}
)

test_that("xegaGpInitGeneFactory sgunknown OK",
 {
 expect_error(
 xegaGpInitGeneFactory(method="sgunknown"),
 "sgunknown")
}
)

