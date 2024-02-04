
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpCross2Gene OK", 
{
 set.seed(17)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpCross2Gene(gene1, gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
d<-xegaGpDecodeGene(gene[[2]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
 expect_identical(identical(c, b), FALSE)
 expect_identical(identical(a, d), FALSE)
 expect_identical(identical(c, d), FALSE)
}
)

test_that("xegaGpCrossGene OK", 
{
 set.seed(17)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpCrossGene(gene1, gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpCrossoverFactory CrossGene OK",
 {
 f<-xegaGpCrossoverFactory(method="CrossGene")
 expect_identical(body(f), body(xegaGpCrossGene))
}
)

test_that("xegaGpCrossoverFactory Cross2Gene OK",
 {
 f<-xegaGpCrossoverFactory(method="Cross2Gene")
 expect_identical(body(f), body(xegaGpCross2Gene))
}
)

test_that("xegaGpCrossoverFactory sgunknown OK",
 {
 expect_error(
 xegaGpCrossoverFactory(method="sgunknown"),
 "sgunknown")
}
)

