
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("xegaGpAllCross2Gene OK", 
{
 set.seed(1)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
gene<-xegaGpAllCross2Gene(gene1, gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
d<-xegaGpDecodeGene(gene[[2]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
 expect_identical(identical(c, b), FALSE)
 expect_identical(identical(a, d), FALSE)
 expect_identical(identical(c, d), FALSE)
}
)

test_that("xegaGpFilterCross2Gene OK", 
{
 set.seed(13)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
gene<-xegaGpFilterCross2Gene(gene1, gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
d<-xegaGpDecodeGene(gene[[2]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
 expect_identical(identical(c, b), FALSE)
 expect_identical(identical(a, d), FALSE)
 expect_identical(identical(c, d), FALSE)
}
)

test_that("xegaGpAllCrossGene OK", 
{
 set.seed(16)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
gene<-xegaGpAllCrossGene(gene1, gene2, lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
 expect_identical(identical(a, b), TRUE)
}
)

test_that("xegaGpAllCrossGene OK", 
{
 set.seed(17)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
gene<-xegaGpAllCrossGene(gene1, gene2, lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpAllCrossGene OK", 
{
set.seed(18)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
gene<-xegaGpAllCrossGene(gene1, gene2, lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
c<-xegaGpDecodeGene(gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
expect_identical(identical(a, b), TRUE)
}
)

test_that("xegaGpFilterCrossGene OK", 
{
 set.seed(17)
gene1<-xegaGpInitGene(lFxegaGpGene)
gene2<-xegaGpInitGene(lFxegaGpGene)
a<-xegaGpDecodeGene(gene1, lFxegaGpGene)
gene<-xegaGpFilterCrossGene(gene1, gene2, lFxegaGpGene)
b<-xegaGpDecodeGene(gene[[1]], lFxegaGpGene)
 expect_identical(identical(a, b), FALSE)
}
)

test_that("xegaGpCrossoverFactory Cross2Gene OK",
 {
 f<-xegaGpCrossoverFactory()
 expect_identical(body(f), body(xegaGpAllCross2Gene))
}
)

test_that("xegaGpCrossoverFactory Cross2Gene OK",
 {
 f<-xegaGpCrossoverFactory(method="Cross2Gene")
 expect_identical(body(f), body(xegaGpAllCross2Gene))
}
)

test_that("xegaGpCrossoverFactory AllCross2Gene OK",
 {
 f<-xegaGpCrossoverFactory(method="AllCross2Gene")
 expect_identical(body(f), body(xegaGpAllCross2Gene))
}
)

test_that("xegaGpCrossoverFactory AllCrossGene OK",
 {
 f<-xegaGpCrossoverFactory(method="AllCrossGene")
 expect_identical(body(f), body(xegaGpAllCrossGene))
}
)

test_that("xegaGpCrossoverFactory FilterCross2Gene OK",
 {
 f<-xegaGpCrossoverFactory(method="FilterCross2Gene")
 expect_identical(body(f), body(xegaGpFilterCross2Gene))
}
)

test_that("xegaGpCrossoverFactory FilterCrossGene OK",
 {
 f<-xegaGpCrossoverFactory(method="FilterCrossGene")
 expect_identical(body(f), body(xegaGpFilterCrossGene))
}
)

test_that("xegaGpCrossoverFactory sgunknown OK",
 {
 expect_error(
 xegaGpCrossoverFactory(method="sgunknown"),
 "sgunknown")
}
)

