
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("EnvXOR OK", 
{
           envXOR<-newEnvXOR() 
           a2<-"OR(OR(D1, D2), (AND(NOT(D1), NOT(D2))))"
           a3<-"OR(OR(D1, D2), AND(D1, D2))"
           a4<-"AND(OR(D1,D2),NOT(AND(D1,D2)))"
           gp4<-"(AND(AND(OR(D2,D1),NOT(AND(D1,D2))),(OR(D2,D1))))"
	   expect_identical(envXOR$name(), "envXOR")
           expect_identical(envXOR$f(a2), 2)
           expect_identical(envXOR$f(a3), 3)
           expect_identical(envXOR$f(a4), 4)
           expect_identical(envXOR$f(gp4), 4)
           expect_gt(envXOR$f(gp4, gene=xegaGpInitGene(lFxegaGpGene)), 3)
}
)

