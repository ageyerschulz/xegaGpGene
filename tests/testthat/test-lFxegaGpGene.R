
library(testthat)
library(xegaSelectGene)
library(xegaGpGene)

test_that("lFxegaGpGene OK", 
{
           expect_identical(lFxegaGpGene$penv$name(), "envXOR")
           expect_equal(lFxegaGpGene$replay(), 0)
           expect_equal(lFxegaGpGene$verbose(), 4)
           expect_equal(lFxegaGpGene$MaxDepth(), 5)
           expect_equal(lFxegaGpGene$MaxTrials(), 5)
           expect_equal(lFxegaGpGene$MutationRate(), 0.05)
           expect_equal(lFxegaGpGene$CrossRate(), 0.2)
           expect_equal(lFxegaGpGene$Max(), 1)
           expect_equal(lFxegaGpGene$Offset(), 1)
           expect_equal(lFxegaGpGene$Eps(), 0.01)
           expect_identical(lFxegaGpGene$Elitist(), TRUE)
           expect_equal(lFxegaGpGene$TournamentSize(), 2)
}
)

