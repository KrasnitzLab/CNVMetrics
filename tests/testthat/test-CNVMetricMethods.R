### Unit tests for CNVMetricsMethods.R functions

library(CNVMetrics)

### Tests prepareInformation() results

context("prepareInformation() results")

test_that("prepareInformation() must return error when segDirectory is not a string", {
    
    error_message <- "chrInfo must be a Seqinfo object."
    
    expect_error(prepareInformation(segDirectory = "test", chrInfo = 33), 
                 error_message)
    
})