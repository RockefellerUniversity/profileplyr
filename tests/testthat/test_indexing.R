library(profileplyr)
library(soGGi)
context("Test subsetting and indexing")
data(chipExampleBig)
p <- as_profileplyr(chipExampleBig,names = "name")
pRowFilt <- p[1:200,]
pColumnFilt <- p[,1:100]
pSampleFilt <- p[,,3:4]

expect_true(is(pRowFilt,"profileplyr"))
expect_true(is(pColumnFilt,"profileplyr"))
expect_true(is(pSampleFilt,"profileplyr"))