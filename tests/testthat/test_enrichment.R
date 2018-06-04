context("Enrichment test outputs")
library(genesetr)

testthat::test_that("fastFET() is good approximation of R implementation test.fisher()", {
  n=1000
  #generate 1000 random contingency tables
  a = round(sample(0:1000, n))
  b = round(sample(0:5000, n))
  c = round(sample(0:5000, n))
  d = round(sample(20000:30000, n))

  test = fastFET(a,b,c,d, alternative = "greater")
  true = numeric(length(a))
  for(i in 1:length(a)) true[i] =
    fisher.test(matrix(c(a[i],c[i],b[i],d[i]),ncol = 2),
      alternative = 'greater')$p.value
  diff = abs(true-test)
  expect_lt(max(diff),1e-8)
})

testthat::test_that("OR computation is accurate",{

})

testthat::test_that("Jaccard distance computation is accurate",{

})

testthat::test_that("Jaccard similarity computation is accurate",{

})

testthat::test_that("Euclidean distance computation is accurate",{

})




