#test reading in Str file
#unit test 
test_that("readStr file", {
  result <- readStr(strfile)
  expected <- strfile
  expect_equal(result,expected)
})

#contigency test 1
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(strfile)
  expect_s3_class(result, expected_data_type)
})

#test reading in indiv file
#unit test
test_that("readIndivq", {
  expected <- indivqfile[,6:ncol(indivqfile)]
  result <- readIndivq(indivqfile)
  expect_equal(result, expected)
})

#contigency test 1: the sum of the values after the fifth column adds up to 1
test_sum_of_indivqfile <- function(indivqfile){
  result <- sum(indivqfile[, 6:ncol(indivqfile)]) == 1
  expected <- TRUE
  if (!identical(result, expected)){
    stop("Test Failed")
  }
}

#contigency test 2
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(indivqfile)
  expect_s3_class(result, expected_data_type)
})

#test reading in pkla file
test_that("readPkla", {
  expected <- pklafile[,3:ncol(pklafile)]
  result <- readPkla(pklafile)
  expect_equal(result, expected)
})

#contigency test 1
test_that("check that datatset is a dataframe", {
  expected_data_type <- "data.frame"
  result <- readStr(pklafile)
  expect_s3_class(result, expected_data_type)
})

#contigency test 2/ check that across columns after the second for one locus adds up to 1
sum_of_unique_values <- function(pklafile, k) {
  unique_value <- unique(pklafile[, "V1"])
  result <- data.frame(V1 = unique_value)
  for (j in 2:k) {
    colmn_sum <- numeric()
    for (i in 1:length(unique_value)) {
      colmn_sum[i] <- sum(pklafile[pklafile[, "V1"] == unique_value[i], j])
    }
    result[, j] <- colmn_sum
  }
  return(result)
}

test_that("sum of unique values in V1 is equal to 1", {
  expected <- sum_of_unique_values(pklafile)
  result <-
})

#test for getting unique individuals from the strfile
#unit test
test_that("get unique indiv", {
  result <- getIndiv(strfile)
  expected <- as.numeric(unique(strfile$V1))
  expect_equal(result,expected)
})

#test for getting K subpopulation from indinvq file
#unit test
test_that("get k subpop", {
  result <- getSubpop(indivqfile)
  expected <- indivqfile[,6:ncol(indivqfile)]
  expect_equal(result,expected)
})

#check that the file is a dataframe
test_that("check dataframe", {
  result <- getLoci(indivqfile)
  expect_true(is.data.frame(indivqfile))
})

#test get unique alleles from pklafile
#unit test
test_that("get unique alleles", {
  result <- getAlleles(pklafile)
  expected <- aggregate(x = pklafile$V2, by = list(pklafile$V1), FUN = function(x) length(unique(x)))
  expect_equal(result,expected)
})

#check that the file is a dataframe
test_that("check dataframe", {
  result <- getAlleles(pklafile)
  expect_true(is.data.frame(pklafile))
})

#test for getting loci from strfile
#unit test
test_that("get loci loci", {
  result <- getLoci(strfile)
  expected <- strfile[,3:ncol(strfile)]
  expect_equal(result,expected)
})
 
#check that the file is a dataframe
test_that("check dataframe", {
  result <- getLoci(strfile)
  expect_true(is.data.frame(strfile))
})

#test functions for Fst function from Pegas
# test is using Jaquar dataset from R package
#test to check for diploid input
test_that("check for non-diploid input", {
  x <- jaguar
  expect_error(Fst(x), "Fst() requires diploid data input")
})

#check for population assignment
test_that("check for population assignment", {
  x <- jaguar
  expect_equal(Fst(x, pop = as.factor(4)), Fst(x, pop = as.factor(4)))
})

#test function for computenumpairs
#indiv function use in computenumpairs
getIndiv <- function(strfile){
  individual <- as.numeric(unique(strfile$V1));
  return(individual)
}

#unit test
test_that("compute number of pairs function properly", {
  result <- computenumpairs(as.numeric(indivqfile))
  expected <- getIndiv(as.numeric(indivqfile*(indivqfile-1)/2))
  expect_equal(result,expected)
})
