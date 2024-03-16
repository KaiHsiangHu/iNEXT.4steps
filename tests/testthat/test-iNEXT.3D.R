context("iNEXT.4steps")
test_that("iNEXT4steps for abundance-based data", {
  
  # Test input by a demo data
  data("Data_spider")
  out <- iNEXT4steps(Data_spider$Closed, datatype = "abundance", q = c(0,1))
  expect_equal(length(out), 2)
  expect_equal(length(out$figure), 6)
  
})


test_that("iNEXT4steps for sampling-unit-based incidence raw data", {
  
  # Test input by a demo data
  data("Data_woody_plant")
  out <- iNEXT4steps(Data_woody_plant$Lowland, q = c(1,2), datatype = "incidence_raw")
  expect_equal(length(out), 2)
  expect_equal(length(out$figure), 6)
  
})


test_that("Completeness in a single assemblage", {
  
  # Test input by a demo data
  data("Data_spider")
  out <- Completeness(Data_spider$Open, datatype = "abundance", q = 1)
  expect_equal(nrow(out), 1)
  
  # Test input by a demo data
  data("Data_woody_plant")
  out <- Completeness(Data_woody_plant$Monsoon, q = 2, datatype = "incidence_raw")
  expect_equal(nrow(out), 1)
  
})


test_that("Evenness in a single assemblage", {
  
  # Test input by a demo data
  data("Data_spider")
  out <- Evenness(Data_spider$Open, datatype = "abundance", q = 1, E.class = 2, method = "Observed")
  expect_equal(length(out), 1)
  expect_equal(nrow(out[[1]]), 1)
  
  data("Data_spider")
  out <- Evenness(Data_spider$Closed, datatype = "abundance", q = 2, E.class = 2, method = "Estimated", SC = 0.99)
  expect_equal(length(out), 1)
  expect_equal(nrow(out[[1]]), 1)
  
  
  # Test input by a demo data
  data("Data_woody_plant")
  out <- Evenness(Data_woody_plant$Upper_cloud, q = 0, datatype = "incidence_raw", E.class = 2, method = "Observed")
  expect_equal(length(out), 1)
  expect_equal(nrow(out[[1]]), 1)
  
  data("Data_woody_plant")
  out <- Evenness(Data_woody_plant$Monsoon, q = 2, datatype = "incidence_raw", E.class = 2, method = "Estimated", SC = 0.99)
  expect_equal(length(out), 1)
  expect_equal(nrow(out[[1]]), 1)
  
})

