library(testthat)
library(gfpop)
context("rule")

test_that("Edge has an integer rule column, NA by default", {
  e <- Edge("a", "b", "std")
  expect_true("rule" %in% names(e))
  expect_identical(names(e)[length(names(e))], "rule")
  expect_type(e$rule, "integer")
  expect_true(is.na(e$rule))
})

test_that("Edge(rule=) sets the rule id", {
  e <- Edge("a", "b", "std", rule = 2)
  expect_identical(e$rule, 2L)
})

test_that("Edge is vectorized over rule", {
  e <- Edge(c("a", "c"), c("b", "d"), "std", rule = 1)
  expect_identical(e$rule, c(1L, 1L))
})

test_that("Edge rejects non-positive rule ids", {
  expect_error(Edge("a", "b", "std", rule = 0), "rule must be a positive integer")
  expect_error(Edge("a", "b", "std", rule = -1), "rule must be a positive integer")
})

test_that("StartEnd and Node carry a rule column matching Edge", {
  se <- StartEnd("a", "a")
  nd <- Node("a", 0, 1)
  e <- Edge("a", "a", "std")
  expect_identical(names(se), names(e))
  expect_identical(names(nd), names(e))
  expect_true(is.na(se$rule[1]))
  expect_true(is.na(nd$rule[1]))
})

test_that("predefined graphs have an all-NA rule column", {
  g <- graph(type = "updown", penalty = 5)
  expect_true("rule" %in% names(g))
  expect_true(all(is.na(g$rule)))
})

test_that("custom graph preserves per-edge rule values", {
  g <- graph(
    Edge("a", "a", "null", rule = 1),
    Edge("a", "a", "std", penalty = 5, rule = 2))
  expect_identical(g$rule, c(1L, 2L))
})

test_that("graphReorder preserves rule values with their edges", {
  g <- graph(
    Edge("a", "a", "null", rule = 1),
    Edge("a", "a", "std", penalty = 5, rule = 2),
    StartEnd("a", "a"))
  reordered <- gfpop:::graphReorder(g)$graph
  expect_true("rule" %in% names(reordered))
  expect_identical(reordered$rule[reordered$type == "null"], 1L)
  expect_identical(reordered$rule[reordered$type == "std"], 2L)
})

test_that("gfpop still solves a standard problem with the rule column present", {
  set.seed(42)
  x <- c(rnorm(40, 0), rnorm(40, 10))
  g <- graph(type = "std", penalty = 2 * log(80))
  fit <- gfpop(x, g, type = "mean")
  expect_s3_class(fit, "gfpop")
  ## last changepoint is the data length
  expect_identical(fit$changepoints[length(fit$changepoints)], 80L)
  ## one internal change near the true location (index 40)
  internal <- fit$changepoints[-length(fit$changepoints)]
  expect_true(any(internal >= 35 & internal <= 45))
  ## parameters and changepoints are aligned in length
  expect_identical(length(fit$parameters), length(fit$changepoints))
})

test_that("rule column is read by the solver but is inert (no effect yet)", {
  set.seed(7)
  x <- c(rnorm(30, 0), rnorm(30, 8))
  p <- 2 * log(60)
  ## same graph, only the rule column differs
  g_na <- graph(
    Edge("a", "a", "null"),
    Edge("a", "a", "std", penalty = p))
  g_r1 <- graph(
    Edge("a", "a", "null", rule = 1),
    Edge("a", "a", "std", penalty = p, rule = 1))
  fit_na <- gfpop(x, g_na, type = "mean")
  fit_r1 <- gfpop(x, g_r1, type = "mean")
  expect_identical(fit_na$changepoints, fit_r1$changepoints)
  expect_equal(fit_na$parameters, fit_r1$parameters)
  expect_equal(fit_na$globalCost, fit_r1$globalCost)
})

## rule 1 = change allowed, rule 2 = stay only
rule_graph <- function(p) {
  graph(
    Edge("a", "a", "null"),                       # NA: always active
    Edge("a", "a", "std", penalty = p, rule = 1), # rule 1: change allowed
    Edge("a", "a", "null", rule = 2))             # rule 2: stay only
}

test_that("no rule argument equals all-active (rule = 1 everywhere)", {
  set.seed(11)
  x <- c(rnorm(40, 0), rnorm(40, 10))
  g <- rule_graph(2 * log(length(x)))
  fit_norule <- gfpop(x, g, type = "mean")
  fit_all1   <- gfpop(x, g, type = "mean", rule = rep(1L, length(x)))
  expect_identical(fit_norule$changepoints, fit_all1$changepoints)
  expect_equal(fit_norule$parameters, fit_all1$parameters)
  expect_equal(fit_norule$globalCost, fit_all1$globalCost)
})

test_that("rule filtering changes the segmentation", {
  set.seed(11)
  x <- c(rnorm(40, 0), rnorm(40, 10))
  g <- rule_graph(2 * log(length(x)))
  ## rule 1: std edge can fire
  fit_change <- gfpop(x, g, type = "mean", rule = rep(1L, length(x)))
  ## rule 2: stay only -> one segment
  fit_block  <- gfpop(x, g, type = "mean", rule = rep(2L, length(x)))
  expect_true(length(fit_change$changepoints) >= 2)
  expect_identical(length(fit_block$changepoints), 1L)
  expect_identical(fit_block$changepoints, length(x))
})

test_that("gfpop validates the rule argument", {
  x <- c(rnorm(20, 0), rnorm(20, 5))
  g <- rule_graph(2 * log(length(x)))
  expect_error(gfpop(x, g, type = "mean", rule = rep(1L, length(x) - 1)),
               "different sizes")
  bad.na <- rep(1L, length(x)); bad.na[3] <- NA
  expect_error(gfpop(x, g, type = "mean", rule = bad.na), "missing values")
  expect_error(gfpop(x, g, type = "mean", rule = rep(0L, length(x))),
               "positive integers")
  expect_error(gfpop(x, g, type = "mean", rule = rep(9L, length(x))),
               "not present in the graph")
})
