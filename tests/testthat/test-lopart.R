library(testthat)
library(gfpop)
context("lopart")

## LOPART as a time-dependent gfpop graph: 2 states (normal, noChange) and
## 4 rules picking the active edges for each label context.
lopart_graph <- function(lambda) {
  graph(
    ## rule 1: unlabeled
    Edge("normal",   "normal",   "null", rule = 1),
    Edge("normal",   "normal",   "std",  rule = 1, penalty = lambda),
    ## rule 2: in positive label, at most one change
    Edge("normal",   "normal",   "null", rule = 2),
    Edge("noChange", "noChange", "null", rule = 2),
    Edge("normal",   "noChange", "std",  rule = 2, penalty = lambda),
    ## rule 3: end of positive label
    Edge("normal",   "normal",   "std",  rule = 3, penalty = lambda),
    Edge("noChange", "normal",   "null", rule = 3),
    ## rule 4: in negative label, no change
    Edge("normal",   "normal",   "null", rule = 4),
    StartEnd(start = "normal", end = c("normal", "noChange")))
}

## labels data.frame (start, end, changes) -> per-point rule vector
lopart_rule_vec <- function(n, labels) {
  rule <- rep(1L, n)                          # unlabeled
  for (i in seq_len(nrow(labels))) {
    s <- labels$start[i]
    e <- labels$end[i]
    if (labels$changes[i] == 1) {             # positive: interior = 2, end = 3
      rule[s:(e - 1)] <- 2L
      rule[e] <- 3L
    } else {                                  # negative: no change
      rule[s:e] <- 4L
    }
  }
  rule
}

## positive label over a real change, negative label over a flat region
make_scenario <- function() {
  set.seed(1)
  x <- c(rnorm(25, 0), rnorm(25, 5), rnorm(25, 5), rnorm(25, 0))
  labels <- data.frame(start = c(20, 55), end = c(30, 70), changes = c(1, 0))
  list(x = x, labels = labels, lambda = 2 * log(100))
}

test_that("gfpop with a rule vector matches LOPART::LOPART()", {
  skip_if_not_installed("LOPART")
  s <- make_scenario()
  lo <- LOPART::LOPART(s$x, s$labels,
                       penalty_unlabeled = s$lambda, penalty_labeled = s$lambda)
  gf <- gfpop(s$x, lopart_graph(s$lambda), type = "mean",
              rule = lopart_rule_vec(length(s$x), s$labels))

  ## same segment ends and means
  expect_equal(as.integer(gf$changepoints), as.integer(lo$segments$end))
  expect_equal(as.numeric(gf$parameters), as.numeric(lo$segments$mean),
               tolerance = 1e-6)
  ## gfpop globalCost = sum((x-mu)^2); LOPART omits the constant sum(x^2)
  expect_equal(gf$globalCost, sum(s$x^2) + lo$loss$total_loss, tolerance = 1e-6)
})

test_that("rule filtering enforces the label structural invariants", {
  s <- make_scenario()
  gf <- gfpop(s$x, lopart_graph(s$lambda), type = "mean",
              rule = lopart_rule_vec(length(s$x), s$labels))
  internal <- gf$changepoints[-length(gf$changepoints)]
  ## positive label: one change
  expect_equal(sum(internal >= 20 & internal <= 30), 1)
  ## negative label: no change
  expect_equal(sum(internal >= 55 & internal <= 70), 0)
})

test_that("an all-unlabeled rule vector reproduces graph(type = std)", {
  s <- make_scenario()
  g_std <- gfpop(s$x, graph(type = "std", penalty = s$lambda), type = "mean")
  g_lab <- gfpop(s$x, lopart_graph(s$lambda), type = "mean",
                 rule = rep(1L, length(s$x)))
  expect_equal(as.integer(g_lab$changepoints), as.integer(g_std$changepoints))
  expect_equal(as.numeric(g_lab$parameters), as.numeric(g_std$parameters),
               tolerance = 1e-6)
})
