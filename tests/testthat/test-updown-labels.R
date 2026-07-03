library(testthat)
library(gfpop)
context("updown-labels")

## Up-down peak detection with labels: 4 states (Dw, Up, noChangeUp, noChangeDw)
## and 6 rules. peakStart holds one up, peakEnd one down; noPeaks is flat.
## No StartEnd, so the rule-1 subset is exactly graph(type = "updown").
updown_labels_graph <- function(lambda) {
  graph(
    ## rule 1: unlabeled = standard up-down
    Edge("Dw", "Dw", "null", rule = 1),
    Edge("Up", "Up", "null", rule = 1),
    Edge("Dw", "Up", "up",   rule = 1, penalty = lambda),
    Edge("Up", "Dw", "down", rule = 1, penalty = lambda),
    ## rule 2: noPeaks, hold the level
    Edge("Dw", "Dw", "null", rule = 2),
    Edge("Up", "Up", "null", rule = 2),
    ## rule 3: in peakStart, at most one up
    Edge("Dw", "Dw", "null", rule = 3),
    Edge("Dw", "noChangeUp", "up", rule = 3, penalty = lambda),
    Edge("noChangeUp", "noChangeUp", "null", rule = 3),
    ## rule 4: end of peakStart, force the up into Up
    Edge("Dw", "Up", "up", rule = 4, penalty = lambda),
    Edge("noChangeUp", "Up", "null", rule = 4),
    ## rule 5: in peakEnd, at most one down
    Edge("Up", "Up", "null", rule = 5),
    Edge("Up", "noChangeDw", "down", rule = 5, penalty = lambda),
    Edge("noChangeDw", "noChangeDw", "null", rule = 5),
    ## rule 6: end of peakEnd, force the down into Dw
    Edge("Up", "Dw", "down", rule = 6, penalty = lambda),
    Edge("noChangeDw", "Dw", "null", rule = 6))
}

## labels data.frame (start, end, type in {peakStart, peakEnd, noPeaks})
## -> per-point rule vector
updown_labels_rule_vec <- function(n, labels) {
  rule <- rep(1L, n)                          # unlabeled
  for (i in seq_len(nrow(labels))) {
    s <- labels$start[i]
    e <- labels$end[i]
    if (labels$type[i] == "peakStart") {      # interior = 3, end = 4
      rule[s:(e - 1)] <- 3L
      rule[e] <- 4L
    } else if (labels$type[i] == "peakEnd") { # interior = 5, end = 6
      rule[s:(e - 1)] <- 5L
      rule[e] <- 6L
    } else {                                  # noPeaks: no change
      rule[s:e] <- 2L
    }
  }
  rule
}

## synthetic peak: baseline, rise, plateau, fall back to baseline.
## peakStart brackets the rise, peakEnd the fall, noPeaks the flat plateau.
make_peak_scenario <- function() {
  set.seed(1)
  x <- c(rnorm(25, 0), rnorm(25, 5), rnorm(25, 5), rnorm(25, 0))
  labels <- data.frame(
    start = c(20, 40, 70),
    end   = c(30, 60, 80),
    type  = c("peakStart", "noPeaks", "peakEnd"),
    stringsAsFactors = FALSE)
  list(x = x, labels = labels, lambda = 2 * log(100))
}

test_that("an all-unlabeled rule vector reproduces graph(type = updown)", {
  s <- make_peak_scenario()
  g_ud  <- gfpop(s$x, graph(type = "updown", penalty = s$lambda), type = "mean")
  g_lab <- gfpop(s$x, updown_labels_graph(s$lambda), type = "mean",
                 rule = rep(1L, length(s$x)))
  expect_equal(as.integer(g_lab$changepoints), as.integer(g_ud$changepoints))
  expect_equal(as.numeric(g_lab$parameters), as.numeric(g_ud$parameters),
               tolerance = 1e-6)
  expect_equal(g_lab$globalCost, g_ud$globalCost, tolerance = 1e-6)
})

test_that("rule filtering enforces the peak label invariants", {
  s <- make_peak_scenario()
  gf <- gfpop(s$x, updown_labels_graph(s$lambda), type = "mean",
              rule = updown_labels_rule_vec(length(s$x), s$labels))
  internal <- gf$changepoints[-length(gf$changepoints)]
  ## peakStart [20, 30]: exactly one change (the rise)
  expect_equal(sum(internal >= 20 & internal <= 30), 1)
  ## peakEnd [70, 80]: exactly one change (the fall)
  expect_equal(sum(internal >= 70 & internal <= 80), 1)
  ## noPeaks [40, 60]: no change
  expect_equal(sum(internal >= 40 & internal <= 60), 0)
  ## and the fit is a peak: the highest segment mean is interior, not an end
  expect_gt(which.max(gf$parameters), 1)
  expect_lt(which.max(gf$parameters), length(gf$parameters))
})
