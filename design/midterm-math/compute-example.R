## compute-example.R
## Generates every number used in midterm-math.tex for the 10-step worked
## example. Writes LaTeX fragments into design/midterm-math/tables/.
##
## Sources:
##   - a self-contained functional-pruning DP (reference implementation of the
##     recursion described in the document) -> the Q-table and the optimal path
##   - gfpop(..., rule = sigma)             -> final segmentation (cross-check)
##   - LOPART::LOPART()                     -> oracle (cross-check)
## The three are asserted mutually consistent before anything is written.

suppressMessages(library(gfpop))
stopifnot("rule" %in% names(formals(gfpop::gfpop)))
here <- "design/midterm-math/tables"
dir.create(here, showWarnings = FALSE, recursive = TRUE)

## ---- the locked dataset ---------------------------------------------------
x      <- c(0.2, -0.1, 0.1, 0.0, 5.1, 4.9, 5.0, 5.2, 4.8, 5.0)
n      <- length(x)
lambda <- 2 * log(n)
## labels: one positive region [3,6] (exactly one change), one negative
## region [8,10] (no change).
labels <- data.frame(start = c(3, 8), end = c(6, 10), changes = c(1, 0))

## label context + rule vector (the sigma map)
context <- rep("unlabeled", n)
sigma   <- rep(1L, n)
for (i in seq_len(nrow(labels))) {
  s <- labels$start[i]; e <- labels$end[i]
  if (labels$changes[i] == 1) {
    sigma[s:(e - 1)] <- 2L; sigma[e] <- 3L
    context[s:(e - 1)] <- "in-pos"; context[e] <- "end-pos"
  } else {
    sigma[s:e] <- 4L; context[s:e] <- "in-neg"
  }
}

## ---- the graph (states normal, noChange; 4 rules) -------------------------
lopart_graph <- function(lam, ends = c("normal", "noChange")) gfpop::graph(
  Edge("normal",   "normal",   "null", rule = 1),
  Edge("normal",   "normal",   "std",  rule = 1, penalty = lam),
  Edge("normal",   "normal",   "null", rule = 2),
  Edge("noChange", "noChange", "null", rule = 2),
  Edge("normal",   "noChange", "std",  rule = 2, penalty = lam),
  Edge("normal",   "normal",   "std",  rule = 3, penalty = lam),
  Edge("noChange", "normal",   "null", rule = 3),
  Edge("normal",   "normal",   "null", rule = 4),
  StartEnd(start = "normal", end = ends))

## edge table used by the reference DP (from, to, type, rule, penalty)
edges <- data.frame(
  from    = c("normal","normal","normal","noChange","normal","normal","noChange","normal"),
  to      = c("normal","normal","normal","noChange","noChange","normal","normal","normal"),
  type    = c("null","std","null","null","std","std","null","null"),
  rule    = c(1L,1L,2L,2L,2L,3L,3L,4L),
  penalty = c(0,lambda,0,0,lambda,lambda,0,0),
  stringsAsFactors = FALSE)
states <- c("normal", "noChange")

## ---- functional-pruning reference DP --------------------------------------
## A cost function is a set of quadratics q(mu) = A mu^2 + B mu + C, its value
## the lower envelope min_k q_k(mu). Each quadratic carries a backpointer.
minval <- function(A, B, C) if (length(A) == 0) Inf else min(C - B^2 / (4 * A))

## store[[k]] : data.frame of quadratics live at time k
store <- vector("list", n)
## k = 1: only the start state, Q_1^normal(mu) = (x_1 - mu)^2
store[[1]] <- data.frame(A = 1, B = -2 * x[1], C = x[1]^2,
                         state = "normal", parent = NA_integer_, ctype = "init",
                         stringsAsFactors = FALSE)

Qmin <- matrix(Inf, nrow = n, ncol = 2, dimnames = list(NULL, states))
Qmin[1, "normal"] <- minval(store[[1]]$A, store[[1]]$B, store[[1]]$C)

for (k in 2:n) {
  prev <- store[[k - 1]]
  rows <- list()
  active <- edges[edges$rule == sigma[k], ]
  for (j in seq_len(nrow(active))) {
    e <- active[j, ]
    src <- which(prev$state == e$from)
    if (length(src) == 0) next                    # source unreachable -> skip
    if (e$type == "null") {                        # continue segment (copy + point)
      for (r in src) {
        rows[[length(rows) + 1]] <- data.frame(
          A = prev$A[r] + 1, B = prev$B[r] - 2 * x[k], C = prev$C[r] + x[k]^2,
          state = e$to, parent = r, ctype = "null", stringsAsFactors = FALSE)
      }
    } else {                                       # std: collapse source envelope
      vals <- prev$C[src] - prev$B[src]^2 / (4 * prev$A[src])
      best <- src[which.min(vals)]
      m <- min(vals) + e$penalty
      rows[[length(rows) + 1]] <- data.frame(       # new segment: {0,0,m}+point
        A = 1, B = -2 * x[k], C = m + x[k]^2,
        state = e$to, parent = best, ctype = "std", stringsAsFactors = FALSE)
    }
  }
  store[[k]] <- if (length(rows)) do.call(rbind, rows) else
    prev[0, , drop = FALSE]
  for (s in states) {
    idx <- store[[k]]$state == s
    Qmin[k, s] <- minval(store[[k]]$A[idx], store[[k]]$B[idx], store[[k]]$C[idx])
  }
}

## final answer + backtrack the optimal path
fin  <- store[[n]]
fval <- fin$C - fin$B^2 / (4 * fin$A)
star <- which.min(fval)
refCost <- min(fval)

path <- character(n); changes <- integer(0)
row <- star
for (k in n:1) {
  path[k] <- store[[k]]$state[row]
  if (store[[k]]$ctype[row] == "std") changes <- c(changes, k - 1)  # seg starts at k
  row <- store[[k]]$parent[row]
  if (is.na(row)) break
}
refChangepoints <- as.integer(sort(unique(c(changes, n))))

## ---- gfpop + LOPART (cross-checks) ----------------------------------------
fit <- gfpop::gfpop(x, lopart_graph(lambda), type = "mean", rule = sigma)
gfCp   <- as.integer(fit$changepoints)
gfMean <- as.numeric(unlist(fit$parameters))

lo <- LOPART::LOPART(x, labels, penalty_unlabeled = lambda, penalty_labeled = lambda)
loCp   <- as.integer(lo$segments$end)
loMean <- as.numeric(lo$segments$mean)

## gfpop reports globalCost as residual SSE (no penalties); the DP objective
## adds K*lambda for the K changes. Reconcile before comparing.
nChanges <- length(gfCp) - 1
penObj   <- fit$globalCost + nChanges * lambda

## ---- consistency assertions (fail loudly) ---------------------------------
stopifnot(
  isTRUE(all.equal(refCost, penObj)),                   # reference DP == gfpop objective
  identical(refChangepoints, gfCp),                     # same changepoints
  identical(gfCp, loCp),                                # gfpop == LOPART (oracle)
  isTRUE(all.equal(gfMean, loMean, tolerance = 1e-6)),
  isTRUE(all.equal(fit$globalCost, sum(x^2) + lo$loss$total_loss)) # loss identity
)
cat("All cross-checks passed. globalCost =", fit$globalCost, "\n")
cat("optimal path:", paste(path, collapse = ","), "\n")

## ---- emit LaTeX fragments -------------------------------------------------
fnum <- function(v, d = 2) formatC(v, format = "f", digits = d)
finf <- function(v, d = 2) ifelse(is.infinite(v), "$+\\infty$", fnum(v, d))
w <- function(file, lines) writeLines(lines, file.path(here, file))
## table bodies are wrapped in a macro: \input inside a tabular misbehaves,
## but an expanded \newcommand does not.
wmac <- function(file, macro, rows)
  w(file, c(paste0("\\newcommand{\\", macro, "}{"), rows, "}"))

## scalar macros
segRanges <- {
  ends <- gfCp; starts <- c(1, head(ends, -1) + 1)
  paste0("[", starts, ",", ends, "]")
}
ctxAbbr <- c(unlabeled = "unlabeled", "in-pos" = "in pos.", "end-pos" = "end pos.",
             "in-neg" = "in neg.")
w("scalars.tex", c(
  paste0("\\newcommand{\\dataN}{", n, "}"),
  paste0("\\newcommand{\\dataLambda}{", fnum(lambda, 3), "}"),
  paste0("\\newcommand{\\dataX}{", paste(fnum(x, 1), collapse = ",\\ "), "}"),
  paste0("\\newcommand{\\dataSigma}{", paste(sigma, collapse = ",\\ "), "}"),
  paste0("\\newcommand{\\sse}{", fnum(fit$globalCost, 3), "}"),
  paste0("\\newcommand{\\penObj}{", fnum(penObj, 3), "}"),
  paste0("\\newcommand{\\lopartLoss}{", fnum(lo$loss$total_loss, 3), "}"),
  paste0("\\newcommand{\\lopartLossAbs}{", fnum(abs(lo$loss$total_loss), 3), "}"),
  paste0("\\newcommand{\\sumsq}{", fnum(sum(x^2), 3), "}"),
  paste0("\\newcommand{\\numChanges}{", nChanges, "}"),
  paste0("\\newcommand{\\pathStates}{", paste(path, collapse = ","), "}")
))

## sigma / active-edge table
edgeShort <- c("1" = "$n\\!\\to\\!n$ null,\\ $n\\!\\to\\!n$ std",
               "2" = "$n\\!\\to\\!n$ null,\\ $c\\!\\to\\!c$ null,\\ $n\\!\\to\\!c$ std",
               "3" = "$n\\!\\to\\!n$ std,\\ $c\\!\\to\\!n$ null",
               "4" = "$n\\!\\to\\!n$ null")
wmac("sigma-table.tex", "sigmaRows", c(
  vapply(1:n, function(t) sprintf(
    "%d & %s & %s & %d & %s \\\\", t, fnum(x[t], 1),
    ctxAbbr[[context[t]]], sigma[t], edgeShort[[as.character(sigma[t])]]),
    character(1))
))

## Q-table (per-state functional minima; +infinity where unreachable)
wmac("q-table.tex", "qRows", c(
  vapply(1:n, function(t) sprintf(
    "%d & %s & %s & %s \\\\", t, fnum(x[t], 1),
    finf(Qmin[t, "normal"]), finf(Qmin[t, "noChange"])),
    character(1))
))

## final segmentation
wmac("segmentation.tex", "segRows", c(
  vapply(seq_along(gfCp), function(i) sprintf(
    "%d & %s & %s \\\\", i, segRanges[i], fnum(gfMean[i], 3)),
    character(1))
))

## oracle comparison
ll  <- lo$loss$total_loss
mid <- if (ll < 0) sprintf("%s - %s", fnum(sum(x^2), 3), fnum(abs(ll), 3)) else sprintf("%s + %s", fnum(sum(x^2), 3), fnum(ll, 3))
wmac("oracle.tex", "oracleRows", c(
  sprintf("changepoints & $\\{%s\\}$ & $\\{%s\\}$ \\\\",
          paste(gfCp, collapse = ","), paste(loCp, collapse = ",")),
  sprintf("segment means & $%s$ & $%s$ \\\\",
          paste(fnum(gfMean, 3), collapse = ",\\,"),
          paste(fnum(loMean, 3), collapse = ",\\,")),
  sprintf("residual SSE & $%s$ & $%s = %s$ \\\\",
          fnum(fit$globalCost, 3), mid, fnum(sum(x^2) + ll, 3))
))

## plot coordinates as macros (expand cleanly inside TikZ "plot coordinates")
w("points.tex", paste0("\\newcommand{\\pointsCoords}{",
   paste(sprintf("(%d,%s)", 1:n, fnum(x, 3)), collapse = ""), "}"))
## fitted mean as a staircase (horizontal per segment, vertical riser at joins)
pts <- sprintf("(%.2f,%s)", 0.5, fnum(gfMean[1], 3))
for (i in seq_along(gfCp)) {
  b <- gfCp[i] + 0.5
  pts <- c(pts, sprintf("(%.2f,%s)", b, fnum(gfMean[i], 3)))
  if (i < length(gfCp)) pts <- c(pts, sprintf("(%.2f,%s)", b, fnum(gfMean[i + 1], 3)))
}
w("fit-step.tex", paste0("\\newcommand{\\fitCoords}{", paste(pts, collapse = ""), "}"))

cat("Wrote fragments to", here, "\n")
