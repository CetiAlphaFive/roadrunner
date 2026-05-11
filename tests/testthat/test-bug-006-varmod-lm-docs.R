# BUG-006 regression test (v0.0.0.9029).
#
# The original docstring at R/ares.R promised "simple mean-dependent
# heteroscedasticity" for varmod = "lm", which the audit showed was an
# over-promise: only YHAT-dependent residual scale is captured, and
# x-driven heteroscedasticity (variance that varies with a covariate
# orthogonal to the mean structure) is missed entirely. The fix is doc-
# only; this test guards the docstring so future copy-edits don't
# silently regress the language back to over-promising.

test_that("varmod='lm' docstring acknowledges the yhat-only limitation", {
  rd_path <- system.file("help", "AnIndex", package = "roadrunner")
  # Read the source Rd file (the rendered help text isn't easy to assert
  # against in-package; the rd source is the authoritative copy).
  rd_src <- file.path(testthat::test_path("..", ".."), "man", "ares.Rd")
  if (!file.exists(rd_src)) skip("man/ares.Rd not found")
  rd_text <- paste(readLines(rd_src, warn = FALSE), collapse = "\n")
  # Two assertions:
  #   1. "yhat-dependent" appears at least once (positive statement).
  #   2. The phrase "does not capture" (or equivalent disclaimer) appears.
  expect_match(rd_text, "yhat-dependent")
  expect_match(rd_text, "x-driven|orthogonal to the mean")
})
