# Skip helper for slow tests.
# By default (env var unset) heavy autotune/CV tests are skipped to keep
# `devtools::test()` under ~30s during iterative development. Set
# `ARES_FULL_TESTS=1` to run the full suite (CI, pre-release, CRAN prep).
skip_if_quick <- function() {
  if (!nzchar(Sys.getenv("ARES_FULL_TESTS"))) {
    testthat::skip("ARES_FULL_TESTS not set; skipping slow test")
  }
}
