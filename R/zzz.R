# Package load hooks. We do not change global thread defaults at load time —
# `ares()` sets RcppParallel thread options per-call.
.onLoad <- function(libname, pkgname) {
  invisible()
}
