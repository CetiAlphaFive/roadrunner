# BUG-004 — OOV factor level silently drops rows from predict() output

**Severity:** robustness (medium).
**File:** `R/predict.R:99-117`.

## Symptom

When `newdata` contains a factor (or character) column with a level not
seen during training, `predict()` returns FEWER rows than `nrow(newdata)`
with no warning.

## Root cause

```r
for (jname in names(fi$xlevels)) {
  col <- newdata[[jname]]
  if (is.character(col) || is.factor(col)) {
    newdata[[jname]] <- factor(col, levels = fi$xlevels[[jname]])
    # ^ OOV values become NA
  }
}
xnew <- stats::model.matrix(~ ., data = newdata)
# model.matrix default na.action = na.omit -> silently drops rows with NA in
# any factor column.
```

The code at predict.R:111-115 already detects missing *columns* and errors,
but `model.matrix` drops *rows* with NA factor levels before that check
runs.

## Repro

```r
library(roadrunner)
set.seed(1)
df <- data.frame(num = runif(50),
                 cat = factor(sample(c("aaa", "bb-bb"), 50, replace = TRUE)))
yy <- df$num + (df$cat == "bb-bb") + rnorm(50, sd = 0.3)
fit <- ares(df, yy, nthreads = 2)

# Now simulate OOV in newdata:
df2 <- df
df2$cat <- as.character(df2$cat)
df2$cat[1:5] <- "OOV"      # 5 rows with unseen level

p <- predict(fit, df2)
length(p)   # 45, not 50 -- SILENT DROP, no warning
```

## Impact

Calling code that assumes `length(predict(fit, newdata)) == nrow(newdata)`
will get a length mismatch with no indication. Downstream alignment
errors (e.g. `data.frame(newdata, pred = predict(fit, newdata))`) crash
or silently misalign rows.

## Suggested fix

Before `model.matrix(~ ., newdata)`, detect OOV rows explicitly:

```r
oov_rows <- integer(0)
for (jname in names(fi$xlevels)) {
  col <- newdata[[jname]]
  if (is.character(col) || is.factor(col)) {
    new_col <- factor(col, levels = fi$xlevels[[jname]])
    oov_rows <- union(oov_rows, which(is.na(new_col)))
    newdata[[jname]] <- new_col
  }
}
if (length(oov_rows)) {
  warning("ares: ", length(oov_rows), " row(s) with out-of-vocabulary",
          " factor levels; predictions for those rows will be NA.",
          call. = FALSE)
  # ... return a length-nrow(newdata) prediction with NA at oov_rows
}
```

Then call `model.matrix` with `na.action = na.pass` so the row count
stays consistent, and fill OOV-row predictions with NA at the end.
