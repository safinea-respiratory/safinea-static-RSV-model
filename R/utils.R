# Round a numeric vector to integers while preserving its total sum.
# Ensures the rounded marginals passed to r2dtable still sum to the same
# reconciled grand total despite per-element floor() truncation.
# The largest fractional remainders absorb any leftover units (largest-
# remainder method).
round_preserve_sum <- function(x) {
  sx   <- sum(x, na.rm = TRUE)
  rx   <- floor(x)
  diff <- as.integer(round(sx - sum(rx)))
  if (diff > 0) {
    o <- order(x - rx, decreasing = TRUE, na.last = NA)
    rx[o[seq_len(diff)]] <- rx[o[seq_len(diff)]] + 1L
  }
  rx
}
