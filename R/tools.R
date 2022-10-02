#' @import zeallot
na_interp_Q <- function(Q, NA_mode = NULL) {
  # Q = daten[, 2]
  n = length(Q)

  lgl_na <- is.na(Q) * 1 # true: 1, false: 0
  rl <- rle(lgl_na)
  to_int <- which(rl$lengths <= NA_mode & rl$values == 1)
  if (length(to_int) >= 1) {
    for (i in 1:length(to_int)) {
      i_beg = sum(rl$lengths[seq(1, to_int[i] - 1)])
      i_end = sum(rl$lengths[seq(1, to_int[i])]) + 1
      inds <- i_beg:i_end
      Q[inds] <- seq(Q[inds[1]], Q[inds[length(inds)]], length.out = length(inds))
    }
  }
  if (Messages) cat("Interpolated", length(to_int), "NA-gaps\n")

  lgl_na <- is.na(Q) * 1
  spl <- rep(NA, n)
  lastNA <- FALSE
  cc <- 1
  for (i in 1:n) {
    if (lgl_na[i] == 0 & lastNA == TRUE) {
      cc <- cc + 1
      spl[i] <- cc
      lastNA <- FALSE
    } else if (lgl_na[i] == 0 & lastNA == FALSE) {
      spl[i] <- cc
      lastNA <- FALSE
    } else {
      spl[i] <- NA
      lastNA <- TRUE
    }
  }
  list(Q = Q, spl = spl)
}

#' @import magrittr
cumsum2 <- function(x, ...) {
  if (any(is.na(x))) {
    total <- numeric(length(x))
    total[1] <- x[1]
    for (k in 2:length(x)) {
      if (is.na(x[k])) {
        total[k] <- 0
      } else {
        total[k] <- total[k - 1] + x[k]
      }
    }
  } else {
    total <- cumsum(x)
  }
  total
}
