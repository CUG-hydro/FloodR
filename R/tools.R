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

# W: 1e6 m^3/d
Q2W <- function(Q) {
  Q * 86400 / 1e6
}

cal_floodInfo <- function(q, t, pos_start, pos_peak, pos_end, 
  monthlyHQ = NULL, comm = "") 
{
  varnames <- c(
    "Begin", "End", "Peak_date", "DailyMQ", "Volume", "dir_Volume",
    "baseflow_peak", "baseflow_begin", "baseflow_end", "No_Peaks", "HQ", "HQ_dir", "Comments"
  )
  
  basefl        <- approxfun(c(t[pos_start], t[pos_end]), c(q[pos_start], q[pos_end]))
  
  start         <- t[pos_start]
  end           <- t[pos_end]
  peak_MQ       <- max(q[pos_start:pos_end])
  peak_date     <- t[pos_peak]
  Volume        <- sum(q[pos_start:pos_end]) %>% Q2W()
  Baseflow_peak <- basefl(peak_date)
  Base_vol      <- integrate(basefl, lower = t[pos_start], upper = t[pos_end])$value %>% Q2W()
  Vol_dir       <- Volume - Base_vol
  base_start    <- basefl(start)
  base_end      <- basefl(end)

  HQ <- get_HQ(monthlyHQ, peak_date, peak_MQ)
  HQ_dir <- HQ - Baseflow_peak

  data.frame(
    start, end, peak_date, peak_MQ, Volume,
    Vol_dir, Baseflow_peak, base_start, base_end, no_peaks, HQ, HQ_dir, comm
  ) %>% set_names(varnames) # info
}

# get_HQ(monthlyHQ, peak_date, peak_MQ)
get_HQ <- function(monthlyHQ, peak_date, peak_MQ) {
  HQ <- monthlyHQ[which(
    (monthlyHQ[, 1] >= (peak_date - 1)) &
    (monthlyHQ[, 1] <= (peak_date + 1))
  ), 2]
  HQ = ifelse(length(HQ) < 1, NA, max(HQ))
  if (!is.na(HQ) && (HQ < peak_MQ)) HQ <- NA
  HQ
}

# c(base_diff, base_rel) %<-% get_base_rel(t, q, pos_start, pos_end)
get_base_rel <- function(t, q, pos_start, pos_end) {
  basefl1 <- approxfun(c(t[pos_start], t[pos_end]), c(q[pos_start], q[pos_end]))
  base_diff <- (basefl1(t[pos_start:pos_end]) - q[pos_start:pos_end])
  base_rel <- (basefl1(t[pos_end]) - basefl1(t[pos_start])) / (pos_end - pos_start)
  c(base_diff, base_rel)
}

