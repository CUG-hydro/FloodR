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

nansum <- function(x, ...) {
  sum(x, ..., na.rm = TRUE)
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

# events_temp = cal_floodInfo(q, t, pos_peak, pos_start, pos_end,
#   monthlyHQ, comm = "")
cal_floodInfo <- function(q, t, pos_peak, pos_start, pos_end, 
  no_peaks, monthlyHQ = NULL, comm = "", ...,
  t_start = NULL, t_end = NULL) 
{
  varnames <- c(
    "Begin", "End", "Peak_date", "DailyMQ", "Volume", "dir_Volume",
    "baseflow_peak", "baseflow_begin", "baseflow_end", "No_Peaks", "HQ", "HQ_dir", "Comments"
  )
  
  if (is.null(t_start)) t_start <- t[pos_start]
  if (is.null(t_end)) t_end <- t[pos_end]
  
  basefl <- approxfun(c(t_start, t_end), c(q[pos_start], q[pos_end]))
  
  peak_MQ       <- max(q[pos_start:pos_end])
  t_peak        <- t[pos_peak]
  Volume        <- sum(q[pos_start:pos_end]) %>% Q2W()
  Baseflow_peak <- basefl(t_peak)
  Base_vol      <- integrate(basefl, lower = t_start, upper = t_end)$value %>% Q2W()
  Vol_dir       <- Volume - Base_vol
  base_start    <- basefl(t_start)
  base_end      <- basefl(t_end)

  HQ <- get_HQ(monthlyHQ, t_peak, peak_MQ)
  HQ_dir <- HQ - Baseflow_peak

  data.frame(
    t_start, t_end, t_peak, peak_MQ, Volume,
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

  ind = pos_start:pos_end
  base_diff <- (basefl1(t[ind]) - q[ind])
  base_rel <- (basefl1(t[pos_end]) - basefl1(t[pos_start])) / (pos_end - pos_start)
  
  list(base_diff, base_rel)
}


# c(no_max, no_peaks, pos_peak) %<-% get_no_peaks(q, var3d, thvar, pos_start, pos_end)
get_no_peaks <- function(q, var3d, thvar, pos_start, pos_end) {
  no_max <- diff(sign(diff(var3d[pos_start:pos_end])))
  no_peaks <- max(1, sum((no_max < -1) &
    (var3d[(pos_start + 1):(pos_end - 1)] > thvar) &
    (diff(sign(diff(q[pos_start:pos_end]))) < -1)))
  pos_peak <- which.max(q[pos_start:pos_end]) + pos_start - 1
  list(no_max, no_peaks, pos_peak)
}
