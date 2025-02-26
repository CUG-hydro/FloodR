set.seed(1)
dailyMQ <- data.frame(
  Date = seq(
    from = as.Date("01.01.2000", format = "%d.%m.%Y"),
    to = as.Date("01.01.2004", format = "%d.%m.%Y"), by = "days"
  ),
  discharge = rbeta(1462, 2, 20) * 100
)

monthlyHQ <- data.frame(
  Date = seq(
    from = as.Date("01.01.2000", format = "%d.%m.%Y"),
    to = as.Date("01.01.2004", format = "%d.%m.%Y"), by = "months"
  ),
  discharge = dailyMQ$discharge[(0:48) * 12 + 1] + rnorm(49, 5, 1)
)

head(dailyMQ)
head(monthlyHQ)

r_old <- FloodR::eventsep(dailyMQ)
r_new <- eventsep(dailyMQ)
all.equal(r_old, r_new)


r2 = eventsep(dailyMQ, monthlyHQ)
r2
