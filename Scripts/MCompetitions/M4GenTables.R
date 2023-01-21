rm(list = ls())
library(M4comp2018)
library(ACV)
library(forecTheta)
library(data.table)
library(trend)
library(ggplot2)
library(doParallel)
library(seastests)
library(stringr)

targetPath <- "Tables"
FileName <- file.path("Outputs")
SelectedPeriods <- c("Yearly", "Quarterly", "Monthly", "Weekly", "Daily", "Hourly")

rp <- function(x, p) {
  formatC(round(x, p), format = "f", digits = p)
}
stars <- function(pval) {
  if (is.null(pval) | is.na(pval)) {
    "   "
  } else if (pval > 0.05) {
    "   "
  } else if (pval > 0.01) {
    "*  "
  } else if (pval > 0.001) {
    "** "
  } else {
    "***"
  }
}

hs <- c(1, 3, 6)
timingTables <- setNames(c(F, T, T), hs)
h <- 1
for (h in hs) {
  res <- do.call(rbind, lapply(SelectedPeriods, function(SelectedPeriod) {
    readRDS(file.path(FileName, paste0(SelectedPeriod, "_", h, ".RDS")))
  }))
  names(res)

  res[, Trending := ifelse(as.logical(TrendCSPval < 0.05), "T", "F")]
  res[, Seasonal := ifelse(as.logical(SeasonalQSPval < 0.05), "T", "F")]
  res[, ECV1 := (L1 - EstLCV1)^2]
  res[, EACV1 := (L1 - EstLACV1)^2]
  res[, DE1 := EACV1 - ECV1]
  res[, ECV2 := (L2 - EstLCV2)^2]
  res[, EACV2 := (L2 - EstLACV2)^2]
  res[, DE2 := EACV2 - ECV2]
  res[, LossOptimal := ifelse(SelectedOptimal == 2, L2, L1)]
  res[, LossCV := ifelse(SelectedCV == 2, L2, L1)]
  res[, LossACV := ifelse(SelectedACV == 2, L2, L1)]
  res[, LossAIC := ifelse(SelectedAIC == 2, L2, L1)]
  res[, DLossCV := LossCV - LossOptimal]
  res[, DLossACV := LossACV - LossOptimal]
  res[, DLossAIC := LossAIC - LossOptimal]

  p <- 2
  fun <- function(SD) {
    SD[, .(
      N = c(str_c(.N), ""),
      MSECV1 = c(rp(mean(ECV1, na.rm = T), p), str_c("(", rp(sd(ECV1, na.rm = T) / sqrt(.N), p), ")")),
      MSEACV1 = c(rp(mean(EACV1, na.rm = T), p), str_c("(", rp(sd(EACV1, na.rm = T) / sqrt(.N), p), ")")),
      ratio1 = c(str_c(rp((mean(EACV1, na.rm = T) / mean(ECV1, na.rm = T) - 1) * 100, p - 1), stars(2 * (1 - pnorm(abs(mean(DE1, na.rm = T) / (sd(DE1, na.rm = T) / sqrt(.N))))))), ""),
      MSECV2 = c(rp(mean(ECV2, na.rm = T), p), str_c("(", rp(sd(ECV2, na.rm = T) / sqrt(.N), p), ")")),
      MSEACV2 = c(rp(mean(EACV2, na.rm = T), p), str_c("(", rp(sd(EACV2, na.rm = T) / sqrt(.N), p), ")")),
      ratio2 = c(str_c(rp((mean(EACV2, na.rm = T) / mean(ECV2, na.rm = T) - 1) * 100, p - 1), stars(2 * (1 - pnorm(abs(mean(DE2, na.rm = T) / (sd(DE2, na.rm = T) / sqrt(.N))))))), "")
    ), ]
  }
  out3 <- res[, fun(.SD), .(Period, Trending, Seasonal)][order(Period, Trending, Seasonal)]
  out2 <- res[, c(Trending = "", Seasonal = "", fun(.SD)), .(Period)]
  out1 <- res[, c(Period = "All", Trending = "", Seasonal = "", fun(.SD)), ]
  out <- rbind(
    do.call(rbind, lapply(SelectedPeriods, function(SelectedPeriod) {
      rbind(out2[Period == SelectedPeriod], out3[Period == SelectedPeriod], fill = T)
    })),
    out1
  )
  out[, Period := ifelse(Trending == "" & Seasonal == "" & N != "", Period, "")]
  out[, Trending := ifelse(N != "", Trending, "")]
  out[, Seasonal := ifelse(N != "", Seasonal, "")]
  write.csv(out, file.path(targetPath, str_c("M4_MSE_", h, ".csv")), quote = F, row.names = F)


  p <- 3
  fun <- function(SD) {
    SD[, .(
      N = c(str_c(.N), ""),
      SelectedOptimal = c(rp(mean(SelectedOptimal == SelectedOptimal, na.rm = T), p), str_c("(", rp(sd(SelectedOptimal == SelectedOptimal, na.rm = T) / sqrt(.N), p + 1), ")")),
      LossOptimal = c(rp(mean(LossOptimal, na.rm = T), p), str_c("(", rp(sd(LossOptimal, na.rm = T) / sqrt(.N), p), ")")),
      SelectedAIC = c(rp(mean(SelectedAIC == SelectedOptimal, na.rm = T), p), str_c("(", rp(sd(SelectedAIC == SelectedOptimal, na.rm = T) / sqrt(.N), p), ")")),
      LossAIC = c(rp(mean(LossAIC, na.rm = T), p), str_c("(", rp(sd(LossAIC, na.rm = T) / sqrt(.N), p), ")")),
      SelectedCV = c(rp(mean(SelectedCV == SelectedOptimal, na.rm = T), p), str_c("(", rp(sd(SelectedCV == SelectedOptimal, na.rm = T) / sqrt(.N), p), ")")),
      LossCV = c(rp(mean(LossCV, na.rm = T), p), str_c("(", rp(sd(LossCV, na.rm = T) / sqrt(.N), p), ")")),
      SelectedACV = c(rp(mean(SelectedACV == SelectedOptimal, na.rm = T), p), str_c("(", rp(sd(SelectedACV == SelectedOptimal, na.rm = T) / sqrt(.N), p), ")")),
      LossACV = c(rp(mean(LossACV, na.rm = T), p), str_c("(", rp(sd(LossACV, na.rm = T) / sqrt(.N), p), ")")),
      LossRelGainAIC = c(str_c(rp((mean(DLossACV, na.rm = T) / mean(DLossAIC, na.rm = T) - 1) * 100, p - 2), stars(2 * (1 - pnorm(abs(mean(LossACV - LossAIC, na.rm = T) / (sd(LossACV - LossAIC, na.rm = T) / sqrt(.N))))))), ""),
      LossRelGainCV = c(str_c(rp((mean(DLossACV, na.rm = T) / mean(DLossCV, na.rm = T) - 1) * 100, p - 2), stars(2 * (1 - pnorm(abs(mean(LossACV - LossCV, na.rm = T) / (sd(LossACV - LossCV, na.rm = T) / sqrt(.N))))))), "")
    ), ]
  }
  out2 <- res[, fun(.SD), Period]
  out1 <- res[, c(Period = "All", fun(.SD))]
  out <- rbind(
    do.call(rbind, lapply(SelectedPeriods, function(SelectedPeriod) {
      out2[Period == SelectedPeriod]
    })),
    out1
  )
  out[, Period := ifelse(N != "", Period, "")]
  out <- out[, -3]
  write.csv(out, file.path(targetPath, str_c("M4_Selection_", h, ".csv")), quote = F, row.names = F)

  if (timingTables[as.character(h)]) {
    out2 <- res[, lapply(.SD, function(x) rp(mean(x), p)), Period, .SDcols = c("LengthObs", "LengthOut", "TimePhiObs1", "TimeEstLCV1", "TimeEstLACV1", "TimePhiObs2", "TimeEstLCV2", "TimeEstLACV2")]
    out1 <- res[, c(Period = "All", lapply(.SD, function(x) rp(mean(x), p))), .SDcols = c("LengthObs", "LengthOut", "TimePhiObs1", "TimeEstLCV1", "TimeEstLACV1", "TimePhiObs2", "TimeEstLCV2", "TimeEstLACV2")]
    out <- rbind(out2, out1)
    write.csv(out, file.path(targetPath, str_c("M4_Timing_", h, ".csv")), quote = F, row.names = F)
  }

  print(h)
  print(res[, .N / 100000, .(Trending, Seasonal)][(Trending == "T") | (Seasonal == "T"), sum(V1)])
  print(res[, .N / 100000, .(Trending, Seasonal)][(Trending == "T"), sum(V1)])
  print(res[, .N / 100000, .(Trending, Seasonal)][(Seasonal == "T"), sum(V1)])
  print(res[, .N / 100000, .(Trending, Seasonal)][(Trending == "T") & (Seasonal == "T"), sum(V1)])
}
