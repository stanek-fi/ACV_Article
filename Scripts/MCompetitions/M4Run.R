# install.packages("https://github.com/carlanetto/M4comp2018/releases/download/0.2.0/M4comp2018_0.2.0.tar.gz",repos=NULL)
rm(list = ls())
library(M4comp2018)
library(ACV)
library(forecTheta)
library(data.table)
library(trend)
library(ggplot2)
library(doParallel)
library(seastests)
data(M4)

Periods <- sapply(sapply(M4, "[", "period"), "[", 1)
table(Periods)

SelectedPeriods <- c("Yearly", "Quarterly", "Monthly", "Weekly", "Daily", "Hourly")


# LossFunctions -----------------------------------------------------------

lossFunction <- function(y, yhat) {
  2 * abs(y - yhat) / (abs(y) + abs(yhat)) * 100
} # sMAPE
# lossFunction <- function(y, yhat) {
#   abs(y - yhat)
# } # sMAPE
# lossFunction <- function(y, yhat) {
#   (y - yhat)^2
# } # MSE

# Algorithms --------------------------------------------------------------

algorithm1 <- function(y) {
  ets(y, model = "ZZZ")
}


algorithm2 <- function(y) {
  auto.arima(y, approximation = T, max.P = 1, max.Q = 1)
}

rhoLimit <- 0.99

# MainLoop ----------------------------------------------------------------
set.seed(1)
SelectedPeriod <- SelectedPeriods[1]
for (SelectedPeriod in SelectedPeriods) {
  Ids <- which(Periods == SelectedPeriod)
  print(length(Ids))

  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(M4comp2018))
  clusterEvalQ(cl, library(ACV))
  clusterEvalQ(cl, library(forecTheta))
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(trend))
  clusterEvalQ(cl, library(seastests))
  clusterExport(cl, c("algorithm1", "algorithm2", "lossFunction"))


  IdsChunks <- split(Ids, ceiling(seq_along(Ids) / length(Ids) * min(100, ceiling(length(Ids) / numCores))))
  resrows <- vector(mode = "list", length = length(IdsChunks))

  Chunk <- 1
  i <- 1
  start <- Sys.time()
  for (Chunk in seq_along(IdsChunks)) {
    M4Subset <- M4[IdsChunks[[Chunk]]]
    clusterExport(cl, c("M4Subset"))

    resrows[[Chunk]] <- do.call(rbind, foreach(i = seq_along(IdsChunks[[Chunk]])) %dopar% {
      resrow <- data.table(
        Chunk = Chunk,
        i = i,
        Id = IdsChunks[[Chunk]][i],
        Period = as.character(NA),
        LengthObs = as.numeric(NA),
        LengthOut = as.numeric(NA),
        TrendRsq = as.numeric(NA),
        TrendPval = as.numeric(NA),
        TrendTval = as.numeric(NA),
        TrendCSPval = as.numeric(NA),
        SeasonalQSPval = as.numeric(NA),
        EstLCV1 = as.numeric(NA),
        EstLACV1 = as.numeric(NA),
        L1 = as.numeric(NA),
        Rho1 = as.numeric(NA),
        AIC1 = as.numeric(NA),
        AICc1 = as.numeric(NA),
        BIC1 = as.numeric(NA),
        EstLCV2 = as.numeric(NA),
        EstLACV2 = as.numeric(NA),
        L2 = as.numeric(NA),
        Rho2 = as.numeric(NA),
        AIC2 = as.numeric(NA),
        AICc2 = as.numeric(NA),
        BIC2 = as.numeric(NA),
        RhoDiff = as.numeric(NA),
        SelectedACV = as.numeric(NA),
        SelectedACVNaive = as.numeric(NA),
        SelectedCV = as.numeric(NA),
        SelectedAIC = as.numeric(NA),
        SelectedAICc = as.numeric(NA),
        SelectedBIC = as.numeric(NA),
        SelectedOptimal = as.numeric(NA)
      )

      try({

        # Preparation and Statistics ----------------------------------------------

        yObs <- M4Subset[[i]]$x
        yOut <- M4Subset[[i]]$xx
        yAll <- ts(c(yObs, yOut), start = start(yObs), frequency = frequency(yObs))
        location <- mean(yAll)
        scale <- var(yAll)^(.5)
        yAll <- (yAll - location) / scale
        mAll <- length(yObs)
        mObs <- length(yObs) - length(yOut)
        yObs <- window(yAll, time(yAll)[1], time(yAll)[mAll])
        yOut <- window(yAll, time(yAll)[mAll + 1], time(yAll)[length(yAll)])
        h <- length(yOut)
        v <- h
        mIdentical <- T

        resrow[1, Period := as.vector(M4Subset[[i]]$period)]
        resrow[1, LengthObs := length(yObs)]
        resrow[1, LengthOut := length(yOut)]
        TrendModel <- summary(lm(y ~ trend, data.frame(y = yAll, trend = seq_along(yAll))))
        resrow[1, TrendRsq := TrendModel$r.squared]
        resrow[1, TrendPval := TrendModel$coefficients["trend", "Pr(>|t|)"]]
        resrow[1, TrendTval := TrendModel$coefficients["trend", "t value"]]
        resrow[1, TrendCSPval := cs.test(yAll)$p.value]
        resrow[1, SeasonalQSPval := qs(yObs)$Pval]

        # Algorithm 1 ----------------------------------------------

        m <- mObs
        PhiObs1 <- tsACV(yObs, algorithm1, m = m, h = h, v = v, lossFunction = lossFunction)
        Est <- estimateL(Phi = PhiObs1, method = "regular")
        resrow[1, EstLCV1 := Est$estimate]
        Est <- estimateL(Phi = PhiObs1, method = "augmented", rhoLimit = rhoLimit)
        resrow[1, EstLACV1 := Est$estimate]
        resrow[1, Rho1 := Est$rho]
        model <- algorithm1(yObs)
        resrow[1, AIC1 := model$aic]
        resrow[1, AICc1 := model$aicc]
        resrow[1, BIC1 := model$bic]

        if (mIdentical) {
          m <- mObs
          PhiAll1 <- tsACV(window(yAll, time(yAll)[length(yOut) + 1], time(yAll)[length(yAll)]), algorithm1, m = m, h = h, v = v, lossFunction = lossFunction)
        } else {
          m <- mAll
          PhiAll1 <- tsACV(yAll, algorithm1, m = m, h = h, v = v, lossFunction = lossFunction)
        }
        resrow[1, L1 := estimateL(Phi = PhiAll1, method = "regular")$estimate]

        # Algorithm 2 ----------------------------------------------

        m <- mObs
        PhiObs2 <- tsACV(yObs, algorithm2, m = m, h = h, v = v, lossFunction = lossFunction)
        Est <- estimateL(Phi = PhiObs2, method = "regular")
        resrow[1, EstLCV2 := Est$estimate]
        Est <- estimateL(Phi = PhiObs2, method = "augmented", rhoLimit = rhoLimit)
        resrow[1, EstLACV2 := Est$estimate]
        resrow[1, Rho2 := Est$rho]
        model <- algorithm2(yObs)
        resrow[1, AIC2 := model$aic]
        resrow[1, AICc2 := model$aicc]
        resrow[1, BIC2 := model$bic]

        if (mIdentical) {
          m <- mObs
          PhiAll2 <- tsACV(window(yAll, time(yAll)[length(yOut) + 1], time(yAll)[length(yAll)]), algorithm2, m = m, h = h, v = v, lossFunction = lossFunction)
        } else {
          m <- mAll
          PhiAll2 <- tsACV(yAll, algorithm2, m = m, h = h, v = v, lossFunction = lossFunction)
        }
        resrow[1, L2 := estimateL(Phi = PhiAll2, method = "regular")$estimate]

        # Algorithm Difference  ----------------------------------------------

        PhiObsDiff <- PhiObs1 - PhiObs2
        Est <- estimateL(Phi = PhiObsDiff, method = "augmented", rhoLimit = rhoLimit)
        resrow[1, SelectedACV := ifelse(Est$estimate > 0, 2, 1)]
        resrow[1, RhoDiff := Est$rho]
        resrow[1, SelectedACVNaive := ifelse(EstLACV1 - EstLACV2 > 0, 2, 1)]
        resrow[1, SelectedCV := ifelse(EstLCV1 - EstLCV2 > 0, 2, 1)]
        resrow[1, SelectedAIC := ifelse(AIC1 - AIC2 > 0, 2, 1)]
        resrow[1, SelectedAICc := ifelse(AICc1 - AICc2 > 0, 2, 1)]
        resrow[1, SelectedBIC := ifelse(BIC1 - BIC2 > 0, 2, 1)]
        resrow[1, SelectedOptimal := ifelse(L1 - L2 > 0, 2, 1)]
      })
      return(resrow)
    })

    print(paste0("Period: ", SelectedPeriod, " ", Chunk / length(IdsChunks) * 100, "% time:", Sys.time()))
  }
  Sys.time() - start
  stopCluster(cl)
  res <- do.call(rbind, resrows)
  saveRDS(res, file.path("Outputs", paste0(SelectedPeriod, ".RDS")))

  print(apply(res, 2, function(x) mean(is.na(x))))
  print(res[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = c("L1", "L2")])
  print(res[, lapply(.SD, function(x) mean((x - L1)^2, na.rm = T)), .SDcols = c("EstLCV1", "EstLACV1")][, .(EstLCV1, EstLACV1, ratio = EstLACV1 / EstLCV1)])
  print(res[, lapply(.SD, function(x) mean((x - L2)^2, na.rm = T)), .SDcols = c("EstLCV2", "EstLACV2")][, .(EstLCV2, EstLACV2, ratio = EstLACV2 / EstLCV2)])
  print(res[, lapply(.SD, function(x) c(mean(x - L1, na.rm = T), var(x, na.rm = T))), .SDcols = c("EstLCV1", "EstLACV1")])
  print(res[, lapply(.SD, function(x) c(mean(x - L2, na.rm = T), var(x, na.rm = T))), .SDcols = c("EstLCV2", "EstLACV2")])
  res[, Diff1 := (L1 - EstLACV1)^2 - (L1 - EstLCV1)^2]
  print(res[, mean(Diff1, na.rm = T) / (sd(Diff1, na.rm = T) / sqrt(.N))])
  res[, Diff2 := (L2 - EstLACV2)^2 - (L2 - EstLCV2)^2]
  print(res[, mean(Diff2, na.rm = T) / (sd(Diff2, na.rm = T) / sqrt(.N))])


  cols <- c("SelectedACV", "SelectedCV", "SelectedAIC", "SelectedAICc", "SelectedBIC", "SelectedOptimal")
  out <- res[, lapply(.SD, function(x) mean(ifelse(x == 2, L2, L1), na.rm = T)), .SDcols = cols]
  print(out)
  print(out / out$SelectedOptimal)

  out <- res[, lapply(.SD, function(x) mean(x == SelectedOptimal, na.rm = T)), .SDcols = cols]
  print(out)

  print(res[, lapply(.SD, function(x) mean(x, na.rm = T)), .SDcols = cols])
}