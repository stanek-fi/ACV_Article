library(ACV)
library(data.table)
library(ggplot2)
library(stringr)
library(doParallel)
rm(list = ls())
targetPath <- "Outputs"
suppressWarnings(dir.create(targetPath))



ms <- 100
ns <- c(10, 20, 50, 100, 300)
v <- 1
hs <- c(1, 3, 6)
R <- 10000
rhos <- 1
Ha <- "!=0"
sigma <- 1
decay <- 0.5

algorithm1 <- function(yInSample, yOutSample, h, omitLast = 0) {
  return(list(
    yhatInSample = rep(0, length(yInSample)),
    yhatOutSample = rep(0, length(yOutSample))
  ))
}

algorithm2 <- function(yInSample, yOutSample, h, omitLast = 0) {
  meanEst <- mean(yInSample[1:(length(yInSample) - omitLast)])
  return(list(
    yhatInSample = rep(meanEst, length(yInSample)),
    yhatOutSample = rep(meanEst, length(yOutSample))
  ))
}

tests <- c("DM", "ADM", "IM", "AIM")
res <- array(NA, dim = c(R, length(hs), length(ms), length(ns), length(rhos), length(tests)))
dimnames(res) <- list(r = 1:R, h = hs, m = ms, n = ns, rho = rhos, test = tests)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)
clusterEvalQ(cl, library(ACV))
start <- Sys.time()
for (h in hs) {
  for (m in ms) {
    for (n in ns) {
      for (rho in rhos) {
        temp <- foreach(r = 1:R) %dopar% {
          set.seed(r)

          thetas <- sapply(seq_len(h) - 1, function(j) decay^(j))
          gammas <- sapply(seq_len(h) - 1, function(j) sigma^2 * sum(thetas[1:(h - j)] * thetas[(1 + j):h]))
          const <- sqrt(rho * (gammas[1] + 1 / m * gammas[1] + 2 * sum(do.call(c, lapply(seq_len(h - 1), function(j) (m - j) / m^2 * gammas[j + 1])))) - gammas[1])
          y <- arima.sim(n = n + m + h - 1, list(ma = thetas[-1]), sd = sigma) + const
          Phi1 <- tsACV(y, algorithm1, m + (h - 1), h = 1, v = 1, omitLast = h - 1)
          Phi2 <- tsACV(y, algorithm2, m + (h - 1), h = 1, v = 1, omitLast = h - 1)
          Phi <- Phi1 - Phi2

          list(
            DM = testL(Phi = Phi, method = "regular", test = "Diebold-Mariano", Ha = Ha)$pval,
            ADM = testL(Phi = Phi, method = "augmented", test = "Diebold-Mariano", Ha = Ha)$pval,
            IM = testL(Phi = Phi, method = "regular", test = "Ibragimov-Muller", Ha = Ha, groups = 2)$pval,
            AIM = testL(Phi = Phi, method = "augmented", test = "Ibragimov-Muller", Ha = Ha, groups = 2)$pval
          )
        }
        res[, as.character(h), as.character(m), as.character(n), as.character(rho), "DM"] <- sapply(temp, "[[", "DM")
        res[, as.character(h), as.character(m), as.character(n), as.character(rho), "ADM"] <- sapply(temp, "[[", "ADM")
        res[, as.character(h), as.character(m), as.character(n), as.character(rho), "IM"] <- sapply(temp, "[[", "IM")
        res[, as.character(h), as.character(m), as.character(n), as.character(rho), "AIM"] <- sapply(temp, "[[", "AIM")

        print(str_c("h=", h, " m=", m, ", n=", n, ", rho=", rho, ", time:", Sys.time()))
      }
    }
  }
}
Sys.time() - start
stopCluster(cl)

saveRDS(res, "res.RDS")
res <- readRDS("res.RDS")

mres <- as.data.table(res, value.name = "pval")
levels <- c(0.01, 0.05, 0.1)

rp <- function(x, p) {
  formatC(round(x, p), format = "f", digits = p)
}
p <- 3
amres <- mres[, .(rejectionProb_001 = rp(mean(pval < 0.01, na.rm = T), p), rejectionProb_005 = rp(mean(pval < 0.05, na.rm = T), p), rejectionProb_010 = rp(mean(pval < 0.10, na.rm = T), p)), .(test, n, h)]
amres[, h := factor(h, levels = hs)]
amres[, n := factor(n, levels = ns)]
amres[, test := factor(test, levels = tests)]
out <- dcast(amres, h + n ~ test, value.var = c("rejectionProb_001", "rejectionProb_005", "rejectionProb_010"))
out[, h := do.call(c, lapply(hs, function(h) c(h, rep("", length(ns) - 1))))]
write.csv(out, file.path(targetPath, str_c("Level_001005010.csv")), quote = F, row.names = F)
