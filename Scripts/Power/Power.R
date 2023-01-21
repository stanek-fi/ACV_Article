library(ACV)
library(data.table)
library(ggplot2)
library(stringr)
library(doParallel)
rm(list = ls())
targetPath <- "Outputs"
suppressWarnings(dir.create(targetPath))
ArticleHeight <- 12
ArticleWidth <- 10
PresentationHeight <- 10
PresentationWidth <- 12

ms <- 100
ns <- c(10, 20, 50, 100, 300)
v <- 1
hs <- c(1, 3, 6)
R <- 2000
rhos <- c(1, 1.03125, 1.0625, 1.125, 1.25, 1.375, 1.5, 1.75, 2)
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
mres[, rho := as.numeric(rho)]
level <- 0.05
amres <- mres[, .(rejectionProb = mean(pval < level, na.rm = T), SE = sqrt(var(pval < level, na.rm = T) / sum(!is.na(pval)))), .(test, n, rho, h)]
amres[, testLong := ifelse(test %in% c("ADM", "DM"), "Diebold-Mariano", "Ibragimov-Muller")]
amres[, method := factor(ifelse(test %in% c("ADM", "AIM"), "ADM, AIM", "DM, IM"), levels = c("DM, IM", "ADM, AIM"))]
amres[, nLong := factor(n, levels = ns, labels = str_c("n=", ns))]
amres[, hLong := factor(h, levels = hs, labels = str_c("h=", hs))]

for (htemp in hs) {
  ggplot(amres[h == htemp], aes(x = rho, y = rejectionProb, colour = method)) +
    geom_hline(yintercept = level, linetype = "dotted") +
    geom_line(alpha = 0.5) +
    geom_errorbar(aes(ymin = rejectionProb + qnorm(0.025) * SE, ymax = rejectionProb + qnorm(0.975) * SE), width = 0.03) +
    scale_x_continuous(breaks = rhos) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    theme(panel.grid.minor = element_blank()) +
    facet_grid(nLong ~ testLong) +
    ylab("Rejection Probability") +
    xlab(expression(varsigma)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "bottom", legend.direction = "vertical") +
    labs(color = "Testing Method")
  ggsave(file.path(targetPath, str_c("Power_h=", htemp, "_Article.pdf")), height = ArticleHeight, width = ArticleWidth)
  ggsave(file.path(targetPath, str_c("Power_h=", htemp, "_Presentation.pdf")), height = PresentationHeight, width = PresentationWidth)
}
