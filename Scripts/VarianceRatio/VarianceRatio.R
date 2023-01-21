library(ACV)
library(data.table)
library(ggplot2)
library(stringr)
library(doParallel)
library(latex2exp)

rm(list = ls())
targetPath <- "Outputs"
suppressWarnings(dir.create(targetPath))
ArticleHeight <- 7
ArticleWidth <- 10
PresentationHeight <- 10
PresentationWidth <- 12

algorithm <- function(yInSample, yOutSample, h) {
  xInSample <- as.matrix(yInSample[1:(length(yInSample) - 1)])
  model <- .lm.fit(x = xInSample, y = yInSample[2:length(yInSample)])
  return(list(
    yhatInSample = c(0, as.vector(xInSample %*% model$coefficients)),
    yhatOutSample = as.vector(as.matrix(c(yInSample[length(yInSample)], yOutSample[1:(length(yOutSample) - 1)])) %*% model$coefficients)
  ))
}

algorithmContrasts <- function(yInSample, yOutSample) {
  temp <- algorithm(yInSample, yOutSample, NA)
  (c(temp$yhatInSample, temp$yhatOutSample) - c(yInSample, yOutSample))^2
}

reconstructPhi <- function(x, v) {
  is <- seq(1, ncol(x), v)
  sapply(is, function(i) {
    c(rep(NA, i - 1), x[, i], rep(NA, max(is) - i))
  })[1:(nrow(x) - v + v * (length(is) - 1)), ]
}

TInf <- 1e+5
ms <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
ns <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300)
vs <- ns
R <- 1000
rhoLimit <- 0.99

m <- ms[1]
res <- do.call(rbind, lapply(ms, function(m) {
  iInf <- TInf - (m + max(vs))
  set.seed(1)
  yInf <- arima.sim(list(ar = c(0.9)), TInf, sd = 1)
  messageBreaks <- round(seq(0, iInf, length.out = 11))
  C <- matrix(1, nrow = m + max(vs), ncol = iInf + 1)
  for (i in 0:iInf) {
    if (i %in% messageBreaks) {
      print(str_c("m=", m, ", C matrix at ", (which(i == messageBreaks) - 1) / (length(messageBreaks) - 1) * 100, "%, time:", Sys.time()))
    }
    yInSample <- yInf[1:m + i]
    yOutSample <- yInf[(m + 1):(m + max(vs)) + i]
    C[, i + 1] <- algorithmContrasts(yInSample, yOutSample)
  }


  Cdm <- C - apply(C, 1, mean)
  VphiBlocks <- lapply(c(0, ns), function(di) {
    print(str_c("di=", di, ", time:", Sys.time()))
    Cdm[, 0:(iInf - di) + 1] %*% t(Cdm[, (0 + di):iInf + 1]) / length(0:(iInf - di) + 1)
  })

  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(ACV))

  ni <- 2
  resns <- vector("list", length(ns))
  for (ni in seq_along(ns)) {
    n <- ns[ni]
    v <- n
    Vphi <- rbind(
      cbind(VphiBlocks[[1]][1:(m + v), 1:(m + v)], VphiBlocks[[1 + which(vs == v)]][1:(m + v), 1:(m + v)]),
      cbind(t(VphiBlocks[[1 + which(vs == v)]][1:(m + v), 1:(m + v)]), VphiBlocks[[1]][1:(m + v), 1:(m + v)])
    )[1:(2 * m + v), 1:(2 * m + v)]
    B <- cbind(diag(m + v), diag(m + v)[, 1:m])
    b <- c(rep(0, m), rep(1 / v, v))
    Vphii <- solve(Vphi)
    lambdaOptimal <- as.vector(Vphii %*% t(B) %*% solve(B %*% Vphii %*% t(B)) %*% b)
    stepsize <- floor((iInf - n - 1) / (R - 1))

    clusterExport(cl, c("m", "reconstructPhi", "C", "rhoLimit"))
    res <- do.call(rbind, foreach(r = 1:R) %dopar% {
      res <- data.frame(n = rep(n, 1), m = rep(m, 1), LCV = rep(NA, 1), LACV = rep(NA, 1), LACVopt = rep(NA, 1))
      # Phi <- reconstructPhi(C[1:(m + v), 1:(n + 1) + (r - 1) * (n + 1)], v)
      Phi <- reconstructPhi(C[1:(m + v), 1:(n + 1) + (r - 1) * stepsize], v)
      estCV <- estimateL(Phi = Phi, method = "regular")
      res[1, "LCV"] <- estCV$estimate
      estACV <- estimateL(Phi = Phi, method = "augmented", rhoLimit = rhoLimit)
      res[1, "LACV"] <- estACV$estimate
      res[1, "LACVopt"] <- sum(na.omit(c(Phi)) * lambdaOptimal)
      return(res)
    })

    print(str_c("m=", m, ", n=", n, ", time:", Sys.time()))
    resns[[ni]] <- res
  }
  res <- as.data.table(do.call(rbind, resns))

  stopCluster(cl)
  print(res[, .(ratio = var(LACV) / var(LCV), ratioopt = var(LACVopt) / var(LCV), ratiocheck = var(LACVopt) / var(LACV), .N), .(m, n)])
  return(res)
}))

saveRDS(res, str_c("res.RDS"))

mres <- res[, .(varRatio = var(LACV) / var(LCV), optimality = (var(LCV) - var(LACV)) / (var(LCV) - var(LACVopt))), .(m, n)]
rules <- do.call(rbind, lapply(c(3, 5, 10), function(r) {
  data.table(
    rule = str_c("1/", r),
    m = ms,
    n = ms / (r - 1)
  )
}))
rules[, rule := factor(rule, levels = unique(rule))]

ggplot(mres, aes(x = m, y = n)) +
  geom_tile(aes(fill = varRatio), alpha = 1) +
  geom_text(aes(label = str_c(format(round(varRatio, 3), nsmall = 3), "\n", "{", format(round(optimality, 2), nsmall = 2), "}")), alpha = 0.5, size = 3.5) +
  scale_fill_gradientn(colours = terrain.colors(50), limits = c(0, 1.0001), name = TeX("$Var(\\widehat{L}_{ACV}) \\, / \\, Var(\\widehat{L}_{CV})$")) +
  scale_x_continuous(name = "m", breaks = ms) +
  scale_y_continuous(name = "n", breaks = ns) +
  geom_line(data = rules[n >= min(ns)], aes(x = m, y = n, linetype = rule), alpha = 0.5) +
  geom_point(data = rules[n >= min(ns)], aes(x = m, y = n, shape = rule), alpha = 0.5) +
  guides(
    fill = guide_colourbar(barwidth = 15, barheight = 1, title.position = "top", title.hjust = 0.5),
    linetype = guide_legend(keywidth = 4, keyheight = 1.3, title.position = "top", title.hjust = 0.5),
    shape = guide_legend(keywidth = 4, keyheight = 1.3, title.position = "top", title.hjust = 0.5)
  ) +
  theme(legend.position = "bottom")
ggsave(file.path(targetPath, str_c("VarianceRatio_Article.pdf")), height = ArticleHeight, width = ArticleWidth)
ggsave(file.path(targetPath, str_c("VarianceRatio_Presentation.pdf")), height = PresentationHeight, width = PresentationWidth)