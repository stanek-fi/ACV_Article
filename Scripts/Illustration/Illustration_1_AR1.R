library(ACV)
library(data.table)
library(ggplot2)
library(stringr)
rm(list = ls())
targetPath <- "Outputs"
suppressWarnings(dir.create(targetPath))
ArticleHeight <- 6
ArticleWidth <- 6
PresentationHeight <- 12
PresentationWidth <- 6

set.seed(1)
Tobs <- 20
order <- c(1, 0, 0)
y <- arima.sim(list(ar = c(0.2)), Tobs, sd = 1)
m <- 16

algorithm <- function(yInSample, yOutSample, h) {
  model <- arima(yInSample, order = order)
  modelForecast <- arima(c(yInSample, yOutSample), order, fixed = model$coef)
  return(list(
    yhatInSample = fitted(model),
    yhatOutSample = fitted(modelForecast)[-seq_along(yInSample)]
  ))
}

LambdaList <- list()

v <- 1
h <- v
Phi <- tsACV(y = y, algorithm = algorithm, m = m, h = h, v = v)
Lambda <- Phi
Lambda[!is.na(Phi)] <- estimateL(Phi = Phi, method = "regular")$lambda
LambdaList[["v1_CV"]] <- Lambda
Lambda <- Phi
Lambda[!is.na(Phi)] <- estimateL(Phi = Phi, method = "augmented")$lambda
LambdaList[["v1_ACV"]] <- Lambda

v <- Tobs - m
h <- v
Phi <- tsACV(y = y, algorithm = algorithm, m = m, h = h, v = v)
Lambda <- Phi
Lambda[!is.na(Phi)] <- estimateL(Phi = Phi, method = "regular")$lambda
LambdaList[["vn_CV"]] <- Lambda
Lambda <- Phi
Lambda[!is.na(Phi)] <- estimateL(Phi = Phi, method = "augmented")$lambda
LambdaList[["vn_ACV"]] <- Lambda

for (i in seq_along(LambdaList)) {
  Lambda <- as.data.table(LambdaList[[i]])
  Lambda[, t := 1:.N]
  mLambda <- melt(Lambda, id.vars = "t", variable.name = c("i"))
  mLambda[, t := factor(t, levels = sort(unique(t), decreasing = T))]
  mLambda <- na.omit(mLambda)
  ggplot(data = mLambda, aes(x = i, y = t, fill = value)) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = format(round(value, 3), nsmall = 3))) +
    # scale_fill_gradient2(midpoint=0, low="red",high="blue",limits=c(-0.21,0.21))+
    scale_fill_gradient2(midpoint = 0, low = "red", high = "blue", limits = c(-0.26, 0.26)) +
    labs(fill = expression(lambda)) +
    scale_x_discrete(position = "top") +
    guides(fill = guide_colourbar(barwidth = 15, barheight = 1, title.position = "top", title.hjust = 0.5)) +
    theme(legend.position = "bottom") +
    coord_fixed(ratio = 0.25) +
    xlab("") +
    ylab("")
  ggsave(file.path(targetPath, str_c("Illustration_", names(LambdaList)[i], "_Article.pdf")), height = ArticleHeight, width = ArticleWidth)
  ggsave(file.path(targetPath, str_c("Illustration_", names(LambdaList)[i], "_Presentation.pdf")), height = PresentationHeight, width = PresentationWidth)
}
