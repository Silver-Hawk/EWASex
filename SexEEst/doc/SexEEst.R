## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results='asis', fig.cap='prediction output'-----------------
precompiled_preds <- read.csv("allPredictionsForVignette.csv")
knitr::kable(head(precompiled_preds, 10))

## ---- echo=FALSE, results='asis', fig.cap='prediction output'-----------------
knitr::kable(table(precompiled_preds[,c('predictedGender', 'Sex')]), row.names = T)

## ----echo = FALSE, message=FALSE, fig.align='center', fig.cap='Prediction plot of the genders'----
knitr::include_graphics("PredictionPlot1024_1.png")

