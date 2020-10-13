## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results='asis', fig.cap='prediction output'-----------------
precompiled_preds <- read.csv("allPredictionsForVignette.csv")
knitr::kable(head(precompiled_preds, 10))

