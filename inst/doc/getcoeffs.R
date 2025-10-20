## ----setup, include=FALSE-----------------------------------------------------
library("clifford")
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(echo = TRUE)
knit_print.function <- function(x, ...){dput(x)}
registerS3method(
  "knit_print", "function", knit_print.function,
  envir = asNamespace("knitr")
)

## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))

## ----showgetcoeffs------------------------------------------------------------
getcoeffs

## ----use_getcoeffs------------------------------------------------------------
set.seed(0)
(a <- rcliff())
getcoeffs(a,list(1:2, 0, c(2,5), c(1,5,6), c(2,6), 1:2))

## ----showcoeffs---------------------------------------------------------------
coeffs(a)

