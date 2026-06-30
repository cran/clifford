## ----setup, include=FALSE-----------------------------------------------------------------------------------
set.seed(0)
library("clifford")
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(echo = TRUE)
knit_print.function <- function(x, ...){dput(x)}
registerS3method(
  "knit_print", "function", knit_print.function,
  envir = asNamespace("knitr")
)

## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE------------------------------
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/lorentz.png", package = "lorentz"))

## ----showsigdef---------------------------------------------------------------------------------------------
signature

## ----label=sig34--------------------------------------------------------------------------------------------
signature(1,2)

## ----showp--------------------------------------------------------------------------------------------------
c(drop(e(1)^2),drop(e(2)^2),drop(e(3)^2))

## ----showg4-------------------------------------------------------------------------------------------------
c(drop(e(4)^2),drop(e(100)^2))

## ----usemaxdim----------------------------------------------------------------------------------------------
options(maxdim = 4)

## ----showerror, error=TRUE----------------------------------------------------------------------------------
try({
c(drop(e(1)^2),drop(e(2)^2),drop(e(3)^2),drop(e(4)^2))
e(5)
})

## ----showsig------------------------------------------------------------------------------------------------
signature()

## ----showinfinitesig----------------------------------------------------------------------------------------
options(maxdim=NULL)
signature(Inf)
signature()

## ----dputinfinitesig----------------------------------------------------------------------------------------
dput(signature())

