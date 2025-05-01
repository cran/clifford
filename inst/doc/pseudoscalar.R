## ----setup, include=FALSE-----------------------------------------------------------------------------------
set.seed(0)
library("clifford")
library("permutations")
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(echo = TRUE)
knit_print.function <- function(x, ...){dput(x)}
registerS3method(
  "knit_print", "function", knit_print.function,
  envir = asNamespace("knitr")
)


## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE------------------------------
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/permutations.png", package = "permutations"))

## ----showI--------------------------------------------------------------------------------------------------
pseudoscalar

## ----error=TRUE---------------------------------------------------------------------------------------------
try({
pseudoscalar()
})

## ----settoseven---------------------------------------------------------------------------------------------
options(maxdim=7)
pseudoscalar()

## ----R3-----------------------------------------------------------------------------------------------------
options(maxdim=3)
signature(3)       # Cl(3,0)
(I <- pseudoscalar())
drop(I^2)

## ----mink---------------------------------------------------------------------------------------------------
options(maxdim=4)
signature(3,1)       # Cl(3,1)
I <- pseudoscalar()
drop(I^2)

## ----Isquaredplusone----------------------------------------------------------------------------------------
options(maxdim=4)
signature(2,2)       # Cl(2,2)
(I <- pseudoscalar())
drop(I^2)

## ----cl5----------------------------------------------------------------------------------------------------
options(maxdim=5)
signature(5)
I <- pseudoscalar()
ai <- list(); for(i in 1:5){ai[[i]] <- as.1vector(rnorm(5))}
ai[[1]] # the other 5 look very similar
Reduce(`^`,ai)

## ----permutevec---------------------------------------------------------------------------------------------
(p <- permutation("(12)(345)"))
is.even(p)

## ----permuteai----------------------------------------------------------------------------------------------
c(drop(Reduce(`^`,ai)),drop(Reduce(`^`,ai[as.word(p)])))

## ----showbits-----------------------------------------------------------------------------------------------
options(maxdim=7)   
signature(7)
(I <- pseudoscalar())
(a <- as.1vector(sample(1:10,5)))
(A <- rcliff())

## ----showverif----------------------------------------------------------------------------------------------
LHS <- cliffdotprod(a, A*I) # Usual idiom would be "a %.% (A*I)"
RHS <- (a^A)*I
LHS - RHS

## ----restoredefaults,echo=FALSE-----------------------------------------------------------------------------
options(maxdim=NULL) # restore defaults

