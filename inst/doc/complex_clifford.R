## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")

## ----out.width='20%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))

## ----label=coercionfuncs------------------------------------------------------
signature(0,1)
options(maxdim=1) # paranoid-level safety measure
complex_to_clifford <- function(z){Re(z) + e(1)*Im(z)}
clifford_to_complex <- function(C){const(C) + 1i*getcoeffs(C,1)}
clifford_to_complex <- function(C){const(C) + 1i*coeffs(Im(C))}


## ----label=numversetup--------------------------------------------------------
z1 <- 35 + 67i
z2 <- -2 + 12i

## ----label=examplecoerce------------------------------------------------------
z1
complex_to_clifford(z1)

## ----label=numverdo-----------------------------------------------------------
complex_to_clifford(z1) * complex_to_clifford(z2) == complex_to_clifford(z1*z2)

## ----label=otherway-----------------------------------------------------------
(C1 <- 23 + 7*e(1))
clifford_to_complex(C1)
C2 <- 2  - 8*e(1)
clifford_to_complex(C1)*clifford_to_complex(C2) == clifford_to_complex(C1*C2)

## ----label=secondcoercions----------------------------------------------------
options(maxdim=2)  # paranoid-level safety measure
signature(2)
complex_to_clifford <- function(z){Re(z) + e(1:2)*Im(z)}
clifford_to_complex <- function(C){const(C) + 1i*coeffs(Im(C))}

## ----label=thennumver---------------------------------------------------------
z1 <- 35 + 67i
z2 <- -2 + 12i
complex_to_clifford(z1) * complex_to_clifford(z2) == complex_to_clifford(z1*z2)

## ----label=gobackwards--------------------------------------------------------
C1 <- 23 + 7*e(1:2)
C2 <- 2  - 8*e(1:2)
clifford_to_complex(C1)*clifford_to_complex(C2) == clifford_to_complex(C1*C2)

## ----label=useothersig--------------------------------------------------------
signature(0,2)
c(
complex_to_clifford(z1)*complex_to_clifford(z2) == complex_to_clifford(z1*z2),
clifford_to_complex(C1)*clifford_to_complex(C2) == clifford_to_complex(C1*C2)
)

## ----label=restoresig---------------------------------------------------------
options(maxdim=NULL)
signature(Inf)

