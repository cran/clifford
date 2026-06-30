## ----setup, include=FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")
library("onion")
library("jordan")
library("emulator")

## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE------------------------------
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/jordan.png", package = "jordan"))

## -----------------------------------------------------------------------------------------------------------
s0 <- matrix(c(1,0,0,1),2,2)
sx <- matrix(c(0,1,1,0),2,2)
sy <- matrix(c(0,1i,-1i,0),2,2)
sz <- matrix(c(1,0,0,-1),2,2)

## -----------------------------------------------------------------------------------------------------------
matrix_to_clifford <- function(M){
      (Re(M[1,1] + M[2,2]))/2             + 
      (Re(M[1,1] - M[2,2]))/2*e(c(  3  )) +
      (Im(M[1,1] + M[2,2]))/2*e(c(1,2,3)) + 
      (Im(M[1,1] - M[2,2]))/2*e(c(1,2  )) +

      (Re(M[2,1] + M[1,2]))/2*e(c(1    )) + 
      (Re(M[2,1] - M[1,2]))/2*e(c(1,  3)) +
      (Im(M[2,1] + M[1,2]))/2*e(c(  2,3)) + 
      (Im(M[2,1] - M[1,2]))/2*e(c(  2  ))
}

## -----------------------------------------------------------------------------------------------------------
rmat <- function(...){matrix(rnorm(4),2,2) + 1i*matrix(rnorm(4),2,2)}
M <- rmat()
M
matrix_to_clifford(M)

## -----------------------------------------------------------------------------------------------------------
M1 <- rmat()
M2 <- rmat()

diff <- matrix_to_clifford(M1)*matrix_to_clifford(M2) - matrix_to_clifford(M1 %*% M2)
diff
Mod(diff)

## -----------------------------------------------------------------------------------------------------------
`clifford_to_matrix` <- function(C){
   return(
                          const(C)*s0 + getcoeffs(C,list(1))*sx 
  +           getcoeffs(C,list(2))*sy + getcoeffs(C,list(3))*sz
  + getcoeffs(C,list(c(1,2,3)))*1i*s0 + getcoeffs(C,list(c(  2,3)))*1i*sx 
  - getcoeffs(C,list(c(1,  3)))*1i*sy + getcoeffs(C,list(c(1,2  )))*1i*sz
  )
} 

## -----------------------------------------------------------------------------------------------------------
rc <- function(...){rcliff(100,d=3,g=3)}
C <- 104 + rc()
C
clifford_to_matrix(C)

## -----------------------------------------------------------------------------------------------------------
clifford_to_matrix(matrix_to_clifford(M)) - M 
matrix_to_clifford(clifford_to_matrix(C))- C

## -----------------------------------------------------------------------------------------------------------
C1 <- 222 + rc()
C2 <- 333 + rc()
clifford_to_matrix(C1*C2) - clifford_to_matrix(C1)%*%clifford_to_matrix(C2)

## -----------------------------------------------------------------------------------------------------------
M1 <- as.1matrix(rchm(1,2))
M2 <- as.1matrix(rchm(1,2))
M1
M2
p1 <- (M1 %*% M2 + M2 %*% M1)/2
p1 - ht(p1)  # zero for Hermitian matrices

## -----------------------------------------------------------------------------------------------------------
C1 <- matrix_to_clifford(M1)
C2 <- matrix_to_clifford(M2)
p2 <- (C1 * C2 + C2 * C1)/2
p2

## -----------------------------------------------------------------------------------------------------------
grades(p2)

