## ----setup, include=FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")
library("onion")

## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE------------------------------
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/onion.png", package = "onion"))

## -----------------------------------------------------------------------------------------------------------
options(maxdim=3)  # paranoid safety measure

## -----------------------------------------------------------------------------------------------------------
q1 <- +1 + 2* -e(c(1,2)) + 3*-e(c(1,3)) + 4*-e(c(2,3))
q1
q2 <- -2 + 1* -e(c(1,2)) - 2*-e(c(1,3)) + 1*-e(c(2,3))
q2
q1*q2

## ----clifftoquat--------------------------------------------------------------------------------------------
`clifford_to_quaternion` <- function(C){
    C <- as.clifford(C)
    tC <- disordR::elements(terms(C))
    stopifnot(all(c(tC,recursive=TRUE) <= 3))
    jj <- unlist(lapply(tC,length))
    stopifnot(all(jj <= 2))    # safety check
    stopifnot(all(jj%%2 == 0)) # safety check
    out <- matrix(c(const(C), -getcoeffs(C,list(c(1,2),c(1,3),c(2,3)) )))
    as.quaternion(out)
}

## ----defineqtoc---------------------------------------------------------------------------------------------
`quaternion_to_clifford` <- function(Q){
  Q <- as.numeric(Q)
  stopifnot(length(Q)==4)
  clifford(list(numeric(0),c(1,2),c(1,3),c(2,3)),c(Q[1],-Q[2:4]))
}

## ----defqh--------------------------------------------------------------------------------------------------
q1 <- +1 + 2* -e(c(1,2)) + 3*-e(c(1,3)) + 4*-e(c(2,3))
q2 <- -2 + 1* -e(c(1,2)) - 2*-e(c(1,3)) + 1*-e(c(2,3))
H1 <- as.quaternion(c(3,-5,2,1),single=TRUE)
H2 <- as.quaternion(c(1,2,-2,2),single=TRUE)

## ----verifyinverse------------------------------------------------------------------------------------------
c(  # check they are inverses of one another
q1 == quaternion_to_clifford(clifford_to_quaternion(q1)),
q2 == quaternion_to_clifford(clifford_to_quaternion(q2)),
H1 == clifford_to_quaternion(quaternion_to_clifford(H1)),
H2 == clifford_to_quaternion(quaternion_to_clifford(H2))
)

## -----------------------------------------------------------------------------------------------------------
c(
q1*q2 == quaternion_to_clifford(clifford_to_quaternion(q1)*clifford_to_quaternion(q2)),
H1*H2 == clifford_to_quaternion(quaternion_to_clifford(H1)*quaternion_to_clifford(H2))
)

## -----------------------------------------------------------------------------------------------------------
signature(0,3)
cliff2quat <- function(C){
  out <- getcoeffs(C,list(numeric(0),c(2,3),c(1,3),c(1,2)))
  out[2] <- -out[2]
  as.quaternion(out,single=TRUE)
}

quat2cliff <- function(H){
  jj <- c(as.matrix(H))
  jj[2] <- -jj[2]
  clifford(list(numeric(0),c(2,3),c(1,3),c(1,2)),jj)
}

## -----------------------------------------------------------------------------------------------------------
c(
  cliff2quat(quat2cliff(H1)) == H1,
  cliff2quat(quat2cliff(H2)) == H2,
  quat2cliff(cliff2quat(q1)) == q1,
  quat2cliff(cliff2quat(q2)) == q2,
  cliff2quat(q1*q2) == cliff2quat(q1) * cliff2quat(q2),
  quat2cliff(H1*H2) == quat2cliff(H1) * quat2cliff(H2)
)

## ----echo=FALSE---------------------------------------------------------------------------------------------
# restore default
options(maxdim = NULL)
signature(Inf)

