## ----out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/onion.png", package = "onion"))

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")
library("onion")
library("disordR")

## ----echo=FALSE---------------------------------------------------------------
M <- noquote(matrix(c(
 " 1",   "  i", "  j",  "  k", "  ε", " εi", " εj", " εk",
 " i",  " -1",  "  k",  " -j", " εi", " -ε", " εk", "-εj",
 " j",  " -k",  " -1",  "  i", " εj", "-εk", " -ε", " εi",
 " k",   "  j", " -i",  " -1", " εk", " εj", "-εi", " -ε",
 " ε",  " εi",  " εj",  " εk", "  0", "  0", "  0", "  0",
"εi",  " -ε",   " εk",  "-εj", "  0", "  0", "  0", "  0",
"εj", "-εk",    " -ε",  " εi", "  0", "  0", "  0", "  0",
"εk",  " εj",   "-εi",  "  -ε","  0", "  0", "  0", "  0"
),8,8,byrow=TRUE))
rownames(M) <- c(" 1"," i"," j"," k"," ε","εi","εj","εk")
colnames(M) <- c(" 1","  i","  j","  k","  ε"," εi"," εj"," εk")

## ----echo=FALSE,print=TRUE----------------------------------------------------
M

## -----------------------------------------------------------------------------
DQ_example <- list(as.quaternion(c(5,8,-3,3),single=TRUE), as.quaternion(c(-1,2,1,12),single=TRUE))
DQ_example

## -----------------------------------------------------------------------------
DQ_prod_DQ <- function(DQ1,DQ2){
  A <- DQ1[[1]] ; B <- DQ1[[2]] ; C <- DQ2[[1]] ; D <- DQ2[[2]]
  list(A*C, A*D+B*C) 
}

## -----------------------------------------------------------------------------
signature(3,0)
e(4)*e(4)

## ----definecoercion-----------------------------------------------------------
`cliff_to_DQ` <- function(C){  # terms such as e_3 and e_123 and e_34 are silently discarded
    quat <- getcoeffs(C,list(numeric(0),c(1,2),c(1,3),c(2,3)))
    quat[-1] <- -quat[-1]

    epsi <- getcoeffs(C,list(4,c(1,2,4),c(1,3,4),c(2,3,4)))
    epsi[-1] <- -epsi[-1]

    return(list(as.quaternion(quat,single=TRUE),as.quaternion(epsi,single=TRUE)))
} 


`DQ_to_cliff` <- function(DQ){ # DQ is a two-element list of quaternions
  jj1 <- c(as.matrix(DQ[[1]]))
  jj1[-1] <- -jj1[-1]
  jj2 <- c(as.matrix(DQ[[2]]))
  jj2[-1] <- -jj2[-1]
  
 clifford(list(numeric(0),c(1,2),c(1,3),c(2,3),4,c(1,2,4),c(1,3,4),c(2,3,4)),c(jj1,jj2))
}

## -----------------------------------------------------------------------------
DQ1 <- list(as.quaternion(c(5,6,2,-7),single=TRUE),as.quaternion(c(-3,1,4,8),single=TRUE))
DQ2 <- list(as.quaternion(c(-1,3,1,4),single=TRUE),as.quaternion(c(1,9,-7,4),single=TRUE))
LHS <- DQ_to_cliff(DQ1) * DQ_to_cliff(DQ2)
RHS <- DQ_to_cliff(DQ_prod_DQ(DQ1,DQ2))
LHS == RHS

## -----------------------------------------------------------------------------
C1 <- clifford(list(numeric(0),c(1,3),c(1,2,4),c(1,3,4)),c(3,7,11,13))
C2 <- clifford(list(numeric(0),c(1,2),c(1,3),c(1,3,4)),c(2,5,6,17))

LHS <- cliff_to_DQ(C1*C2)
RHS <- DQ_prod_DQ(cliff_to_DQ(C1),cliff_to_DQ(C2))
identical(LHS,RHS)

## -----------------------------------------------------------------------------
identical(DQ1,cliff_to_DQ(DQ_to_cliff(DQ1)))
identical(C1,DQ_to_cliff(cliff_to_DQ(C1)))

