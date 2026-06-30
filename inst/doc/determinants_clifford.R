## ----out.width='20%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))

## ----formalshow---------------------------------------------------------------
suppressMessages(library("clifford"))
set.seed(0)
(M <- matrix(rnorm(9),3,3))
o <- as.1vector(M[,1]) ^ as.1vector(M[,2]) ^ as.1vector(M[,3])
Adag <- rev(e(seq_len(3)))
c(drop(Adag %.% o), det(M))

## ----directwedge--------------------------------------------------------------
as.1vector(M[,1]) ^ as.1vector(M[,2]) ^ as.1vector(M[,3])

## ----usecoeffs----------------------------------------------------------------
coeffs(as.1vector(M[,1]) ^ as.1vector(M[,2]) ^ as.1vector(M[,3]))

## ----detofi3------------------------------------------------------------------
coeffs(e(1) ^ e(2) ^ e(3))

## ----largermatrix-------------------------------------------------------------
cliff_det <- function(M){
  o <- as.1vector(M[,1])
  for(i in 2:nrow(M)){
    o <- o ^ as.1vector(M[,i])
  }
  return(coeffs(o))
}

## ----usenewfunction-----------------------------------------------------------
M <- matrix(rnorm(100),10,10)
LHS <- det(M)
RHS <- cliff_det(M)
c(LHS,RHS,LHS-RHS)

