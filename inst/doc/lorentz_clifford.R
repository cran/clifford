## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")
library("quadform")
library("lorentz")
library("jordan")

## ----label=hexstickerplotter, out.width='15%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))
knitr::include_graphics(system.file("help/figures/jordan.png", package = "jordan"))
knitr::include_graphics(system.file("help/figures/quadform.png", package = "quadform"))

## ----define4vec---------------------------------------------------------------
(fourvec <- c(1,5,3,2))  # a four-vector
u <- c(0.2,0.3,0.4)  # a three-velocity

## ----b1b2---------------------------------------------------------------------
(Bmat <- boost(u))  # Bmat = "B-matrix"

## ----matmultb-----------------------------------------------------------------
Bmat %*% fourvec

## ----cliff1-------------------------------------------------------------------
(scliff <- as.1vector(fourvec))

## ----squaredinterval----------------------------------------------------------
M <- diag(c(1,-1,-1,-1))
t(fourvec) %*% M %*% fourvec

## ----usequadform--------------------------------------------------------------
quad.form(M,fourvec)

## ----usescalprod--------------------------------------------------------------
signature(1,3)
scalprod(scliff,scliff)

## -----------------------------------------------------------------------------
phi <- 2.1234534   # just a made-up random value
B <- cosh(phi/2) + sinh(phi/2)*e(1:2) 
Binv <- rev(B) # cosh(phi/2)- sinh(phi/2)*e(1:2)
B*Binv

## -----------------------------------------------------------------------------
B <- function(phi){cosh(phi/2) + sinh(phi/2)*e(1:2)}
B(0.26) * B(1.33)
B(0.26 + 1.33) # should match

## -----------------------------------------------------------------------------
B3 <- function(phi,k){cosh(phi/2) + (
     +k[1]*sinh(phi/2)*e(c(1,2))
     +k[2]*sinh(phi/2)*e(c(1,3))
     +k[3]*sinh(phi/2)*e(c(1,4))
   )}
k <- function(kx,ky){c(kx, ky, sqrt(1-kx^2-ky^2))}
kx <- +0.23
ky <- -0.38


k1 <- k(kx=0.23, ky=-0.38)
sum(k1^2) # verify; should be = 1
zap(B3(0.3,k1)*B3(1.9,k1))  # zap() kills terms with small coefficients
zap(B3(0.3+1.9,k1)) # should match previous line (up to numerical accuracy)

## -----------------------------------------------------------------------------
k2 <- k(-0.5,0.1)
zap(B3(2.4,k1) * B3(1.9,k2))

## ----speedtocliff-------------------------------------------------------------
f <- function(u){
	phi <- acosh(gam(u))               # rapidity
	k <- cosines(u)                    # direction cosines
 	return(
	       cosh(phi/2)                 # t
	+ k[1]*sinh(phi/2)*basis(c(1,2))   # x
	+ k[2]*sinh(phi/2)*basis(c(1,3))   # y
	+ k[3]*sinh(phi/2)*basis(c(1,4))   # z
	)
}

## ----dasadx-------------------------------------------------------------------
u <- as.3vel(-c(0.2,0.3,0.4))  # negative (passive transform)
options(digits=5)
(B <- f(u))

## ----checkrev-----------------------------------------------------------------
B*rev(B)

## ----dsdasdfds----------------------------------------------------------------
zap(rev(B)*scliff*B)

## ----comparelorentz-----------------------------------------------------------
Bmat %*% fourvec

## ----checkinterval------------------------------------------------------------
jj <- rev(B)*scliff*B
scalprod(jj,jj)

## ----defv---------------------------------------------------------------------
u <- as.3vel(c(0.2, 0.3,  0.4))
v <- as.3vel(c(0.5, 0.0, -0.4))
w <- as.3vel(c(0.0, 0.7,  0.1))
Buvw <- f(u)*f(v)*f(w)
zap(Buvw)

## ----checkbuvw----------------------------------------------------------------
zap(Buvw*rev(Buvw))

## ----testiton1000-------------------------------------------------------------
n <- as.1vector(c(1,0,0,0))
zap(rev(Buvw) * n * Buvw)

## ----algecliff----------------------------------------------------------------------------------------------
signature(1,3)
L <- list(
    C     = basis(numeric()),
    e12   = basis(c(1,2)), e13 = basis(c(1,3)),
    e14   = basis(c(1,4)), e23 = basis(c(2,3)),
    e24   = basis(c(2,4)), e34 = basis(c(3,4)),
    e1234 = basis(1:4)
) 

out <- noquote(matrix("",8,8))
rownames(out) <- names(L)
colnames(out) <- names(L)
for(i in 1:8){
  for(j in 1:8){
    out[i,j] <- gsub('[_ ]','',as.character(L[[i]]*L[[j]]))
  }
}
options("width" = 110)
out

## ----echo=FALSE---------------------------------------------------------------------------------------------
signature(Inf) # restore default, to avoid interference with other vignettes

