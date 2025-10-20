## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("clifford")

## ----out.width='20%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))

## ----setdims------------------------------------------------------------------
dimension <- 3
options("maxdim" = dimension+2)  # paranoid safety measure
signature(dimension + 1,1)
eplus <- basis(dimension+1)
eminus <- basis(dimension + 2)

e0 <-  (eminus - eplus)/2
einf <- eminus + eplus
E <- e0 ^ einf

## ----showe0einfE--------------------------------------------------------------
e0
einf
E

## ----definepointfunc----------------------------------------------------------
point <- function(x){ as.1vector(x) + sum(x^2)*einf/2 + e0 }

## ----pointexample-------------------------------------------------------------
a <- c(1,2,5)
b <- c(2,2,2)
point(a)
point(b)

## ----verifypointfunc----------------------------------------------------------
c(conformal=drop(point(a) %.% point(b)), Euclidean = -sum((a-b)^2)/2)

## ----definesphere-------------------------------------------------------------
sphere <- function(x,r){ point(x) - r^2*einf/2}
S <- sphere(1:3,5)  # center (1,2,3) radius 5:
S

## ----verifyradius-------------------------------------------------------------
drop(S^2)   # 5^2 = 25

## ----calcsandichprod----------------------------------------------------------
S*einf*S

## ----usescaling---------------------------------------------------------------
-S*einf*S/2

## ----recog123-----------------------------------------------------------------
point(1:3)

## ----fourpoints---------------------------------------------------------------
origin <- point(c(0,0,0))
px <- point(c(1,0,0))
py <- point(c(0,1,0))
pz <- point(c(0,0,1))

(S <- origin ^ px ^ py ^ pz)

## ----definespherestar---------------------------------------------------------
spherestar <- function(...){Reduce(`^`,list(...))}
spherestar(origin, px, py, pz)

## ----checksphere--------------------------------------------------------------
p <- point(c(1,1,1+sqrt(3))/2)
Mod(p ^ S)

## ----defineplanefunc----------------------------------------------------------
plane <- function(n,d){ as.1vector(n/sqrt(sum(n^2))) + d*einf}

## ----generatend---------------------------------------------------------------
n <- c(1,2,5)
nhat <- n/sqrt(sum(n^2))
d <- 7
Pi <- plane(nhat,d)
Pi

## ----generateuvp1p2p3---------------------------------------------------------
u <- c(2,-1,0)
v <- c(5,0,-1)
P1 <- point(d*nhat)                
P2 <- point(d*nhat + 1.3*u + 3.44*v)
P3 <- point(d*nhat - 6.1*u + 1.02*v)

## ----verifyplaneipns----------------------------------------------------------
c(drop(Pi %.% P1),drop(Pi %.% P2),drop(Pi %.% P3))

## ----plane3points-------------------------------------------------------------
Pi2 <- P1 ^ P2 ^ P3 ^ einf

## ----checkonplane-------------------------------------------------------------
p4 <- point(d*nhat + 7.6*u - 9.23*v)
Mod(Pi2 ^ p4)

## ----definecirclefun----------------------------------------------------------
circle <- function(S1,S2){  # IPNS
    S1 ^ S2
}

## ----usecirclefun-------------------------------------------------------------
circle(sphere(1:3,5),sphere(c(1.1,2.1,3.4),6))

## ----definecirclestar---------------------------------------------------------
circlestar <- function(...){  # OPNS; A^B^C
	jj <- list(...)
    stopifnot(length(jj) == dimension)
    Reduce(`^`,lapply(jj,point))
}

(CIRC <- circlestar(c(1,2,3),c(5,6,3),c(8,8,-2)))

## ----usenumericstofindapointoncircle,echo=FALSE,cache=TRUE--------------------
M <- rbind (c(1,2,3),c(5,6,3),c(8,8,-2))
badness_center <- function(pos){
  bad1 <-  sd(c(
  sum((pos-M[1,])^2),
  sum((pos-M[2,])^2),
  sum((pos-M[3,])^2)))   # pos should be equidistant from rows of M
  bad2 <- abs(det(sweep(M,2,pos))) # pos should be collinear with rows of M
  return(bad1+bad2)
  }

center <- nlm(badness_center,c(0,0,0))$estimate
badness_pointoncircle <- function(pos){ # pos is a point [that should be] on the circle
bad1 <-  var(c(sum((center-M[1,])^2), sum((center-M[2,])^2), sum((center-M[3,])^2),   sum((center-pos  )^2)))
bad2 <-  abs(det(sweep(M,2,pos))) # pos should be coplanar with M[i,]
return(bad1+bad2)
}

poc <- point(nlm(badness_pointoncircle,c(0,0,0))$estimate)

## ----verifycircleusingnum-----------------------------------------------------
poc  # point on circle, found numerically [chunk omitted]
poc ^ CIRC

## ----defline------------------------------------------------------------------
line <- function(P1,P2){ P1 ^ P2 }

## ----defpp--------------------------------------------------------------------
pointpair <- function(S1,S2,S3){ S1 ^ S2 ^ S3 }

## ----verifypointpairassoc-----------------------------------------------------
S1 <- sphere(c(3,2,4),3)
S2 <- sphere(c(3,1,4),4)
S3 <- sphere(c(1,3,3),3)
(S1^S2)^S3
(S1^S2)^S3 == S1^(S2^S3)

## ----restoredefault, echo=FALSE-----------------------------------------------
options("maxdim" = NULL) # restore default

