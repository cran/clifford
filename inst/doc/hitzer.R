## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----loadcliffordlibrary,echo=TRUE,print=FALSE,message=FALSE------------------
library(clifford)

## ----label=equation2.6--------------------------------------------------------
(M <- rcliff())
a1 <- cliffconj(M)
a2 <- gradeinv(rev(M))
a3 <- rev(gradeinv(M))
is.zero(a2-a1) & is.zero(a3-a1)

## ----label=equation6.2--------------------------------------------------------
(x <- rcliff(d=3,g=3))
x*cliffconj(x)

## ----label=defineequation6.3--------------------------------------------------
f <- function(x){
	jj <- x*cliffconj(x)
	is.real(jj*rev(jj))
}

## ----label=executeequation6.3-------------------------------------------------
signature(0,3)
f(rcliff(d=3,g=3))
signature(1,2)
f(rcliff(d=3,g=3))
signature(2,1)
f(rcliff(d=3,g=3))
signature(3,0)
f(rcliff(d=3,g=3))

## ----label=defineequation6.5--------------------------------------------------
RI3 <- function(x){ # right inverse
	jj <- cliffconj(x)*gradeinv(x)*rev(x)
	return(jj/drop(x*jj))
}


## -----------------------------------------------------------------------------
a <- 5+rcliff(d=3,g=3)
a
RI3(a)
zap(a*RI3(a))
zap(RI3(a)*a)

## ----label=defineequation7.7--------------------------------------------------
f77 <- function(x){
	jj <- cliffconj(x)*neg(x*cliffconj(x),3:4)
	return(jj/drop(x*jj))
}

f78 <- function(x){
	jj <- neg(cliffconj(x)*x,3:4)*cliffconj(x)
	return(jj/drop(jj*x))
}

a <- 3 + rcliff(d=4)
a
f77(a)
zap(a*f77(a))
zap(f77(a)*a)

## -----------------------------------------------------------------------------
set.seed(0)
sigs <- 0:4
left <- rep(NA,5)
right <- rep(NA,5)
diff <- rep(NA,5)
for(i in seq_along(sigs)){
	signature(sigs[i])
	a <- sample(1:9,1) + rcliff(d=4)
	left[i]  <- Mod(a*f77(a) -1)
	right[i] <- Mod(f77(a)*a -1)
	diff[i] <- Mod(f77(a)-f78(a))
}
left
right
diff

## ----label=section8-----------------------------------------------------------
f822 <- function(x){
	jj <- cliffconj(x)*gradeinv(x)*rev(x)
	jj <- jj*neg(x*jj,c(1L,4L))
	jj/drop(zap(x*jj))
}

## -----------------------------------------------------------------------------
a <- 7+clifford(list(1,3,5,1:2,c(1,5),c(3,4),1:3,2:4,c(2,3,5),1:4,2:5,c(1,2,3,5),1:5),1:13)
a
f822(a)
zap(a*f822(a))
zap(f822(a)*a)

## -----------------------------------------------------------------------------
sigs <- 0:6
diffl <- rep(NA,5)
diffr <- rep(NA,5)
for(i in seq_along(sigs)){
	signature(sigs[i])
	a <- sample(1:9,1) + rcliff(d=5)
	diffl[i] <- Mod(a*f822(a)-1)
	diffr[i] <- Mod(f822(a)*a-1)
}
diffl
diffr

## -----------------------------------------------------------------------------
f824 <- function(x){  # left inverse
	jj <- rev(x)*gradeinv(x)*cliffconj(x)
	jj <- neg(jj*x,c(1L,4L))*jj
	jj/drop(zap(jj*x))
}

## -----------------------------------------------------------------------------
zap(f824(x)*x)
zap(f822(x)*x)

## -----------------------------------------------------------------------------
signature(0,5)
Mod(f822(x) - f824(x))
signature(1,4)
Mod(f822(x) - f824(x))
signature(2,3)
Mod(f822(x) - f824(x))
signature(3,2)
Mod(f822(x) - f824(x))
signature(4,1)
Mod(f822(x) - f824(x))

## -----------------------------------------------------------------------------
a <- rcliff(d=7)   # Cl(4,3)
b <- rcliff(d=7)   # Cl(4,3)
signature(4,3)     # e1^2 = e2^2 = e3^2 = e4^2 = +1; e5^2 = ... = -1
ab <- a*b          # multiplication in Cl(4,3)

signature(0,7)   # e1^2 = ... = -1
cartan(a)*cartan(b) == cartan(ab) # multiplication in Cl(0,7)

## -----------------------------------------------------------------------------
cartan_inverse(cartan(a) * cartan(b)) == ab  # precalculated product!

## -----------------------------------------------------------------------------
signature(5,2); ab_sig5 <- a*b

signature(1,7)
cartan(a,2) * cartan(b,2)  == cartan(ab_sig5,2)
cartan_inverse(cartan(a,2) * cartan(b,2),2)  == ab_sig5

