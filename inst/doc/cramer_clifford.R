## ----out.width='20%', out.extra='style="float:right; padding:10px"',echo=FALSE----
knitr::include_graphics(system.file("help/figures/clifford.png", package = "clifford"))

## ----label=loadlib,results="hide",include=FALSE-------------------------------
library("clifford",quietly=TRUE)  # document requires package version 1.0-9 or above
set.seed(0)

## ----label=useone,eval=TRUE---------------------------------------------------
a <- as.1vector(runif(3))
b <- as.1vector(runif(3))
c <- as.1vector(runif(3))

(x <- as.1vector(1:3))

options(maxdim = 3)  # needed to drop() pseudoscalars
abc <- drop(a ^ b ^ c)

alpha <- drop(x ^ b ^ c)/abc
beta  <- drop(a ^ x ^ c)/abc
gamma <- drop(a ^ b ^ x)/abc

c(alpha,beta,gamma)

alpha*a + beta*b + gamma*c
Mod(alpha*a + beta*b + gamma*c-x)

## ----label=showdrop-----------------------------------------------------------
y <- a*1 + b*2 + c*3
c(
    drop(y ^ b ^ c)/abc,
    drop(a ^ y ^ c)/abc,
    drop(a ^ b ^ y)/abc
)

## ----arbdim-------------------------------------------------------------------
n <- 5                                                 # dimensionality of space
options(maxdim=5)                                      # safety precaution
x <- as.1vector(seq_len(n))                            # target vector
x
L <- replicate(n,as.1vector(rnorm(n)),simplify=FALSE)  # spanning vectors
subst <- function(L,n,x){L[[n]] <- x; return(L)}       # list substitution
coeff <- function(n,L,x){
      drop(Reduce(`^`,subst(L,n,x))/Reduce(`^`,L))
}

## ----showcoeffs---------------------------------------------------------------
(alpha <- sapply(seq_len(n),coeff,L,x))

## ----reconvec-----------------------------------------------------------------
out <- as.clifford(0)
f <- function(i){alpha[i]*L[[i]]}
for(i in seq_len(n)){
      out <- out + f(i)
}
Mod(out-x)  # zero to numerical precision

## ----somewhatslicker----------------------------------------------------------
Reduce(`+`,sapply(seq_len(n),f,simplify=FALSE))

## ----conversely---------------------------------------------------------------
coeffs <- 15:11
x <- 0
for(i in seq_len(5)){x <- x + coeffs[i]*L[[i]]}
x

## ----findcoeffs---------------------------------------------------------------
sapply(seq_len(n),coeff,L,x)

## ----restoredefaults, echo=FALSE----------------------------------------------
options(maxdim=NULL) # restore default for maxdim: options persist
                     # between vignettes, and leaving maxdim set
                     # causes problems for other vignettes

