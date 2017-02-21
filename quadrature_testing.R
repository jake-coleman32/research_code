##Quadrature rules

trapezoid <- function(fun, a, b, n=100) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b-a)/n
  x <- seq(a, b, by=h)
  y <- fun(x)
  s <- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
  return(s)
}

trap_me1 <- function(x,y){
  
}

trap_me <- function(x,y){
  n = length(y)
  idx = 2:(n-1)
  s <- c(y[1]*(x[2]-x[1])/2, y[idx]*(x[idx+1]-x[idx-1])/2, y[n]*(x[n]-x[n-1])/2)
  return(s)
}

log_dat <- log(mgCl)

mean_first <- trap_me(mgCl_x[[1]],apply(log_dat,2,mean))

trapz <- function (x, y){
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
}

##Simpson rules are terrible

simpson <- function(fun, a, b, n=100) {
  # numerical integral using Simpson's rule
  # assume a < b and n is an even positive integer
  h <- (b-a)/n
  x <- seq(a, b, by=h)
  if (n == 2) {
    s <- fun(x[1]) + 4*fun(x[2]) +fun(x[3])
  } else {
    s <- fun(x[1]) + fun(x[n+1]) + 2*sum(fun(x[seq(2,n,by=2)])) + 4 *sum(fun(x[seq(3,n-1, by=2)]))
  }
  s <- s*h/3
  return(s)
}

simpson_v3 <- function(x,y){
  y[is.nan(y)]=0
  s <- c(y[1], y[n+1], 2*y[seq(2,n,by=2)], 
    4 *y[seq(3,n-1, by=2)])
  s <- s*h/3
  return(sum(s))
}

phi <- function(x) exp(-x^2/2)/sqrt(2*pi)  ##normal distribution
simpson_v2(phi, -Inf, 3, n=100)
pnorm(3)
