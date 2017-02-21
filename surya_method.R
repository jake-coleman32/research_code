#Trying to replicate Surya's method
#Algorithm 1
library(dplyr)
library(VGAM)
library(matrixcalc)
library(truncnorm)
library(caTools)



cov_mat <- function(x,x_pred,ell,lambda, nugget=0.){
  inds <- 1:length(x)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(x,x_pred)))[inds,-inds]^2/ell^2)
  
  if(length(x)==length(x_pred)){
    return(out_mat + nugget*diag(length(x)))
  }
  
  return(out_mat)
}


sqrt_inv <- function(mat){
  #Take inverse first
  inv <- chol2inv(chol(mat))
  
  #Take sqrt next
  e <- eigen(inv,symmetric=TRUE)
  out_mat <- e$vectors%*%diag(sqrt(e$values))%*%t(e$vectors)
  return(out_mat)
}

#Not using currently
gp_pred <- function(sqrt_inv_mat,knots,x,y,ell,lambda,nugget=0){
  #Need sqrt_inv function
  A_gam <- sqrt_inv_mat%*%cov_mat(x=knots,x_pred=y,ell=ell,lambda=lambda)
  
  #may not be t()
  return(t(matrix(x))%*%A_gam)
}


density_pred <- function(sqrt_inv_mat,knots,x,y,ell,lambda,nugget=0){
  int_grid <- seq(0,1,length=101)
  
  #Knot preds
  A_gam_knot = sqrt_inv_mat%*%cov_mat(x=knots,x_pred=y,ell=ell,lambda=lambda)
  knot_vals = exp(t(matrix(x))%*%A_gam_knot)
  
  #All preds
  A_gam_all =  sqrt_inv_mat%*%cov_mat(x=knots,x_pred=int_grid,ell=ell,
                                      lambda=lambda)
  all_vals <- exp(t(matrix(x))%*%A_gam_all)

  est_int <- trapz(int_grid,all_vals)
  
  return(list(dens_est = knot_vals/est_int, f_est = all_vals/est_int))
}

#test <- density_pred(sig_gam_sqrt_inv_cur,knots,x_cur,y_dat,ell_cur,lambda_cur)


alpha = 7
beta = 3

m <- 21
n <- 100
knots <- seq(0,1,length=m)

#Data comes from Beta(alpha,beta) distribution
y_dat <- rbeta(n,alpha,beta)
hist(y_dat)
plot(grid,dbeta(grid,alpha,beta))

#Initial values for X
init <- mvrnorm(1,rep(0,m),diag(m))

#Tuning parameter for MH moves
kap_x <- rep(1,m)
kap_lam <-1
kap_ell <- 0.5

delta_ell <- 0.3
delta_lam <- 0.3
  
#Setup output values
N = 1E4
burnin <- 1E3

f_out <- matrix(0,N,101)
x_out <- matrix(0,N,m)
ell_out <- lambda_out <- numeric(N)

lam_ratio <- ell_ratio <- 0
x_ratios <- numeric(m)


#Draw from covariance parameters
ell_loc <- 0
ell_scale <- 1
ell_shape <- 0.5


lam_a <- 3
lam_b <- 3

ell_a <- 1
ell_b <- 1

#draws from priors
(lambda_cur <- rgamma(1,lam_a,lam_b))
#(ell_cur <- rfrechet(1,ell_loc,ell_scale,ell_shape))
ell_cur <- rgamma(1,ell_a,ell_b)

sig_gam_sqrt_inv_cur <- sqrt_inv(cov_mat(x=knots,x_pred=knots,ell=ell_cur,
                                         lambda=lambda_cur,
                                         nugget=1E-5))

x_cur <- init
vals_cur <- density_pred(sig_gam_sqrt_inv_cur,knots,x_cur,y_dat,
                         ell_cur,lambda_cur)
log_f_cur <- log(vals_cur$dens_est)

start_t <- proc.time()
for(i in 1:(N+burnin)){
  #update X
  
  for(j in 1:m){
    xj_st <- rnorm(1,x_cur[j],kap_x)
    x_st <- x_cur
    x_st[j] <- xj_st
    vals_st <- density_pred(sig_gam_sqrt_inv_cur,knots, x_st,y_dat,
                            ell_cur,lambda_cur)
    log_f_st <- log(vals_st$dens_est)
    
    xj_cur= x_cur[j]
    #For each x, MH step
    #If accept x_j, update all current vals
    if((dnorm(xj_st,log=TRUE) + sum(log_f_st) - dnorm(xj_cur,log=TRUE) - 
        sum(log_f_cur)) > log(runif(1))){#Er, check this
      x_cur[j] <- xj_st
      vals_cur <- vals_st
      log_f_cur <- log_f_st
      x_ratios[j] <- x_ratios[j] + 1
    }
  }
  #############
  #Update ell
  #############
  #propose new ell - truncated normal
#  ell_st <- rtruncnorm(1,a=0,mean = ell_cur,sd=kap_ell)
  z_ell <- rnorm(1)
  ell_st <- ell_cur*exp(delta_ell*z_ell)
  #adjust_ratio <- pnorm(ell_cur,log=TRUE) - pnorm(ell_st,log=TRUE)
  adjust_ratio <- delta_ell*z_ell
  
  #get likelihood from ell_st
  sig_gam_sqrt_inv_st <- sqrt_inv(cov_mat(x=knots,x_pred=knots,ell=ell_st,
                                          lambda=lambda_cur,nugget=1E-5))
  vals_st <- density_pred(sig_gam_sqrt_inv_st,knots,x_cur,y_dat, ell = ell_st,
                          lambda = lambda_cur)
  log_f_st <- log(vals_st$dens_est)
  
  #Do accept/reject step
#   if((dfrechet(ell_st,ell_loc,ell_scale,ell_shape, log=TRUE) + sum(log_f_st) -
#       dfrechet(ell_cur,ell_loc,ell_scale,ell_shape, log=TRUE) - sum(log_f_cur) +
#        adjust_ratio + rexp(1)) > 0){#Er, check this
  if((dgamma(ell_st,ell_a,ell_b, log=TRUE) + sum(log_f_st) -
      dgamma(ell_cur,ell_a,ell_b, log=TRUE) - sum(log_f_cur) +
      adjust_ratio + rexp(1)) > 0){#Er, check this
    
    ell_cur <- ell_st
    vals_cur <- vals_st
    log_f_cur <- log_f_st
    ell_ratio <- ell_ratio + 1
  }
  
  ############
  #Update lambda
  #############
  
  #propose new lambda - truncated normal
  #lambda_st <- rtruncnorm(1,a=0,mean = lambda_cur,sd=kap_lam)
  z_lam <- rnorm(1)
  lambda_st <- lambda_cur*exp(delta_lam*z_lam)
  #adjust_ratio <- pnorm(lambda_cur,log=TRUE) - pnorm(lambda_st,log=TRUE)
  adjust_raito <- delta_lam*z_lam
  
  #get likelihood from lambda_st
  sig_gam_sqrt_inv_st <- sqrt_inv(cov_mat(x=knots,x_pred = knots, ell=ell_cur,
                                          lambda=lambda_st, nugget = 1E-5))
  vals_st <- density_pred(sig_gam_sqrt_inv_st,knots,x_cur,y_dat, ell = ell_cur,
                          lambda = lambda_st)
  log_f_st <- log(vals_st$dens_est)
  
  #Do accept/reject step
  if((dgamma(lambda_st,lam_a,lam_b, log=TRUE) + sum(log_f_st) -
      dgamma(lambda_cur,lam_a,lam_b, log=TRUE) - sum(log_f_cur) +
      adjust_ratio + rexp(1)) > 0){#Er, check this
    lambda_cur <- lambda_st
    vals_cur <- vals_st
    log_f_cur <- log_f_st
    lam_ratio <- lam_ratio + 1
  }
  
  sig_gam_sqrt_inv_cur <- sqrt_inv(cov_mat(x=knots,x_pred=knots,ell=ell_cur,
                                           lambda=lambda_cur, nugget=1E-5))
  
  
  #Auto-tune?
  #if past burnin, keep current things
  if(i>burnin){
    x_out[i-burnin,] <- x_cur
    f_out[i-burnin,] <- vals_cur$f_est
    lambda_out[i-burnin] <- lambda_cur
    ell_out[i-burnin] <- ell_cur
  }

}
proc.time() - start_t

(lam_ratio <- lam_ratio/(N+burnin))
(ell_ratio <- ell_ratio/(N+burnin))
(x_ratios <- x_ratios/(N+burnin))

grid <- seq(0,1,length=101)
f_mean <- apply(f_out,2,mean)
plot(grid,f_mean)
plot(lambda_out,type = 'l')
plot(ell_out,type = 'l')
plot(f_out[,1],type = 'l')
plot(x_out[,3],type ='l')

plot(grid,f_mean,main = 'Estimate for Beta(7,3)', type = 'l',
     xlab = 'y',ylab = 'Density',ylim = c(0,4),lwd=2)
lines(grid,apply(f_out,2,quantile,0.025),col = 'blue',lty=2)
lines(grid,apply(f_out,2,quantile,0.975),col = 'blue',lty=2)
lines(grid,dbeta(grid,alpha,beta),col = 'red',lwd=2)
lines(density(y_dat),col = 'green4',lwd = 2)
legend('topleft',c('Posterior Mean','Credible Int','Kernel Density','Truth'),
      col = c('black','blue','green4','red'),lwd = c(2,1,2,2),lty = c(1,2,1,1))
