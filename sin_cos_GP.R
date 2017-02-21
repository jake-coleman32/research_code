
r_pow <- function(r,N){
  if(abs(r)>=1){
    stop("r parameter too big")
  }
  return(r^c(0:N))
}

cos_mult <- function(t,t_p,n){
  return(outer(cos(2*pi*n*t),cos(2*pi*n*t_p)))
}
sin_mult <- function(t,t_p,n){
  return(outer(sin(2*pi*n*t),sin(2*pi*n*t_p)))
}


#This is the same as the super naive, for-loop method
get_sin_cov_open <- function(t,t_prime,r,N,lam = 1,nugget = 0.){
  a = b = r_pow(r,N)
  #a_n = b_n = 1/c(1:N)^2
  
  lambda = lam*(1-r^2)/(2*r^2)
  
  out_mat <- matrix(0,length(t),length(t_prime))
  for(n in 1:N){
    out_mat <- out_mat + lambda*2*a[n]^2*cos_mult(t,t_prime,n) +
      lambda*2*b[n]^2*sin_mult(t,t_prime,n)
  }  
  
  
  if(length(t)==length(t_prime)){
    out_mat = out_mat + nugget*diag(length(t))
  }
  return(out_mat)
}

system.time(test2 <- get_sin_cov_open(x1,x1,.8,50,lam=2) )

##Closed form for a_n = b_n = r^n
crazy_cos <- function(t,r){
  if(abs(r)>=1) stop('r must be less than 1 in abs val')
  cos_a <- cos(2*pi*t)
  return(2*r^2*(cos_a - r^2)/(1-2*r^2*cos_a + r^4))
}

get_sin_cov <- function(t,t_prime,r,lambda = 1,nugget = 0.){
  diff_t_mat <- outer(t,t_prime,"-")
  
  lam_s = lambda*crazy_cos(0,r)
  
  out_mat <- lam_s^(-1)*crazy_cos(diff_t_mat,r)
  
  if(length(t)==length(t_prime)){
    out_mat <- out_mat + nugget*diag(length(t))
  }
  return(out_mat)
}

(test <- get_sin_cov(x1,x1,1))


draw_prior_gp <- function(t,r,lam,nugget = 10E-5){
  
  cov_mat <- get_sin_cov(t,t,r=r,lambda=lam,nugget = nugget)
  # Note: without 1.e-8*diag(length(grid)) nugget term cov_mat may not be positive definite,
  #      which is required to do the following Cholesky decomposition
  L  <- t(chol(cov_mat))
  
  path = L%*%rnorm(length(t))
  plot(t2,path,type = 'l',xlab = 'x',
       main =bquote('Prior Sample Path: '~lambda~'='~.(lam)~', r = '~.(r)))
  return(path)
  
}

system.time(prior_draw <- draw_prior_gp(t2,r=.9,lam=10,nugget = 10E-12))


draw_post_gp <- function(t,t_p,r,N,lam=1,y = NULL, beta1 = 3,beta2 = 4,
                         closed_form = FALSE, ...){
  
  if(!closed_form){
    obs_cov_mat <- get_sin_cov_open(t,t,r,N,lam,nugget = 10E-4)
    pred_cov_mat <- get_sin_cov_open(t_p,t_p,r,N,lam,nugget = 10E-4)
    cross_cov  <- get_sin_cov_open(t,t_p,r,N,lam)
  }else{
    obs_cov_mat <- get_sin_cov(t,t,r,lam,nugget = 10E-4)
    pred_cov_mat <- get_sin_cov(t_p,t_p,r,lam,nugget = 10E-4)
    cross_cov  <- get_sin_cov(t,t_p,r,lam)
  }

  
  if(is.null(y)){
    y = dbeta(t,beta1,beta2)
    title_plot = paste0("Predicting Beta(",beta1,",",beta2,")")
  }else{
    title_plot = "Predicting Distribution"
  }  
  
  cond_mean <- rep(1,length(t_p)) +  
    t(cross_cov) %*% solve(obs_cov_mat) %*% (y - rep(1,length(t)))
  
  cond_var <- pred_cov_mat - t(cross_cov) %*% solve(obs_cov_mat) %*% cross_cov
  
  plot(t,y,pch = 19,ylim = c(0,max(cond_mean + 2*sqrt(diag(cond_var)))), ylab = 'f(t)', 
       main = title_plot)
  lines(t_p,cond_mean,col = 'blue')
  lines(t_p,cond_mean + 1.96*sqrt(diag(cond_var)),col = 'red',lty = 2)
  lines(t_p,cond_mean - 1.96*sqrt(diag(cond_var)),col = 'red',lty = 2)
  abline(h = 0)
  
  L <- t(chol(cond_var))
  return(list(mean = cond_mean,var = cond_var))
}

x1 <- seq(0.1,.99,length = 7)
#plot(x1,y1)

t2 <- seq(0,.999,length= 1001)

system.time(post_draw <- draw_post_gp(x1,t2,r=.2,N=100,lam = 10,closed_form = TRUE))

system.time(post_draw <- draw_post_gp(Aj$V1,t2,y=Aj$V2,r=.4,100,.3))
