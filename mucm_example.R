##Making an emulator
#From first example of MUCM

truth <- function(x){
  return(3*x + cos(5*x))
}

#Covariance function
cov_exp <- function(d,ell,alp){
  return(exp(-(abs(d/ell))^alp))
}
cov_matern <- function(d,ell,nu){
  return(Matern(d,smoothness = nu,range = ell/sqrt(2*nu)))
}

build_cov_mat <- function(x_vals,x_vals_p,ell,lambda,type = "exponential",nu = 3/2){
  return(lambda^{-1}*build_R_mat(x_vals,x_vals_p,ell,type = type,nu = nu))
}

build_R_mat <- function(x_vals,x_vals_p,ell, alpha = 2, nugget = 1E-4,
                        type = "exponential",nu = 3/2){
  out_mat <- matrix(numeric(length(x_vals)*length(x_vals_p)),ncol = length(x_vals_p))
  for(i in 1:length(x_vals)){
    for(j in 1:length(x_vals_p)){
      d = abs(x_vals[i]-x_vals_p[j])
      if(type=="exponential"){
        out_mat[i,j] <- cov_exp(d,ell,alp = alpha)
      }else if(type =="matern"){
        out_mat[i,j] <- cov_matern(d,ell,nu)
      }else{
        stop("Type must be either exponential or matern")
      }
    }
  }
  if(dim(out_mat)[1]==dim(out_mat)[2]){
 #   out_mat <- out_mat + nugget*diag(dim(out_mat)[1])
  }
  return(out_mat)
}

#######################
##Simulating Process###
#######################

sim_process <- function(n,x_vals,ell,lambda, mean_vec = NULL){
  n_x <- length(x_vals)
  Y_mat <- matrix(numeric(n*n_x),ncol = n)
  if(is.null(mean_vec)){mean_vec = rep(0,n_x)}
  for(k in 1:n){
    
    cov_mat <- build_cov_mat(x_vals,x_vals,ell,lambda)
    #S <- chol(cov_mat)
    #y <- mean_vec + t(S)%*%rnorm(n_x,0,1)
    y <- mvrnorm(1,mean_vec, cov_mat)
    Y_mat[,k] <- y
  }
  return(Y_mat)
}


plot_gp <- function(x_vals, Y, ...){
  n = dim(Y)[2]
  for(i in 1:n){
    if(i==1){
      plot(x_vals,Y[,i],type ='l',...)
    }
    else{
      lines(x_vals,Y[,i],...)
    }
  }
}


#Making the first plot
par(mfrow = c(2,2))
x_vals = seq(-1,1,length.out = 100)
ell = .1
lambda = .1
gp_prior <- sim_process(5,x_vals,ell,lambda)
d01_s10 <- plot_gp(x_vals=x_vals,gp_prior,ylab = "f(x)",xlab = "x",
                   main = paste0("ell = ",ell," Lambda = ",lambda),ylim = c(-20,20))

ell = 1
lambda = .1
gp_prior <- sim_process(5,x_vals,ell,lambda)
d1_s10 <- plot_gp(x_vals=x_vals,gp_prior,ylab = "f(x)",xlab = "x",
                  main = paste0("ell = ",ell," Lambda = ",lambda),ylim = c(-20,20))

ell = .1
lambda = .01
gp_prior <- sim_process(5,x_vals,ell,lambda)
d01_s100 <- plot_gp(x_vals=x_vals,gp_prior,ylab = "f(x)",xlab = "x",
                    main = paste0("ell = ",ell," Lambda = ",lambda),ylim = c(-20,20))

ell = 1
lambda = .01
gp_prior <- sim_process(5,x_vals,ell,lambda)
d1_s100 <- plot_gp(x_vals=x_vals,gp_prior, ylab = "f(x)",xlab = "x",
                   main = paste0("ell = ",ell," Lambda = ",lambda),ylim = c(-20,20))





par(mfrow = c(1,1))
###############################################
##Inference, with Uncertainty in Beta, Lambda##
###############################################
x_vals = seq(-1,1,length.out = 101)
d1 = c(-1,sample(x_vals,5),1)
d1 = c(-1,-2/3,-1/3,0,1/3,2/3,1)

Y1 = t(t(truth(d1)))#make it a 7x1 matrix
n = length(d1)
p = 2

R <- build_R_mat(d,d,ell,type = "matern")
beta_hat <- solve(t(D)%*%solve(R)%*%D)%*%t(D)%*%solve(R)%*%Y

#log-likehood ells
ell_mle <- function(d,true_fun,R_type = "exponential",nu = 3/2, alpha = 2,
                    mucm_version = FALSE,plotit = FALSE,...){
  p = 2
  n = length(d)
  Y = t(t(true_fun(d)))
  #tau_vec <- seq(-6,1,by = 0.001)
  ell_vec <- seq(0.001,1, by = 0.0001)
  loglike <- numeric(length(ell_vec))
  D = cbind(rep(1,length(d)),d)
  
  for(k in 1:length(ell_vec)){
    R <- build_R_mat(d,d,ell_vec[k],type = R_type,nu = nu, alpha = alpha)
    beta_hat <- solve(t(D)%*%solve(R)%*%D)%*%t(D)%*%solve(R)%*%Y
    b = t(Y)%*%solve(R)%*%Y - t(Y)%*%solve(R)%*%D%*%beta_hat
    loglike[k] <- -((n-p)/2)*log(b) -
      .5*log(det(R)) - .5*log(det(t(D)%*%solve(R)%*%D))
  }
  #tau <- tau_vec[which.max(loglike)]
  ell <- ell_vec[which.max(loglike)]
  if(plotit){
    plot(ell_vec,loglike,type = 'l', xlab = '\u2113', ylab = 'p(\u2113|Y)', 
         ...)
    abline(v = ell,col = 'red', lwd = 2)
    legend('bottomleft',legend = bquote(hat('\u2113')~ " = " ~ .(ell)),lwd = 2,col = 'red')
  }
  return(ell)
}

GP_emulator1d <- function(x,d,true_fun = truth, ell=NULL,R_type = "exponential", nu = 3/2,alpha = 2,
                          ...){
  p = 2
  n = length(d)
  X = cbind(rep(1,length(x)),x)
  D = cbind(rep(1,length(d)),d)
  Y = t(t(true_fun(d)))
  if(is.null(ell)){
    ell <- ell_mle(d,true_fun = true_fun,R_type = R_type, nu = nu,alpha = alpha)
  }
  
  #To make the end plots
  R <- build_R_mat(d,d,ell,type = R_type, nu = nu,alpha = alpha)
  
  beta_hat <- solve(t(D)%*%solve(R)%*%D)%*%t(D)%*%solve(R)%*%Y #MLE/mean of Beta posterior
  print(beta_hat)
  lam_hat <- as.numeric((n-p)/(t(Y)%*%solve(R)%*%Y - t(Y)%*%solve(R)%*%D%*%beta_hat))
  
  
  lam <- rgamma(1,(n-p)/2,(t(Y)%*%solve(R)%*%Y - t(Y)%*%solve(R)%*%D%*%beta_hat)/2)
  bet <- t(rmvnorm(1,beta_hat,solve(lam*t(D)%*%solve(R)%*%D)))
  
  
  #Covariance 
  sig_11 = build_R_mat(x_vals,x_vals,ell,type = R_type, nu = nu,alpha = alpha)
  sig_12 = build_R_mat(x_vals,d,ell,type = R_type, nu = nu,alpha = alpha)
  sig_21 = build_R_mat(d,x_vals,ell,type = R_type, nu = nu,alpha = alpha)
  
  
  #If you want to "integrate Beta out"
  mu_star <- X%*%bet + sig_12%*%solve(R)%*%(Y - D%*%bet) #lambdas cancel
  
  #If you want to just use Beta MLE
  mu_star_static <- X%*%beta_hat + sig_12%*%solve(R)%*%(Y - D%*%beta_hat) #lambdas cancel
  sig_star <- lam_hat^{-1}*(sig_11-sig_12%*%solve(R)%*%sig_21)
  
  #Plotting things
  plot(x_vals,mu_star_static,type = 'l',ylim = c(-5,4),ylab = 'f(x)',main = "Predictive Distribution",
       xlab = 'x',...)
  points(d,Y,pch = 19)
  lines(x_vals,truth(x_vals),col = 'red')
  
  lines(x_vals, mu_star_static + 1.96*sqrt(round(diag(sig_star),10)),col = 'blue',lty = 2,lwd = 2)
  lines(x_vals, mu_star_static - 1.96*sqrt(round(diag(sig_star),10)),col = 'blue',lty = 2, lwd = 2)
  lines(x_vals,truth(x_vals),col = 'red')
  legend('topleft',c("Truth","Post Mean","95% Cred Interval"),lty = c(1,1,2),
         col = c('black','red','blue'))
  
  ellhat = eval(bquote(hat('\u2113')~ " = " ~ .(ell)))
  param = ifelse(R_type == "matern", eval(bquote(expression(nu ~ "=" ~ .(nu)))),
                                         eval(bquote(expression(alpha ~ "=" ~ .(alpha)))))
  #param = paste("nu: ",nu)
  cov = paste("Cov. Fun:",ifelse(R_type=="exponential","Pwr Exp","Mat\u{E9}rn"))
  
  #cov = paste0("Cov. Fun: ", R_type)
  legend('bottomright',c(cov,param,ellhat))
  return(list(mu_pred = mu_star_static,cov_pred = sig_star))
}

test1d = GP_emulator1d(x_vals,d1,R_type = "matern",nu = 51/2,alpha = 1.5)
plot(diag(test1d$cov_pred),type = 'l')

