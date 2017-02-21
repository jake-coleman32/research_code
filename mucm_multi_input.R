library(fields)

###########LHD Functions#################
latin_hypercube <- function(n,p){
  u_mat <- b_mat <- x_mat <- matrix(numeric(n*p),ncol = n)
  for(i in 1:p){
    u_mat[i,] <- runif(n)
    b_mat[i,] <- sample(0:(n-1))
  }
  x_mat = (b_mat + u_mat)/n
  return(x_mat)
}

opt_latin_hypercube <- function(n,p,iters,crit){
  cubes <- vector("list",iters)
  vals <- numeric(iters)
  for(k in 1:iters){
    cubes[[k]] <- latin_hypercube(n,p)
    vals[k] <- crit(cubes[[k]])
  }
  return(cubes[[which.max(vals)]])
}

maximin <- function(D){
  #D is a p x n matrix
  cur_min = Inf
  for(i in 1:(dim(D)[2]-1)){
    for(j in (i+1):dim(D)[2]){
      cd = ifelse(dim(D)[1]==1,abs(D[,i]-D[,j]),dist(D[,i]-D[,j]))
      cur_min = min(cd,cur_min)
    }
  }
  return(cur_min)
}
#####################################
##Begin 2D input example##########
#################################
#Based on http://mucm.aston.ac.uk/MUCM/MUCMToolkit/index.php?page=ExamCoreGP2Dim.html

x1 = c(0.86,	0.27,	0.57,	0.5,	0.14,	0.23,	0.39,	0.93,	0.98,	0.54,
       0.77,	0.81,	0.18,	0.71,	0.31,	0.04,	0.69,	0.4,	0.63,	0.08,
       0.34,	0.04,	0.83,	0.21,	0.67,	0.47,	0.14,	0.51,	0.75,	0.95)

#x1 = (1420-1370)*x1 + 1370

x2 = c(0.7,	0.97,	0.29,	0.77,	0.73,	0.11,	0.21,	0.45,	0.05,	0.9,
       0.52,	0.33,	0.56,	0.84,	0.49,	0.09,	0.19,	0.37,	0.64,	0.86,
       0.45,	0.58,	0.75,	0.84,	0.37,	0.99,	0.15,	0.66,	0.07,	0.24)
#x2 = (0.4-0.2)*x2 + 0.2


plot(x1,x2,main = "Training Inputs",
     xlab = bquote("Solar Constant (W/"*m^2*")"), ylab = "Albedo (%)", pch = "+")
D1 = cbind(rep(1,length(x1)),x1,x2)
n = dim(D1)[1]
p = dim(D1)[2]

y = c(11.81,	-4.95,	25.51,	5.27,	4.85,	30.28,	27.50,	20.89,	34.50,	0.43,
      17.59,	25.27,	13.39,	3.67,	16.50,	29.96,	29.78,	21.14,	12.78,	-0.98,
      17.99,	12.11,	9.10,	1.03,	22.71,	-4.61,	28.33,	11.60,	34.10,	29.28)
Y1 = t(t(y))


#Covariance function
cov_exp <- function(d,ells,alp){
  return(exp(-(abs(d/ells))^alp))
}
cov_matern <- function(d,ells,nu){
  return(Matern(d,smoothness = nu,range = ells/sqrt(2*nu)))
}


multi_R_exp <- function(X1,X2, ells,type = "exponential",alp = 2,nu = 3/2,nug = 1E-7){
  n1 = dim(X1)[1]
  n2 = dim(X2)[1]
  out_mat <- matrix(numeric(n1*n2),ncol = n2)
  for(i in 1:n1){
    for(k in 1:n2){
      d = abs(X1[i,] - X2[k,])
      if(type =="exponential"){
        out_mat[i,k] = prod(exp(-(abs(d/ells))^alp))
      }else if(type=="matern"){
        out_mat[i,k] = prod(Matern(d,smoothness = nu,range = ells/sqrt(2*nu)))
      }
    }
  }
  if(n1 == n2){
    out_mat <- out_mat + nug*diag(n1)
  }
  return(out_mat)
}

tau_post <- function(params,D,Y,alp=2,nu=3/2,type = "exponential"){
  n = dim(D)[1]
  p = dim(D)[2]
  taus = params
  ells = exp(taus/2)
  R = multi_R_exp(D[,-1],D[,-1],ells,type,alp,nu)
  beta_hat <- solve(t(D)%*%solve(R)%*%D)%*%t(D)%*%solve(R)%*%Y
  b = t(Y)%*%solve(R)%*%Y - t(Y)%*%solve(R)%*%D%*%beta_hat
  neg_log_post = ((n-p)/2)*log(b) + .5*as.numeric(determinant(R, logarithm=TRUE)$modulus) +
    .5*as.numeric(determinant(t(D)%*%solve(R)%*%D,logarithm=TRUE)$modulus)
  return(neg_log_post)
}

GP_emulator <- function(D,yd,pred_X,R_type = "exponential",alp = 2, nu = 3/2){
  if(R_type=="exponential"){
    print(paste("Kernel is exponential with alpha =",alp))
  }else if(R_type=="matern"){
    print(paste("Kernel is Matérn with nu =",nu))
  }else{
    stop("R_type must be either exponential or matern")
  }
  YD = t(t(yd))
  n = dim(D)[1]
  p = dim(D)[2]
  mle_output = optim(par = rep(0,p-1),f = tau_post,D = D, Y = YD,
                       type = R_type, nu = nu,alp = alp)
  ell_hat <- exp(mle_output$par/2)
  
  if(mle_output$convergence){
    print(paste("Convergence failure:",mle_output$convergence))
  }
  
  R <- multi_R_exp(D[,-1],D[,-1],ell_hat, type = R_type,nu = nu,alp = alp)
  (beta_hat <- solve(t(D)%*%solve(R)%*%D)%*%t(D)%*%solve(R)%*%YD)
  (lam_hat <- as.numeric((n-p)/(t(YD)%*%solve(R)%*%YD - t(YD)%*%solve(R)%*%D%*%beta_hat)))
  
  new_X <- cbind(rep(1,dim(pred_X)[1]),pred_X)

  R11  <- multi_R_exp(new_X[,-1],new_X[,-1],ell_hat, type = R_type,nu = nu,alp = alp)
  R12 <-  multi_R_exp(new_X[,-1],D1[,-1],ell_hat, type = R_type,nu = nu,alp = alp)
  R21 <- multi_R_exp(D1[,-1],new_X[,-1],ell_hat, type = R_type,nu = nu,alp = alp)
  b <- t(YD)%*%solve(R)%*%YD - t(YD)%*%solve(R)%*%D%*%beta_hat
  
  mu <- new_X%*%beta_hat + R12%*%solve(R)%*%(YD-D%*%beta_hat)
  cov <- as.numeric(b/(n-p))*(R11 - R12%*%solve(R)%*%R21)
  
  return(list(mu_predict = mu,cov_predict = cov,ell_hat=ell_hat))
}

validate <- function(D,yd,X,R_type = "exponential",alpha = 2, nu = 3/2,val_y,...){
  gp_output <- GP_emulator(D,yd,X,R_type,alpha, nu)
  mu <- gp_output$mu_predict
  cov <- gp_output$cov_predict

  
  errors = (val_y - mu)/sqrt(diag(cov))
  plot(errors,ylim = c(-4,4),xlab = "Validation point index",ylab = "Errors",pch = 4,
       main = "Individual Standardized Errors",...)
  abline(h =c(-1.96,1.96))
  
  param = ifelse(R_type == "matern", eval(bquote(expression(nu ~ "=" ~ .(nu)))),
                 eval(bquote(expression(alpha ~ "=" ~ .(alpha)))))
  cov = paste("Cov. Fun:",ifelse(R_type=="exponential","Pwr Exp","Mat\u{E9}rn"))
  
  legend('topright',c(cov,param))
}


##Validation##
v1 = c(0.00,	0.83,	0.16,	0.37,	0.76,	0.64,	0.58,	0.91,	0.42,	0.21)
#v1 = (1420-1370)*v1 + 1370

v2 = c(0.12,	0.83,	0.51,	0.64,	0.44,	0.70,	0.46,	0.03,	0.24,	0.94)
#v2 = (0.4-0.2)*v2 + 0.2

vy = c(28.66,	4.68,	15.04,	11.63,	20.39,	10.83,	18.78,	34.76,	26.60,	-4.08)
n_v = length(v1)

plot(v1,v2,main = "Training Inputs With Validation",xlab = bquote("Solar Constant (W/"*m^2*")"),
     ylab = "Albedo (%)", pch = as.character(1:n_v), col = 'green4')
legend('topright',c("Training","Validation"),pch = c(3,19), col = c("black","green4"),cex = .8)

validate(D1,y,cbind(v1,v2),R_type = "matern",val_y = vy,alpha = 1.9,nu = 5/2)


#Ind. Standardized Errors
(errors = (vy - mu)/sqrt(diag(cov)))
plot(errors,ylim = c(-4,4),main = "Individual Standardized Errors",xlab = "Validation point index",ylab = "Errors",pch = 4)
abline(h =c(-1.96,1.96))

#Mahalanobis Distance - not sure where MUCM got this
(M <- t(newY - mu)%*%solve(cov)%*%(newY-mu))
(var_M <- 2*n_v*(n_v+n-p-2)/(n-p-4))#does not match up with MUCM.
d1 = n_v
d2 = n-p
2*d2^2*(d1+d2-2)/(d1*(d2-2)^2*(d2-4))


#New Emulator
D2 = rbind(D1,newX)
Y2 = rbind(Y1,newY)

(ell_hat2 <- exp(optim(par = c(-1,-1),f = tau_post,D = D2, Y = Y2)$par/2))
R2 <- multi_R_exp(D2[,-1],D2[,-1],ell_hat)
(beta_hat <- solve(t(D2)%*%solve(R2)%*%D2)%*%t(D2)%*%solve(R2)%*%Y2)

