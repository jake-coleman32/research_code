#GP Histogram toy example - not emulation
library(dplyr) #or lag doesn't work 
library(truncnorm)
library(TruncatedNormal)
library(lpSolve)

#Make the underlying data
dat_alph <- 7
dat_bet <-3
x_dat <- rbeta(10000,dat_alph,dat_bet)

#Make the histogram data
#My own function makes it easier to control bins
n_bins <- 10
my_hist <- function(dat,bounds){
  bins <- numeric(length(alpha)-1)
  for(i in 1:length(bins)){
    bins[i] <- length(which(dat>bounds[i] & dat <=bounds[i+1]))
  }
  return(bins)
}

alpha <- seq(0,1,length=(n_bins+1))
alpha[2:(length(alpha)-1)] <- alpha[2:(length(alpha)-1)] + 
  rnorm(length(alpha)-2,0,1E-3)

(y <- my_hist(x_dat,bounds=alpha))

#Make coefficient matrices
make_coef_mats <- function(N,alpha){
  J <- length(alpha)
  A <- B <- matrix(0,N,(J-1))
  
  for(n in 1:N){
    for(j in 2:J){
      A[n,j-1] <- (sin(2*pi*n*alpha[j]) - sin(2*pi*n*alpha[j-1]))/(2*pi*n)
      B[n,j-1] <- (cos(2*pi*n*alpha[j-1]) - cos(2*pi*n*alpha[j]))/(2*pi*n)
    }
  }
  
  return(list(A=A,B=B))
}
coef_5 <- make_coef_mats(5,alpha = alpha)
A = coef_5$A
B = coef_5$B




mh_hist <- function(r = 0.5,iters = 1E4,kap_z =rep(1E-1,N), kap_w = rep(1E-1,N)){
  #Write some MH sampler
  burnin <- 0.1*iters
  z_cur <- rnorm(N)
  w_cur <- rnorm(N)
  
  p_star <- p_cur <- numeric(J-1)
  p_cur <- sapply(seq_along(1:length(p_cur)),function(j){
    sum(r^(1:N)*(z_cur*A[,j] + w_cur*B[,j])) + 0.1
  })
  l_cur <- sum(y*log(p_cur))
  
  z_mat <- w_mat <- matrix(0,iters,N)
  z_acc <- z_ob <- w_acc <- w_ob <- numeric(N)
  
  for(i in 1:(iters+burnin)){
    #propose new values for Z and W
    
    for(n in 1:N){
      #Now W's
      w_star_n <- rnorm(1,w_cur[n],kap_w[n])
      w_star <- w_cur
      w_star[n] <- w_star_n
      
      p_star <- sapply(seq_along(1:length(p_star)),function(j){
        sum(2*r^(1:N)*(z_cur*A[,j] + w_star*B[,j] + bin_sizes[j]))
      })      
      
      if(sum(p_star<0)){
        l_star <- -Inf
        w_ob[n] <- w_ob[n] + 1
      }else{
        l_star <- sum(y*log(p_star))
      }
      
      ratio <- l_star + dnorm(w_star_n,log = TRUE) - 
        l_cur - dnorm(w_cur[n],log = TRUE)
      
      if(runif(1)<exp(ratio)){
        w_cur <- w_star
        l_cur <- l_star
        w_acc[n] <- w_acc[n] + 1
      }
      #else{
      #  print(l_star)
       # print(dnorm(w_star_n,log = TRUE))
      #  print(l_cur)
       # print(dnorm(w_cur[n],log = TRUE))
      #}
    }
    
    
    for(n in 1:N){
      #Z's first
      z_star_n <- rnorm(1,z_cur[n],kap_z[n])
      z_star <- z_cur
      z_star[n] <- z_star_n
      
      
      p_star <- sapply(seq_along(1:length(p_star)),function(j){
        sum(sqrt(2)*r^(1:N)*(z_star*A[,j] + w_cur*B[,j])) + bin_sizes[j]
      })      
      
      if(sum(p_star<0)){
        l_star = -Inf
        z_ob[n] <- z_ob[n] + 1
      }else{
        l_star <- sum(y*log(p_star))
      }
      ratio <- l_star + dnorm(z_star_n,log = TRUE) - 
        l_cur - dnorm(z_cur[n],log = TRUE)
      
      if(runif(1)<exp(ratio)){
        z_cur <- z_star
        l_cur <- l_star
        z_acc[n] <- z_acc[n] + 1
      }
    }

    
    
    if(i>burnin){
      z_mat[i-burnin,] <- z_cur
      w_mat[i-burnin,] <- w_cur
    }
    
  }
  print(z_acc/(iters+burnin))
  print(w_acc/(iters+burnin))
  print("ob's")
  print(z_ob)
  print(w_ob)
  
  return(list(z_mat = z_mat,w_mat = w_mat,z_acc = z_acc, w_acc = w_acc,
              z_ob=z_ob, w_ob = w_ob))
}

res <- mh_hist(r = 0.5,kap_z = c(rep(1E-10,N-1),1),kap_w = rep(1E-2,N))

plot(res$z_mat[,3],type = 'l')
plot(res$w_mat[,1],type = 'l')

est_dens <- function(z_mat,w_mat,T_out=seq(0,1,length=20)){
  f_mat <-matrix(0,dim(z_mat)[1],length(T_out))
  N <- dim(z_mat)[2]
  for(t in 1:length(T_out)){
    f_mat[,t] <- 2*z_mat%*%matrix(r^(1:N)*cos(2*pi*T_out[t]*(1:N))) + 
      2*w_mat%*%matrix(r^(1:N)*sin(2*pi*T_out[t]*(1:N)))
  }
  return(f_mat)
}

f_est <- est_dens(res$z_mat,res$w_mat)[-c(1:1000),]
mean_est <- apply(f_est,2,mean)
plot(T_out,mean_est,type = 'l')#,ylim = c(0,2))
lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)

plot(T_out,dbeta(T_out,dat_alph,dat_bet))
#















###################
#Doing a truncated normal update
####################

#Making C
coef_5 <- make_coef_mats(5,alpha = alpha)
A = coef_5$A
B = coef_5$B


C <- t(rbind(A,B))
C <- C[,-5]#All zeroes
Nx <- dim(C)[2]

J <- dim(C)[1]
r = 0.5
r_vec <- sqrt(2)*c(r^(1:floor(Nx/2)),r^(1:ceiling(Nx/2)))

non_j
phat <- y/sum(y)
C_trim <- C[-non_j,]
x_mle <- as.numeric(solve(C_trim)%*%(phat[-non_j]-0.1)/r_vec)
(bin_sizes <- (alpha-lag(alpha))[-1])

(y2 <- rbinom(length(y),size = y,p = 1E-2))

mh_hist_trunc <- function(r = 0.5,iters = 1E4, burnin_prop = 0.1, kap =rep(1E-1,Nx),
                          verbose = FALSE,xstart= 'MLE',
                          y = y){
  #Write some MH sampler
  burnin <- floor(burnin_prop*iters)
  p_star <- p_cur <- numeric(J)
  x_mat <-matrix(0,iters,Nx)
  x_acc <- x_ob <- numeric(Nx)
  
  
  if(xstart=="MLE"){x_cur <- x_mle}
  else{x_cur <- rnorm(Nx)}
  

  (p_cur <- as.numeric(C_trim%*%(r_vec*x_cur))+bin_sizes[-non_j])
  if(sum(p_cur<0)){
    stop("Unlucky: try again")
  }
  
  l_cur <- sum(y[-1]*log(p_cur)) + y[1]*log(1-sum(p_cur))
  
  for(i in 1:(iters+burnin)){
    #propose new values for Z and W
    for(n in 1:Nx){
      if(verbose){
        cat("\r i = ", i, ", n = ", n) 
        flush.console()
      }
      
      auto_reject <- FALSE
      constraints <- (-0.1-C[,-n]%*%(r_vec[-n]*x_cur[-n]))/
        (r_vec[n]*C[,n])
      (a_constraints <- constraints[which(round(C[,n],10)>0)])
      (b_constraints <- constraints[which(round(C[,n],10)<0)])
      (bad_constraints <- which(round(C[,n],10)==0))
      
      if(length(bad_constraints)){
        for(j in bad_constraints){
          if((bin_sizes[j]+C[j,-n]%*%(r_vec[-n]*x_cur[-n])) <= 0){
            auto_reject <- TRUE
            l_star <- -Inf
          }
        }
      }
      if(!auto_reject){
        x_star_n <- rtruncnorm(1,a=max(a_constraints),b = min(b_constraints),
                                mean = x_cur[n],sd = kap[n])
        x_star <- x_cur
        x_star[n] <- x_star_n
        
        
        (p_star <- as.numeric(C%*%(r_vec*x_star)) + bin_sizes)
        if(round(sum(p_star),10)!=1){stop(paste(i,", ",n))}
        #l_star <- sum(y[-1]*log(p_star)) + y[1]*log(1-sum(p_star))
        l_star <- sum(y*log(p_star))
       # print(l_star)
        ratio <- l_star + dnorm(x_star_n,log = TRUE) - 
          l_cur - dnorm(x_cur[n],log = TRUE) +
          log(dtruncnorm(x_cur[n],a = max(a_constraints),b = min(b_constraints),
                     mean = x_star_n,sd = kap[n])) -
          log(dtruncnorm(x_star_n,a = max(a_constraints),b = min(b_constraints),
                     mean = x_cur[n],sd = kap[n]))

        if(runif(1)<exp(ratio)){
          x_cur <- x_star
          l_cur <- l_star
          x_acc[n] <- x_acc[n] + 1
        }
      }else{#auto-reject
        x_ob[n] <- x_ob[n] + 1
      }

     # if(sum(p_star<0)){
      #  l_star = -Inf
       # x_ob[n] <- x_ob[n] + 1
      #}else{
      #}

    }
    if(i>burnin){
      x_mat[i-burnin,] <- x_cur
    }
  }
  print(x_acc/(iters+burnin))
  
  return(list(x_mat = x_mat,x_acc = x_acc/(iters + burnin), x_ob=x_ob))
}

res2 <- mh_hist_trunc(iters = 1E4,burnin_prop = 0.3, r = 0.5,#kap = c(2E-3,2E-2,3E-2,2E-2,5E-1,
                                           #       1,1E-2,4E-1,1))
                     kap = c(5E-3,#z1
                             5E-1,#z2
                             1E-1,#z3
                             1E-1,#z4
                             1E-2,#w1 - 5,
                             5E-3,#w2 -6
                             1E-2,#w3 - 7
                             1E-1,#w4 -8
                             1E-1),#w5 -9
                     #kap = rep(0,Nx),
                     verbose = FALSE,xstart = 'MLE',
                     y = y)
plot_traces(save_pics = FALSE)

meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_2_9/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = FALSE
suffix = '_reg_mle'

save(res2,file=paste0(path,'res',suffix,'.Rdata'))

apply(res2$x_mat,2,mean)
x_mle


if(save_pics){
  mean_table <- round(cbind('Post_Mean' = apply(res2$x_mat,2,mean),
                      'MLE' = x_mle),4)
  write.table(mean_table,file = paste0(path,'mean_table',suffix,'.txt'))
}

xtable(mean_table,digits = 3)

plot_traces <- function(save_pics = FALSE){
  for(i in 1:Nx){
    if(i<ceiling(Nx/2)){
      var_type = 'Z'
      index = i
    }
    else{
      var_type = 'W'
      index = i-floor(Nx/2)
    }
    if(save_pics) pdf(paste0(path,var_type,i,suffix,'.pdf'))
    plot(res2$xmat[,i],type = 'l',ylab = bquote(.(var_type)[.(index)]),
         main = bquote('Trace Plot of '~.(var_type)[.(index)]),
         xlab = 'Iteration')
    if(save_pics)dev.off()
  }
}







est_dens2 <- function(x_mat,T_out=seq(0,1,length=20),r = 0.5){
  f_mat <-matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    cos_vec <- cos(2*pi*T_out[t]*1:floor(Nx/2))*r^(1:floor(Nx/2))
    sin_vec <- sin(2*pi*T_out[t]*1:ceiling(Nx/2))*r^(1:ceiling(Nx/2))
    cos_sin_vec <- c(cos_vec,sin_vec)
    
    f_mat[,t] <- sqrt(2)*x_mat%*%matrix(cos_sin_vec) + 1
  }
  return(f_mat)
}

T_out=seq(0,1,length=100)

f_est <- est_dens2(res2$xmat,T_out = T_out)
mean_est <- apply(f_est,2,mean)


if(save_pics) pdf(paste0(meeting_parent,meeting_folder,'mh_dens',suffix,'.pdf'))
plot(T_out,mean_est,type = 'l', main = 'GP Density Estimate',
     ylab = 'Density',xlab = 'y',lwd = 2)#,ylim = c(0,2))
lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)
lines(T_out,dbeta(T_out,dat_alph,dat_bet),type ='l',col = 'red',lwd = 2)
legend('topleft',c('Post Mean','Post 95% Cred','Truth'),
       lwd = 2,lty = c(1,2,1),col = c('black','blue','red'))
if(save_pics) dev.off()


##########################
#Michael's Augmented Gibbs Sampler
###########################
check_rank <- function(C){
  rankifremoved <- sapply(1:ncol(C), function (x) rankMatrix(C[,-x])[1])
  return(rankifremoved)
}
#Get C matrix
coef_6 <- make_coef_mats(6,alpha = alpha)
A <- coef_6$A
B <- coef_6$B

A <- t(A)
View(round(sweep(A,1,A[,6],"/"),4))

C <- t(rbind(coef_6$A,coef_6$B)) %>%
  rbind(rnorm(dim(C)[2],0,5E-2))
#View(round(C,4))
rankMatrix((C))[1]

#Delete column of zeroes, keep track of index
#Add noise to make a square matrix
#Follow package stuff
#Build Gibbsy
A <- coef_5$A
B <- coef_5$B
C <- t(rbind(A,B))
C <- C[,-5]
C_trim <- C[-1,]

gibbs_hist <- function(r=0.5, iters = 1E4){
  
  burnin <- 0.1*iters
  p_star <- p_cur <- numeric(J)
  x_mat <-matrix(0,iters,Nx)
  x_acc <- x_ob <- numeric(Nx)
  
  r_vec <- sqrt(2)*c(r^(1:floor(Nx/2)),r^(1:ceiling(Nx/2)))
  
  #x_cur <- (solve(C_trim))%*%(phat[-1]-0.1)/r_vec
  
  #p_cur <- as.numeric(C%*%diag(r_vec)%*%(x_cur))+bin_sizes
  
  u <- runif(length(phat[-1]),0,phat[-1]^log(y[-1])*y[-1])
  Ap <- C%*%diag(r_vec)
  Ap_trim <- C_trim%*%diag(r_vec)
  App <- t(Ap)%*%Ap
  for(i in 1:iters){
    bounds <- u - bin_sizes[-1]
    #lower <- solve(t(Ap)%*%Ap)%*%t(Ap)%*%bounds
    Y <- mvrandn(l = bounds, u = rep(Inf,Nx),Sig = Ap_trim%*%t(Ap_trim), n = 1)
    x_star <- solve(C_trim)%*%Y
    #x_star <- mvrandn(l = lower, u = rep(Inf,Nx),Sig = diag(Nx), n = 1)
    (p_star <- as.numeric(Ap_trim%*%x_star+bin_sizes[-1]))
    (l_x <- p_star^log(y[-1])*y[-1])
    u <- runif(length(p_star),0,l_x)

    x_mat[i,] <- x_star
  }
  return(x_mat)
}

test <- gibbs_hist()

for(j in 1:9){
  print(lp(direction = 'max',objective.in = Ap[j,],const.mat = Ap[-j,],
     const.dir = '>',const.rhs = bounds[-j]))
}


#Look at posterior estimates for 



