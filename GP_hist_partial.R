#GP Partial Histogram Example
library(dplyr) #or lag doesn't work 
library(truncnorm)
library(TruncatedNormal)
library(lpSolve)

#Make the underlying data
dat_alph <- 3
dat_bet <-7
x_dat <- rbeta(10000,dat_alph,dat_bet)

#Make the histogram data
#My own function makes it easier to control bins
my_hist <- function(dat,bounds){
  bins <- numeric(length(alpha)-1)
  for(i in 1:length(bins)){
    bins[i] <- length(which(dat>bounds[i] & dat <=bounds[i+1]))
  }
  return(bins)
}

t_star = 1

bin_size = 0.1
(alpha <- seq(0,t_star, length = 50))
alpha[2:(length(alpha)-1)] <- alpha[2:(length(alpha)-1)] + 
  rnorm(length(alpha)-2,0,1E-3)

(y <- my_hist(x_dat,bounds=alpha) + 1)#so x_mle doesn't freak out

#Make coefficient matrices
make_coef_mats <- function(N){
  J <- length(alpha)
  A <- B <- matrix(0,N,(J-1))
  
  for(n in 1:N){
    for(j in 2:J){
      A[n,j-1] <- t_star*(sin(2*pi*n*alpha[j]/t_star) - 
                            sin(2*pi*n*alpha[j-1]/t_star))/(2*pi*n)
      B[n,j-1] <- t_star*(cos(2*pi*n*alpha[j-1]/t_star) - 
                           cos(2*pi*n*alpha[j]/t_star))/(2*pi*n)
    }
  }
  
  return(list(A=A,B=B))
}
coef_5 <- make_coef_mats(24)
A = coef_5$A
B = coef_5$B

C <- t(rbind(A,B))
#C <- C[,-3]#All zeroes
(Nx <- dim(C)[2])

(J <- dim(C)[1])
(bin_sizes <- (alpha-lag(alpha))[-1])

rankMatrix(C)[1]
r = 0.5

#non_j = 1
#phat <- (phat <- y/sum(y))
#C_trim <- C[-non_j,]
#x_mle <- as.numeric(solve(C_trim)%*%(phat[-non_j]-0.1)/r_vec)


mh_hist_trunc <- function(r = 0.5,iters = 1E4, burnin_prop = 0.1, kap =rep(1E-1,Nx),
                          verbose = FALSE,
                          gam = pbeta(t_star,dat_alph,dat_bet),
                          y = y){
  #Write some MH sampler
  burnin <- floor(burnin_prop*iters)
  p_star <- p_cur <- numeric(J)
  x_mat <-matrix(0,iters,Nx)
  x_acc <- x_ob <- numeric(Nx)
  
  r_vec <- sqrt(2)*c(r^(1:floor(Nx/2)),r^(1:ceiling(Nx/2)))
  
  x_cur <- rnorm(Nx)
  
  
  (p_cur <- as.numeric(C%*%(r_vec*x_cur))+bin_sizes*gam/t_star)
  if(sum(p_cur<0)){
    stop("Unlucky: try again")
  }
  
  l_cur <- sum(y*log(p_cur))
  
  for(i in 1:(iters+burnin)){
    #propose new values for Z and W
    for(n in 1:Nx){
      if(verbose){
        flush.console()
        cat("\r i = ", i, ", n = ", n) 
      }
      
      auto_reject <- FALSE
      constraints <- (-bin_sizes*gam/t_star-C[,-n]%*%(r_vec[-n]*x_cur[-n]))/
        (r_vec[n]*C[,n])
      (a_constraints <- constraints[which(round(C[,n],10)>0)])
      (b_constraints <- constraints[which(round(C[,n],10)<0)])
      (bad_constraints <- which(round(C[,n],10)==0))
      if(length(bad_constraints)){
        for(j in bad_constraints){
          if((bin_sizes[j]*gam/t_star+C[j,-n]%*%(r_vec[-n]*x_cur[-n])) <= 0){
            auto_reject <- TRUE
            l_star <- -Inf
          }
        }
      }
      if(!auto_reject){
        bot_const <- ifelse(length(a_constraints)>0,max(a_constraints),
                            -Inf)
        top_const <- ifelse(length(b_constraints)>0,min(b_constraints),
                            Inf)
        x_star_n <- rtruncnorm(1,a=bot_const,b = top_const,
                               mean = x_cur[n],sd = kap[n])
        x_star <- x_cur
        x_star[n] <- x_star_n
        
        
        (p_star <- as.numeric(C%*%(r_vec*x_star)) + bin_sizes*gam/t_star)
        if(sum(p_star<0)){stop(paste("neg pstar: i=",i,", n=",n))}
        
        if(round(sum(p_star),10)!=gam){stop(paste("pstar doesn't add up: i=",i,", n=",n))}
        
        l_star <- sum(y*log(p_star))
        ratio <- l_star + dnorm(x_star_n,log = TRUE) - 
          l_cur - dnorm(x_cur[n],log = TRUE) +
          log(dtruncnorm(x_cur[n],a = bot_const,b = top_const,
                         mean = x_star_n,sd = kap[n])) -
          log(dtruncnorm(x_star_n,a = bot_const,b = top_const,
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
      
    }#n loop
    
    if(i>burnin){
      x_mat[i-burnin,] <- x_cur
    }
  }#i loop
  print(x_acc/(iters+burnin))
  
  return(list(x_mat = x_mat,x_acc = x_acc/(iters + burnin), x_ob=x_ob))
}

res2 <- mh_hist_trunc(iters = 1E4,burnin_prop = 0.3, r = 0.5,#kap = c(2E-3,2E-2,3E-2,2E-2,5E-1,
                      #       1,1E-2,4E-1,1))
#                       kap = c(5E-2,#z1
#                               1E-1,#z2
#                               1E-1,#z3
#                               1E-1,#z4
#                               1E-1,#z5
#                               4E-1,#w1 -6
#                               1E-1,#w2 -7
#                               1E-1,#w3 -8
#                               1E-1,#w4-9
#                               1E-1#w5 -10
#                               ),
                      kap = rep(1E-1,Nx),
                      verbose = FALSE,
                      #gam = sum(y/length(x_dat)),
                      y = y)

plot_dens(res2,r=0.5,T_out,save_pics = FALSE)
lines(density(x_dat),col = 'green')

plot_traces(save_pics = FALSE)

meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_2_16/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = FALSE
suffix = '_5'

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
    plot(res2$x_mat[,i],type = 'l',ylab = bquote(.(var_type)[.(index)]),
         main = bquote('Trace Plot of '~.(var_type)[.(index)]),
         xlab = 'Iteration')
    if(save_pics)dev.off()
  }
}







est_dens2 <- function(x_mat, r,T_out=seq(0,t_star,length=20)){
  f_mat <-matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    cos_vec <- cos(2*pi*T_out[t]*(1:floor(Nx/2))/t_star)*r^(1:floor(Nx/2))
    sin_vec <- sin(2*pi*T_out[t]*(1:ceiling(Nx/2))/t_star)*r^(1:ceiling(Nx/2))
    cos_sin_vec <- c(cos_vec,sin_vec)
    
    f_mat[,t] <- sqrt(2)*x_mat%*%matrix(cos_sin_vec) + gam/t_star
  }
  return(f_mat)
}

T_out=seq(0,t_star,length=100)



plot_dens <- function(results,r,T_out, save_pics = FALSE,...){
  f_est <- est_dens2(results$x_mat,r,T_out = T_out)
  mean_est <- apply(f_est,2,mean)
  
  if(save_pics) pdf(paste0(meeting_parent,meeting_folder,'mh_dens',suffix,'.pdf'))
  plot(T_out,mean_est,type = 'l', main = 'GP Density Estimate',
       ylab = 'Density',xlab = 'y',lwd = 2,
       ylim = c(0,max(apply(f_est,2,quantile,0.975))),
       ...)#,ylim = c(0,2))
  lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
  lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)
  lines(T_out,dbeta(T_out,dat_alph,dat_bet),type ='l',col = 'red',lwd = 2)
  legend('topleft',c('Post Mean','Post 95% Cred','Truth'),
         lwd = 2,lty = c(1,2,1),col = c('black','blue','red'))
  if(save_pics) dev.off()
}

summary(apply(f_est,1,function(x){trapz(T_out,x)}))
