#GP Partial Histogram Example
library(dplyr) #or lag doesn't work 
library(truncnorm)
library(TruncatedNormal)
library(lpSolve)
library(Matrix)

sig_to_c <- function(sig,r=0.5){
  return(sig*sqrt(1-r^2)/sqrt(2*r^2))
}
c_to_sig <- function(c,r=0.5){
  return(c*sqrt(2*r^2/(1-r^2)))
}
prior_prob <- function(sig,t_st = 0.5,t_min = 0){
  return(pnorm(-(1/(t_st-t_min))/sig))
}

needed_sig <- function(prob,r,t_st = 0.5,t_min=0){
  return(-1/(qnorm(prob)*(t_st-t_min)))
}



set.seed(47)
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

t_minus = 0.1
t_star = 0.6

#bin_size = 0.1
(alpha <- seq(t_minus,t_star, length = 11))

jitter = FALSE
if(jitter){
  alpha[2:(length(alpha)-1)] <- alpha[2:(length(alpha)-1)] + 
    rnorm(length(alpha)-2,0,1E-3)
}
(y <- my_hist(x_dat,bounds=alpha))# + 1)#so x_mle doesn't freak out



#Make coefficient matrices
#Coef for cos(2*pi*nt), sin(2*pi*nt), 1
coef_mats_2pi <- function(N){
  J <- length(alpha)
  A <- B <- matrix(0,N,(J-1))
  
  for(n in 1:N){
    for(j in 2:J){
      A[n,j-1] <- sqrt(2)*(t_star-t_minus)*(sin(2*pi*n*(alpha[j] - t_minus)/(t_star - t_minus)) - 
                                              sin(2*pi*n*(alpha[j-1] - t_minus)/(t_star - t_minus)))/(2*pi*n)
      B[n,j-1] <- sqrt(2)*(t_star - t_minus)*(cos(2*pi*n*(alpha[j-1] - t_minus)/(t_star - t_minus)) - 
                                                cos(2*pi*n*(alpha[j] - t_minus)/(t_star-t_minus)))/(2*pi*n)
    }
  }
  
  return(list(A=A,B=B))
}

coef_mats_pi <- function(N){
  J <- length(alpha)
  A <- B <- matrix(0,N,(J-1))
  
  for(n in 1:N){
    for(j in 2:J){
      A[n,j-1] <- sqrt(2)*(t_star-t_minus)*(sin(pi*n*(alpha[j] - t_minus)/(t_star - t_minus)) - 
                                              sin(pi*n*(alpha[j-1] - t_minus)/(t_star - t_minus)))/(pi*n)
      B[n,j-1] <- sqrt(2)*(t_star - t_minus)*(cos(pi*n*(alpha[j-1] - t_minus)/(t_star - t_minus)) - 
                                                cos(pi*n*(alpha[j] - t_minus)/(t_star-t_minus)))/(pi*n)
    }
  }
  
  return(list(A=A,B=B))}

coef_mats_cos_only <- function(N){
  J <- length(alpha)
  C <- matrix(0,N,(J-1))
  for(n in 1:N){
    for(j in 2:J){
      C[n,j-1] <- sqrt(2)*(t_star-t_minus)*
        (sin(pi*n*(alpha[j] - t_minus)/(t_star - t_minus)) - 
           sin(pi*n*(alpha[j-1] - t_minus)/(t_star - t_minus)))/(pi*n)
    }
  }
  return(C)
}
coef_mats_cos_only_odd <- function(N){
  J <- length(alpha)
  C <- matrix(0,N,(J-1))
  for(n in 1:N){
    for(j in 2:J){
      C[n,j-1] <- sqrt(2)*(t_star-t_minus)*
        (sin(pi*(n/2 + 1)*(alpha[j] - t_minus)/(t_star - t_minus)) - 
           sin(pi*(n/2 + 1)*(alpha[j-1] - t_minus)/(t_star - t_minus)))/(pi*(n/2 + 1))
    }
  }
  return(C)
  
}

N <- 20

basis_type = "cos_only"
if(basis_type=="2pi"){
  print("Using 2pi basis/covariance")
  coef_mats <- coef_mats_2pi(N)
  C <- t(rbind(coef_mats$A,coef_mats$B))
  
}else if(basis_type =="pi"){
  print("Using pi basis/covariance")
  
  coef_mats <- coef_mats_pi(N)
  C <- t(rbind(coef_mats$A,coef_mats$B))
  
}else if(basis_type=="cos_only"){
  print("Using cosine-only basis/covariance")
  
  coef_mats <- coef_mats_cos_only(N)
  C <- t(coef_mats)
}

(bad_cols <- which(apply(C,2,function(x){sum(abs(x)<1E-5)})==dim(C)[1]))
if(length(bad_cols)) C <- C[,-bad_cols]



dim(C)
rankMatrix(C)[1]


(Nx <- dim(C)[2])

(J <- dim(C)[1])
(bin_sizes <- (alpha-lag(alpha))[-1])

r = 0.95
(sig <- needed_sig(1E-5,t_st = t_star, t_min = t_minus))#Check these functions are right
prior_prob(sig,t_star,t_minus)

(c <- sig_to_c(sig = sig,r=r))

#This won't cause problems for the cos-only vector, because we'll never use Nw
##and bad_cols will only be in 1:N
Nz <- Nw <-  1:N
if(length(bad_cols)){
  bad_z <- bad_cols[which(bad_cols<=N)]
  bad_w <- (bad_cols-N)[which(bad_cols>N)]
  if(length(bad_z)) Nz <- Nz[-bad_z]
  if(length(bad_w)) Nw <- Nw[-bad_w]
}

if(basis_type %in% c("2pi","pi")){
  (r_vec <- c*c(r^Nz,r^Nw))
}else{
  r_vec <- c*r^Nz
}


#non_j = 1
#phat <- (phat <- y/sum(y))
#C_trim <- C[-non_j,]
#(x_mle <- as.numeric(solve(C_trim)%*%(phat[-non_j]-0.1)/r_vec))

(gam = pbeta(t_star,dat_alph,dat_bet) - pbeta(t_minus,dat_alph,dat_bet))

mh_hist_trunc <- function(r = 0.5,iters = 1E4, burnin_prop = 0.1, kap =rep(1E-1,Nx),
                          verbose = FALSE, x_start  = "random",
                          y = y){
  #Write some MH sampler
  burnin <- floor(burnin_prop*iters)
  p_star <- p_cur <- numeric(J)
  x_mat <-matrix(0,iters,Nx)
  x_acc <- x_ob <- numeric(Nx)
  
  
  
  if(x_start == "random") x_cur <- rnorm(Nx)
  else if(x_start== "mle") x_cur <- x_mle
  
  print("x_start is")
  print(x_cur)
  
  
  
  sum((p_cur <- as.numeric(C%*%(r_vec*x_cur))+bin_sizes*gam/(t_star - t_minus)))
  if(sum(p_cur<0)){
    stop("Unlucky: try again")
  }
  
  l_cur <- sum(y*log(p_cur))
  
  time_start <- proc.time()
  for(i in 1:(iters+burnin)){
    #propose new values for Z and W
    if(!i%%(iters*0.1)){
      flush.console()
      cat("\r i = ", i, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
      #print(paste0("t = ", t, ", n = ", n))
    }
    
    for(n in 1:Nx){
      if(verbose){
        if(!(i%%(0.1*iters)))
          flush.console()
        cat("\r i = ", i, ", n = ", n) 
      }
      
      auto_reject <- FALSE
      constraints <- (-bin_sizes*gam/(t_star - t_minus)-C[,-n]%*%(r_vec[-n]*x_cur[-n]))/
        (r_vec[n]*C[,n])
      (a_constraints <- constraints[which(round(C[,n],10)>0)])
      (b_constraints <- constraints[which(round(C[,n],10)<0)])
      (bad_constraints <- which(round(C[,n],10)==0))
      if(length(bad_constraints)){
        for(j in bad_constraints){
          if((bin_sizes[j]*gam/(t_star-t_minus) + C[j,-n]%*%(r_vec[-n]*x_cur[-n])) <= 0){
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
        (x_star_n <- rtruncnorm(1,a=bot_const,b = top_const,
                               mean = x_cur[n],sd = kap[n]))
        x_star <- x_cur
        x_star[n] <- x_star_n
        
        
        (p_star <- as.numeric(C%*%(r_vec*x_star)) + bin_sizes*gam/(t_star - t_minus))
        if(sum(p_star<0)){stop(paste("neg pstar: i=",i,", n=",n))}
        
        #if(round(sum(p_star),10)!=gam){stop(paste("pstar doesn't add up: i=",i,", n=",n))}
        
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
  print("\n")
  print(x_acc/(iters+burnin))
  
  return(list(x_mat = x_mat,x_acc = x_acc/(iters + burnin), x_ob=x_ob))
}

kap_x <- c(1E-1#1
           ,1E-1#2
           ,1E-1#3
           ,1#4
           ,1#5
           ,1#6
           ,1#7
           ,1#8
           ,1#9
           ,1#10,
           ,1#11,
           ,1#12,
           ,1#13,
           ,1#14,
           ,1#15,
           ,1#16,
            ,1#17,
            ,1#18,
           #            5#19,
)

res2 <- mh_hist_trunc(iters = 1E5,burnin_prop = 0.3, 
                      #kap = rep(1,Nx),
                      kap = kap_x,
                      x_start = "random",
                      verbose = FALSE,
                      #gam = sum(y/length(x_dat)),
                      y = y)
res2$x_ob
apply(res2$x_mat,2,mean)
acf(res2$x_mat[,1])

#Thinning if necessary
thin = 20
iters <- 1:dim(res2$x_mat)[1]
x_thin <- res2$x_mat[!iters%%thin,]

#Checking easy-to-check
effectiveSize(x_thin)
acf(x_thin[,1])

T_out=seq(t_minus,t_star,length=100)

plot_dens(x_thin,T_out,save_pics = FALSE,legend_side = 'topright', plot_beta = TRUE,
          normalize = FALSE)
add_hist(y)
abline(h = 0)
lines(density(x_dat),col = 'green')

plot_coefs(x_thin)
plot_traces(x_thin,save_pics = FALSE)


meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_2_23/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = FALSE
suffix = '_49_bins'

save(res2,file=paste0(path,'res',suffix,'.Rdata'))

apply(res2$x_mat,2,mean)
x_mle


if(save_pics){
  mean_table <- round(cbind('Post_Mean' = apply(res2$x_mat,2,mean),
                            'MLE' = x_mle),4)
  write.table(mean_table,file = paste0(path,'mean_table',suffix,'.txt'))
}

xtable(mean_table,digits = 3)

plot_traces <- function(x_mat,save_pics = FALSE){
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
    plot(x_mat[,i],type = 'l',ylab = bquote(X[.(i)]),
         main = bquote('Trace Plot of '~X[.(i)]),
         xlab = 'Iteration')
    if(save_pics)dev.off()
  }
}
plot_coefs <- function(xmat){
  means <- apply(xmat,2,mean)
  plot(means,pch = 19, col = 'blue', xlab = "Component",ylab = 'Value',
       main = 'Components 95% Credible Intervals',
       ylim = c(min(apply(xmat,2,quantile,probs = 0.025)),
                max(apply(xmat,2,quantile,probs = 0.975))))
  abline(h = 0,col = 'red')
  
  arrows(1:dim(xmat)[2], apply(xmat,2,quantile,probs = 0.025),
         1:dim(xmat)[2], apply(xmat,2,quantile,probs = 0.975), 
         length=0.05, angle=90, code=3)
  
}


est_dens2 <- function(x_mat,T_out=seq(0,t_star,length=20),normalize = FALSE){
  f_mat <- matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    
    if(basis_type=="2pi"){
      cos_vec <- sqrt(2)*cos(2*pi*Nz*(T_out[t]-t_minus)/(t_star-t_minus))*r^Nz
      sin_vec <- sqrt(2)*sin(2*pi*Nw*(T_out[t]-t_minus)/(t_star-t_minus))*r^Nw
      cos_sin_vec = c(cos_vec,sin_vec)
    }else if(basis_type=="pi"){
      cos_vec <- sqrt(2)*cos(pi*Nz*(T_out[t]-t_minus)/(t_star-t_minus))*r^Nz
      sin_vec <- sqrt(2)*sin(pi*Nw*(T_out[t]-t_minus)/(t_star-t_minus))*r^Nw
      cos_sin_vec = c(cos_vec,sin_vec)
    }else if(basis_type=="cos_only"){
      cos_sin_vec <- sqrt(2)*cos(pi*Nz*(T_out[t]-t_minus)/(t_star-t_minus))*r^Nz
    }else(stop("You need to pick a basis type, big fella"))
    
    f_mat[,t] <- c*x_mat%*%matrix(cos_sin_vec) + gam/(t_star-t_minus) 
  }
  if(normalize){
    f_mat[which(f_mat<0)] = 0
    f_mat <-  t(apply(f_mat,1,function(x) gam*x/trapz(T_out,x)))
  }
  return(f_mat)
}


plot_dens <- function(x_mat,T_out, save_pics = FALSE,legend_side = 'topright', normalize = FALSE,
                      plot_beta = FALSE,...){
  f_est <- est_dens2(x_mat,T_out = T_out, normalize = normalize)
  mean_est <- apply(f_est,2,mean)
  
  if(save_pics) pdf(paste0(meeting_parent,meeting_folder,'mh_dens',suffix,'.pdf'))
  plot(T_out,mean_est,type = 'l', 
       main = bquote('GP Density Estimate, Bins = '~.(length(alpha)-1)~"& "~N[x]~"= "~.(Nx)),
       # main = "GP Density Estimate",
       ylab = 'Density',xlab = 'y',lwd = 2,
       ylim = c(0,max(apply(f_est,2,quantile,0.975))),
       ...)#,ylim = c(0,2))
  if(plot_beta) lines(T_out,dbeta(T_out,dat_alph,dat_bet),type ='l',col = 'red',lwd = 2)
  
  lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
  lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)
  legend(legend_side,c('Post Mean','Post 95% Cred','Truth'),
         lwd = 2,lty = c(1,2,1),col = c('black','blue','red'))    
  if(save_pics) dev.off()
}
add_hist <- function(Y = Y_new_trunc,add = TRUE){
  #Requires alpha
  
  Y_fake_dat <- c()
  for(j in 1:length(Y)){
    Y_fake_dat <- c(Y_fake_dat, runif(Y[j], alpha[j],alpha[j+1]))
  }
  hist(Y_fake_dat, breaks = alpha,probability = 1,add = add,
       col = rgb(1,0,0,0.1))
}

summary(apply(f_est,1,function(x){trapz(T_out,x)}))

check_rank <- function(C){
  rankifremoved <- sapply(1:ncol(C), function (x) rankMatrix(C[,-x])[1])
  return(rankifremoved)
}

#Scractch
Y_small_bins <- numeric(length(Y_new_trunc)/2)
for(j in 1:(length(Y_new_trunc)/2)){
  k1 = 2*j-1
  print(k1)
  k2 = 2*j
  print(k2)
  Y_small_bins[j] <- Y_new_trunc[k1] + Y_new_trunc[k2]
}

(y <- Y_new_trunc[1:10])

