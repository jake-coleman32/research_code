#Hierarchical GP
#Individual updates
library(TruncatedNormal)
library(mvtnorm)
library(dplyr)



GP_cov <- function(d,lambda,ell,nugget = 0.){
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(d))^2/ell^2) + nugget*diag(length(d))
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell){
  inds <- 1:length(d)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(d,d_star)))[inds,-inds]^2/ell^2)
  
  return(out_mat)
}

GP_cross_cov(d_scaled[i],d_scaled[-i],lam_cur[n],ell_cur[n])

#Data stuff
read_hist_data <- function(hist_loc, q_file){
  q_vals <- read.table(q_file)[,1]
  
  Y_mat <- lapply(as.list(q_vals),function(x){
    read.table(file = paste0(hist_loc,"q_",as.character(x),".dat"))[,2]
  }) %>%
    do.call(rbind,.)
  
  bin_left = read.table(file = paste0(hist_loc,"q_",as.character(q_vals[1]),".dat"))[,1]
  
  colnames(Y_mat) <- as.character(bin_left)
  
  return(Y_mat)
  
}

jet_path_comp <- "/home/grad/jrc71/Documents/Research/Computer_Emulation/JETSCAPE/JETSCAPE-STAT"
jet_path_lap <- "/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/"
hist_folder <- "q_dat_100k/"

current_path = jet_path_lap

q_vals_file <- paste0(current_path,"qhat_vals_100k.dat")


#Data
num_counts <- 100000
holdout <- 7
Y <- read_hist_data(hist_loc = paste0(current_path,hist_folder),
                    q_file = q_vals_file) *num_counts
Y <- Y[-holdout,]
d <- read.table(q_vals_file)[-holdout,1]
(d_scaled <- (d-min(d))/(max(d)-min(d)))
t_star <- 0.5
trunc_cols <- which(as.numeric(colnames(Y))<t_star)
Y_trunc <- Y[,trunc_cols]


I <- dim(Y)[1]
J <- dim(Y_trunc)[2]
alpha <- c(as.numeric(colnames(Y_trunc)),t_star)#Somewhat sketch
jitter = FALSE
if(jitter){
  alpha[-c(1,J+1)] <- alpha[-c(1,J+1)] + rnorm(J-1,0,1E-3)
}
b <- (alpha - lag(alpha))[-1]

#Parameter things
r <- 0.5
N <- (J-1) #One fewer than number of bins: the max number for N
#More in B matrix than A matrix
r_vec <- c(r^(1:floor(N/2)),r^(1:ceiling(N/2)))
R <- diag(r_vec)

make_coefs <- function(){
  A <- matrix(0,floor(N/2),J)
  B <- matrix(0,ceiling(N/2),J)
  for(n in 1:floor((N/2))){
    for(j in 2:(J+1)){
      A[n,j-1] <- t_star*(sin(2*pi*n*alpha[j]/t_star) - sin(2*pi*n*alpha[j-1]/t_star))/(2*pi*n)
      B[n,j-1] <- t_star*(cos(2*pi*n*alpha[j-1]/t_star) - cos(2*pi*n*alpha[j]/t_star))/(2*pi*n)
    }
  }
  if(N%%2){
    n = ceiling(N/2)
    for(j in 2:(J+1)){
      B[n,j-1] <- t_star*(cos(2*pi*n*alpha[j-1]/t_star) - cos(2*pi*n*alpha[j]/t_star))/(2*pi*n)
    }
  }
  
  return(list(A=A,B=B))
}

A <- make_coefs()$A
B <- make_coefs()$B
C <- cbind(t(A),t(B))

show_ranks = TRUE
if(show_ranks){
  print(paste("A",rankMatrix(t(A))[1]))
  print(paste("B",rankMatrix(t(B))[1]))
  print(paste("C",rankMatrix(t(C))[1]))
}

#Starting place for X
#Use fact that the sum of probabilities add to one
#Could be issue if we trim before observations end
#I.e. if gam != 1
j_trim <- 1
C_trim <- C[-j_trim,]
P_hat <- t(Y_trunc[,-j_trim]/num_counts) #P_hat is (J-1)xI

#Different starting points
X_mle <- sqrt(0.5)*solve(R)%*%solve(C_trim)%*%(P_hat - replicate(I,b[-j_trim]/t_star))
X_0 <- replicate(12,rep(0,9))



#Priors
lam_a <- rep(1,N)
lam_b <- rep(1,N)

ell_a <- rep(1,N)
ell_b <- rep(1,N)



hier_gp_mh_i <- function(iters = 1E4, burnin_prop = 0.1,
                       delta_lam = rep(0.3,N),
                       delta_ell = rep(0.3,N),
                       X_kap = replicate(N,rep(1E-1,I)), #column n is diagonal of proposal
                       #for GP n
                       verbose = FALSE
                       ){
  burnin <- iters*burnin_prop
  
  X_array <- array(0,dim = c(iters,N,I)) #Columns are GPs, rows are histograms
  lam_mat <- ell_mat <- matrix(0,iters,N)
  
  #Current values
  X_cur <- X_mle #currently N x I
  ell_cur <- rgamma(N,ell_a,ell_b)
  lam_cur <- rgamma(N,lam_a,lam_b)
  
  p_cur <- t(sqrt(2)*C%*%sweep(X_cur,1,r_vec,"*") + replicate(I,b/t_star))
  l_cur <- sum(Y_trunc*log(p_cur))# + sum(Y_trunc[,j_trim]*log(1-apply(p_cur,1,sum)))
  
  ell_acc <- lam_acc  <- numeric(N)
  x_acc <- matrix(0,N,I)
  
  time_start <- proc.time()
  for(t in 1:(burnin + iters)){
    for(n in 1:N){
      if(!t%%100){
        if(verbose){
          flush.console()
          cat("\r t = ", t, "Elapsed Time: ",as.numeric((proc.time() - time_start)[3])) 
          #print(paste0("t = ", t, ", n = ", n))
        }
      }
      
      #Update X_n
      auto_reject = FALSE
      x_cur_n <- X_cur[n,]
      
      #Need to find upper and lower bounds for each histogram
    #  l <- -Inf*(numeric(I) + 1)
     # u <- Inf*(numeric(I) + 1)
      for(i in 1:I){
        
        
        #All the constraints
        constr <- (-b/t_star - sqrt(2)*C[,-n]%*%(r_vec[-n]*X_cur[-n,i]))/(sqrt(2)*r_vec[n]*C[,n])
        
        #Lower constraints
        (a_constr <- constr[which(round(C[,n],10)>0)])
        
        #Upper constraints
        (b_constr <- constr[which(round(C[,n],10)<0)])
        
        #Bad constraints
        (bad_constr <- which(round(C[,n],10)==0))
        
        #If any bad constraints
        if(length(bad_constr)){
          for(j in bad_constr){
            
            #For each bad constraint, check if the rest of X[i,] take care of it
            #If they don't auto-reject the whole vector
            if(b[j]/t_star + sqrt(2)*C[j,-n]%*%(r_vec[-n]*X_cur[-n,i])<=0){
              print(paste("Bad times: i=",i,"n =",n,"j =",j))
              auto_reject = TRUE
              l_star = -Inf
            }
          }
        }
        
        #Update X_n[i]
        (l <- ifelse(length(a_constr),max(a_constr),-Inf))
        (u <- ifelse(length(b_constr),min(b_constr),Inf))
        x_cur_i <- x_cur_n[i]
        (x_st_i <- rtruncnorm(1,mean = x_cur_i,a = l, b = u,sd = X_kap[i,n]))
        
        (x_cond <- x_cur_n[-i])
        (sig_11 = 1/lam_cur[n])
        (sig_22 <- GP_cov(d_scaled[-i],lam_cur[n],ell_cur[n],nugget = 1E-5))
        (sig_12 <- t(matrix(GP_cross_cov(d_scaled[i],d_scaled[-i],lam_cur[n],ell_cur[n]))))
        
        (prior_mean = sig_12%*%(solve(sig_22)%*%x_cond))
        (prior_var = sig_11 - sig_12%*%solve(sig_22)%*%t(sig_12))
        
        X_star <- X_cur
        X_star[n,i] <- x_st_i
        
        
        p_star <- t(sqrt(2)*C%*%(sweep(X_star,1,r_vec,"*")) + replicate(I,b/t_star))
        if(sum(p_star<0)){
          print(paste("j = ",j,"i = ",i,"t = ",t,"n = ",n))
          stop("Dammit this shouldn't happen")
        }
        
        l_star <- sum(Y_trunc*log(p_star))
        
        adjust <- dtruncnorm(x=x_cur_i,a = l,b=u,mean = x_st_i,sd = X_kap[i,n]) -
          dtruncnorm(x=x_st_i,a = l,b=u,mean = x_cur_i,sd = X_kap[i,n])
        
        ratio <- l_star - l_cur +
          
          #Prior
          dnorm(x_st_i,mean = prior_mean, sd = sqrt(prior_var)) - 
          dnorm(x_cur_i,mean = prior_mean, sd = sqrt(prior_var)) + 
          
          #Adjustment ratio
          adjust
        
        if(length(ratio)>1){stop(paste("wtf n=",n))}
        if(runif(1)<exp(ratio)){
          x_acc[n,i] <- x_acc[n,i] + 1
          l_cur <- l_star
          X_cur[n,i] <- x_st_i
        }
        
        
      }# i loop
      
      #Calculate proposed likelihood based on proposed x_star_n
      
      #############
      #Update ell_n
      #############
      z_ell <- rnorm(1)
      ell_st <- ell_cur[n]*exp(delta_ell[n]*z_ell)
      adjust_ratio <- delta_ell[n]*z_ell
      
      #MH Ratio
      
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[n,],sigma = GP_cov(d_scaled,lam_cur[n],ell_st,nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[n,],sigma = GP_cov(d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(ell_st,ell_a[n],ell_b[n],log = TRUE) -
        dgamma(ell_cur[n],ell_a[n],ell_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
      
      if(ratio + rexp(1)>0){
        ell_acc[n] <- ell_acc[n] + 1
        ell_cur[n] <- ell_st
      }
      
      ############
      #Update lam_n
      ############
      z_lam <- rnorm(1)
      lam_st <- lam_cur[n]*exp(delta_lam[n]*z_lam)
      adjust_ratio <- delta_lam[n]*z_lam
      
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[n,],sigma = GP_cov(d_scaled,lam_st,ell_cur[n],nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[n,],sigma = GP_cov(d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
        #Priors
        dgamma(lam_st,lam_a[n],lam_b[n],log = TRUE)-
        dgamma(lam_cur[n],lam_a[n],lam_b[n],log = TRUE) +
        
        #Proposal Adjustment
        adjust_ratio
      
      
      if(ratio + rexp(1)>0){
        lam_acc[n] <- lam_acc[n] + 1
        lam_cur[n] <- lam_st
      }
      
      
    }# n loop
    
    if(t > burnin){
      X_array[t-burnin,,] <- X_cur
      lam_mat[t-burnin,] <- lam_cur
      ell_mat[t-burnin,] <- ell_cur
    }
  }# t loop
  
  return(list(X = X_array,lam = lam_mat,ell = ell_mat,
              x_acc = x_acc/(iters+burnin),
              lam_acc = lam_acc/(iters+burnin),
              ell_acc = ell_acc/(iters+burnin)))
}

hope_i <- hier_gp_mh_i(iters = 1E3,verbose = TRUE, burnin_prop = 0.3,
                    X_kap = matrix(nrow = I, byrow = FALSE,data = c(
                      rep(1E-3,I),#1
                      rep(1E-3,I),#2
                      rep(1E-3,I), #3
                      rep(1E-2,I), #4
                      rep(1E-3,I),#5
                      rep(1E-3,I),#6
                      rep(1E-2,I),#7
                      rep(1E-1,I),#8
                      rep(1E-1,I))) #9
)

hope_i$x_acc
apply(hope_i$x_acc,1,mean)

hope_i$lam_acc
hope_i$ell_acc

save(hope, file = "run_1k.Rdata")

stop("We're done here")

blah = 1:5
save(blah,file = "didnt_stop.Rdata")
#summaryRprof('Hope_against_hope.out')


X <- hope_i$X

x_means <- apply(X,c(2,3),mean)
x_means - X_mle

thin = 500
iters <- 1:dim(X)[1]
X_thin <- X[!iters%%thin,,]

i <- 4
X_i <- X_thin[,,i]

plot_traces(X_i,save_pics = FALSE)

gam=1
T_out=seq(0,t_star,length=75)
plot_dens_i(X_i)
plot(as.numeric(colnames(Y_trunc)),Y_trunc[i,])


meeting_parent <- '/Users/Jake/Dropbox/Research/Computer_Emulation/meetings/2017/'
meeting_folder <- 'meeting_3_9/'
path <- paste0(meeting_parent,meeting_folder)
save_pics = TRUE
suffix = '_indv_large'




plot_traces <- function(x_hist,save_pics = FALSE){
  
  for(n in 1:N){
    if(n<ceiling(N/2)){
      var_type = 'Z'
      index = n
    }
    else{
      var_type = 'W'
      index = n-floor(N/2)
    }
    if(save_pics) pdf(paste0(path,var_type,n,suffix,'.pdf'))
    plot(x_hist[,n],type = 'l',ylab = bquote(.(var_type)[.(index)]),
         main = bquote('Trace Plot of '~.(var_type)[.(index)])~', Histogram'~.(i),
         xlab = 'Iteration')
    if(save_pics) dev.off()
  }
}

plot_thetas <- function(save_pics = FALSE){
  for(n in 1:N){
    if(save_pics) pdf(file = paste0(path,'lam',n,suffix,'.pdf'))
    plot(hope_i$lam[,n],xlab = "Iteration",ylab = bquote(lambda[.(n)]),type = 'l',
         main = bquote('Trace Plot of '~lambda[.(n)]))
    if(save_pics) dev.off()
    
    if(save_pics) pdf(file = paste0(path,'ell',n,suffix,'.pdf'))
    plot(hope_i$ell[,n],xlab = "Iteration",ylab = bquote('ell'[.(n)]),type = 'l',
         main = bquote('Trace Plot of ell'[.(n)]))
    if(save_pics) dev.off()
  }
}

plot_thetas(save_pics = TRUE)


est_dens_i <- function(x_mat, r){
  f_mat <-matrix(0,dim(x_mat)[1],length(T_out))
  Nx <- dim(x_mat)[2]
  for(t in 1:length(T_out)){
    cos_vec <- cos(2*pi*T_out[t]*(1:floor(Nx/2))/t_star)*r^(1:floor(Nx/2))
    sin_vec <- sin(2*pi*T_out[t]*(1:ceiling(Nx/2))/t_star)*r^(1:ceiling(Nx/2))
    cos_sin_vec <- c(cos_vec,sin_vec)
    
    f_mat[,t] <- sqrt(2)*t(matrix(cos_sin_vec))%*%t(x_mat) + gam/t_star
  }
  return(f_mat)
}


plot_dens_i <- function(x_mat,r=0.5, save_pics = FALSE,legend_side = 'topright',...){
  f_est <- est_dens_i(x_mat,r)
  mean_est <- apply(f_est,2,mean)
  
  if(save_pics) pdf(paste0(meeting_parent,meeting_folder,'mh_dens',suffix,'.pdf'))
  plot(T_out,mean_est,type = 'l', main = 'GP Density Estimate',
       ylab = 'Density',xlab = 'y',lwd = 2,
       ylim = c(0,max(apply(f_est,2,quantile,0.975))),
       ...)#,ylim = c(0,2))
  lines(T_out,apply(f_est,2,quantile,0.025),col = 'blue',lty=2)
  lines(T_out,apply(f_est,2,quantile,0.975),col = 'blue',lty=2)
  # legend(legend_side,c('Post Mean','Post 95% Cred'),
  #       lwd = 2,lty = c(1,2),col = c('black','blue'))
  if(save_pics) dev.off()
}


est_probs_i <- function(X_i,save_pics = FALSE){
  p_out <- sqrt(2)*C%*%sweep(t(X_i),1,r_vec,FUN = "*") + b/t_star
  return(p_out)
}

i = 5
plot(alpha[-J],apply(est_probs_i(X_i),1,mean),type = 'l')
lines(alpha[-J],apply(est_probs_i(X_i),1,quantile,probs = 0.025),lty = 2,col = 'blue')
lines(alpha[-J],apply(est_probs_i(X_i),1,quantile,probs = 0.975),lty = 2,col = 'blue')
points(as.numeric(colnames(Y_trunc)),apply(est_probs_i(i-1),1,mean))
points(as.numeric(colnames(Y_trunc)),Y_trunc[i,]/num_counts)


#Ok let's plot all data points
cols = rainbow(I)

if(save_pics) pdf(file = paste0(path,'all_data',suffix,'.pdf'))
plot(as.numeric(colnames(Y_trunc)),
     Y_trunc[1,]/num_counts,col = cols[1],type = 'l',
     xlab = 'Aj',
     ylab = 'Fraction of Bin',
     main = 'Data Across Input')
for(i in 2:I){
  lines(as.numeric(colnames(Y_trunc)),
        Y_trunc[i,]/num_counts,col = cols[i],type = 'l')
}
if(save_pics) dev.off()

#Now lets plot our estiamtes
#Doesn't vary more than the 95% confidence intervals
if(save_pics) pdf(file = paste0(path,'all_emulations',suffix,'.pdf'))
plot(as.numeric(colnames(Y_trunc)),
    apply(est_probs_i(X[,,1]),1,mean),col = cols[1],type = 'l',
     xlab = 'Aj',
     ylab = 'Probability of Bin',
     main = 'Emulated Values Across Input')
for(i in 2:I){
  lines(as.numeric(colnames(Y_trunc)),
        apply(est_probs_i(X[,,i]),1,mean),col = cols[i],type = 'l')
}
if(save_pics) dev.off() 

##Drawing from prior

#Draw lamda, ell
#Draw N GP realizations
#Keep if none of the P's are negative
##There is no chance of this happening

lam <- rgamma(N,lam_a,lam_b)
ell <- rgamma(N,ell_a,ell_b)

how_many=0
bad_p = TRUE
while(bad_p){
  how_many = how_many + 1
  X_prior <- matrix(0,N,I)
  for(n in 1:N){
    X_prior[n,] <- rmvnorm(1,sigma = GP_cov(d_scaled,lam[n],ell[n]))
  }
  
  P <- sqrt(2)*C%*%sweep(X_prior,1,r_vec) + replicate(I,b/t_star)
  if(!(sum(P)>0)) bad_p=FALSE
}



#To predict new d
#Use draws from X as conditioning values
#For each draw, draw from conditional normal - so N*500k cond. normal draws
#Do this for each X_n?
cov_inv_mats <-  vector("list",N)
inv_vec_mult <- array(0,dim=c(dim(X_thin)))#T x N x I

for(n in 1:N){
  cov_inv_mats[[n]] <- vector("list",dim(X_thin)[1])
  for(t in 1:dim(X_thin)[1]){
    lam_nt <- hope_i$lam[t,n]
    ell_nt <- hope_i$ell[t,n]
    cov_inv_mats[[n]][[t]] <- solve(GP_cov(d_scaled,lambda = lam_nt,ell = ell_nt,nugget = 1E-6))
    inv_vec_mult[t,n,] <- cov_inv_mats[[n]][[t]]%*%X_thin[t,n,]
  }
}

pred_p <- function(X,d_prime, d_cond=d_scaled,lam,ell, sig_22_inv_list = cov_inv_mats,
                   sig_inv_x = inv_vec_mult, verbose = FALSE){
  #Loop through N components, collect values
  #500k x 9 matrix
  
  n_iters <- dim(X)[1]
  pred_components <- matrix(0,n_iters,N)
  
  
  for(n in 1:N){
    pred_mean <- pred_var <- dim(X)[1]
    time_start <- proc.time()
    for(t in 1:n_iters){
      sig_22_inv <- sig_22_inv_list[[n]][[t]]
      sig_12 <- t(GP_cross_cov(d_cond, d_prime,lambda = lam[t,n],ell = ell[t,n]))
      
      pred_mean[t] <- sig_12%*%inv_vec_mult[t,n,]
      #pred_var[t] <- 1/lam[t,n] - sig_12%*%sig_22_inv%*%t(sig_12)
      #if(pred_var[t]<0)stop(paste("pred_var is",pred_var[t],"nt is",n,t)) 
      
    }
    
    
    #pred_components[,n] <- rnorm(n_iters,pred_mean,sqrt(pred_var))
    pred_components[,n] <- pred_mean
    
  }

  return(pred_components)
}

system.time(test <- pred_p(X=X_thin, d_prime = 0.5, d_cond = d_scaled,lam = hope_i$lam, 
               ell = hope_i$ell))
