#Hierarchical GP
#Individual updates
library(truncnorm)
library(mvtnorm)
library(dplyr)
library(caTools)
library(Matrix)


GP_cov <- function(d,lambda,ell,nugget = 0.){
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(d))^2/ell^2) + nugget*diag(length(d))
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell){
  inds <- 1:length(d)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(d,d_star)))[inds,-inds]^2/ell^2)
  
  return(out_mat)
}

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

jet_path_comp <- "/home/grad/jrc71/Documents/Research/Computer_Emulation/JETSCAPE/JETSCAPE-STAT/"
jet_path_lap <- "/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/"
hist_folder <- "q_dat_100k/"

on_comp = TRUE
#Create new directory based on date/time, change to it
if(on_comp){
  model_time = Sys.time()
  new_dir = paste0(jet_path_comp,'/',model_time)
  dir.create(new_dir)
  setwd(new_dir)
}else{
  stop('Fix boolean')
}


current_path = jet_path_comp

q_vals_file <- paste0(current_path,"qhat_vals_100k.dat")


#Data
num_counts <- 100000
holdout <- 7
Y <- read_hist_data(hist_loc = paste0(current_path,hist_folder),
                    q_file = q_vals_file) *num_counts
Y_new <- Y[holdout,]
Y <- Y[-holdout,]
d <- read.table(q_vals_file)[-holdout,1]
d_new <-read.table(q_vals_file)[holdout,1]
(d_scaled <- (d-min(d))/(max(d)-min(d)))
d_new_s <- (d_new-min(d))/(max(d)-min(d))
t_star <- 0.5
trunc_cols <- which(as.numeric(colnames(Y))<t_star)
Y_trunc <- Y[,trunc_cols]
Y_new_trunc <- Y_new[trunc_cols]

#data and parameters
data <- list()
data$num_counts <- num_counts
data$holdout <- holdout
data$Y <- Y
data$Y_new <- Y_new
data$Y_trunc <- Y_trunc
data$Y_new_trunc <- Y_new_trunc
data$t_star <- t_star
data$d_scaled <- d_scaled
data$d_new_s <- d_new_s
data$trunc_cols <-trunc_cols



I <- dim(Y)[1]
J <- dim(Y_trunc)[2]
alpha <- c(as.numeric(colnames(Y_trunc)),t_star)#Somewhat sketch
jitter = FALSE
if(jitter){
  alpha[-c(1,J+1)] <- alpha[-c(1,J+1)] + rnorm(J-1,0,1E-3)
}
b <- (alpha - lag(alpha))[-1]

data$I <- I
data$J <- J
data$alpha <- alpha
data$b <- b

#Save the list of data objects 
save_data = FALSE #Change this if you change something
if(save_data){
  save(data, file = "data_list.Rdata")
}

#Parameter things
params <- list()

r <- 0.5
N <- (J-1) #One fewer than number of bins: the max number for N
#More in B matrix than A matrix
Nz <- 1:floor(N/2)
Nw <- 1:ceiling(N/2)

c <- 1
r^max(Nw)
pnorm(-1/sqrt((2*c*r^2/(1-r^2))))
r_vec <- c(c*r^Nz,c*r^Nw)
#r_vec <- c(Nz*r^Nz,Nw*r^Nw)
#r_vec <- c((r^Nz)/Nz^2,(r^Nw)/Nw^2)
R <- diag(r_vec)

params$c <- c
params$N <- N
params$r <- r
params$Nz <- Nz
params$Nw <- Nw
params$r_vec <- r_vec

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

params$C <- C

#Save the parameters to use later
save(params,file = "params_list.Rdata")

run_description <- "r_vec is c*r^n, with r =0.5,c=1. Lambda is fixed at 1. Checking new code structure "
write(run_description,file="model_description.txt")


#Check the ranks of our coefficient matrices
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
X_0 <- replicate(I,rep(0,N))



#Priors
lam_a <- rep(1,N)
lam_b <- rep(1,N)

ell_a <- rep(1,N)
ell_b <- rep(1,N)



hier_gp_mh_i <- function(iters = 1E4, burnin_prop = 0.1,
                         delta_lam = rep(0,N),
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
  lam_cur <- rep(1,N)
  
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

hope_i <- hier_gp_mh_i(iters = 5E5,verbose = TRUE, burnin_prop = 0.3,
                       X_kap = matrix(nrow = I, byrow = FALSE,data = c(
                         rep(1E-1,I),#1
                         rep(1E-1,I),#2
                         rep(1E-1,I), #3
                         rep(1E-1,I), #4
                         rep(1E-1,I),#5
                         rep(1E-1,I),#6
                         rep(1E-1,I),#7
                         rep(1E-1,I),#8
                         rep(1E-1,I))) #9
)


save(hope, file = "sampler_vals.Rdata")




