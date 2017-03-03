#Hierarchical GP
library(TruncatedNormal)
library(mvtnorm)
library(dplyr)
library(truncnorm)



GP_cov <- function(d,lambda,ell,nugget = 0.){

  out_mat <- lambda^(-1)*exp(-as.matrix(dist(d))^2/ell^2) + nugget*diag(length(d))
  
  return(out_mat)
}

GP_cross_cov <- function(d,d_star,lambda,ell){
  inds <- 1:length(d)
  
  out_mat <- lambda^(-1)*exp(-as.matrix(dist(c(d,d)))[inds,-inds]^2/ell^2)
  
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

jet_path_comp <- "/home/grad/jrc71/Documents/Research/Computer_Emulation/JETSCAPE/JETSCAPE-STAT"
jet_path_lap <- "/Users/Jake/Dropbox/Research/JETSCAPE/JETSCAPE-STAT/"
hist_folder <- "q_dat_100k/"

current_path = jet_path_lap
setwd(current_path)
q_vals_file <- paste0(current_path,"qhat_vals_100k.dat")


#Data
num_counts <- 100000
Y <- read_hist_data(hist_loc = paste0(current_path,hist_folder),
                    q_file = q_vals_file) *num_counts
d <- read.table(q_vals_file)[,1]
(d_scaled <- (d-min(d))/(max(d)-min(d)) )
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
R <- diag(c(r^(1:floor(N/2)),r^(1:ceiling(N/2))))

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



hier_gp_mh <- function(iters = 1E4, burnin_prop = 0.1,
                       delta_lam = rep(0.3,N),
                       delta_ell = rep(0.3,N),
                       X_kap = replicate(N,rep(1E-1,I)), #column n is diagonal of proposal
                                                        #for GP n
                       verbose = FALSE
                       ){
  burnin <- iters*burnin_prop
  
  X_array <- array(0,dim = c(iters,I,N)) #Columns are GPs, rows are histograms
  lam_mat <- ell_mat <- matrix(0,iters,N)

  #Current values
  X_cur <- t(X_mle)
  ell_cur <- rgamma(N,ell_a,ell_b)
  lam_cur <- rgamma(N,lam_a,lam_b)
  
  p_cur <- t(sqrt(2)*C%*%R%*%t(X_cur) + replicate(I,b/t_star))
  l_cur <- sum(Y_trunc*log(p_cur))# + sum(Y_trunc[,j_trim]*log(1-apply(p_cur,1,sum)))
  
  ell_acc <- lam_acc <- x_acc <- numeric(N)
  
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
      x_cur_n <- X_cur[,n]
      
      #Need to find upper and lower bounds for each histogram
      l <- -Inf*(numeric(I) + 1)
      u <- Inf*(numeric(I) + 1)
      for(i in 1:I){
        
        
        #All the constraints
        constr <- (-b/t_star - sqrt(2)*C[,-n]%*%R[-n,-n]%*%X_cur[i,-n])/(sqrt(2)*R[n,n]*C[,n])
        
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
            if(b[j]/t_star + sqrt(2)*C[j,-n]%*%R[-n,-n]%*%X_cur[i,-n]<=0){
              print(paste("Bad times: i=",i,"n =",n,"j =",j))
              auto_reject = TRUE
              l_star = -Inf
            }
          }
        }
        
        l[i] <- ifelse(length(a_constr),max(a_constr),-Inf)
        u[i] <- ifelse(length(b_constr),min(b_constr),Inf)
        
      }# i loop
      
      #Draw truncated normal proposal
      #See help file for instruction
      x_star_n <- mvrandn(l = l - x_cur_n,
                          u = u - x_cur_n, 
                          Sig = diag(X_kap[,n]),
                          n = 1) + 
        x_cur_n
      
      #Calculate proposed likelihood based on proposed x_star_n
      X_star <- X_cur
      X_star[,n] <- x_star_n
      
      
      p_star <- t(sqrt(2)*C%*%R%*%t(X_star) + replicate(I,b/t_star))
      if(sum(p_star<0)){
        print(paste("j = ",j,"i = ",i,"t = ",t,"n = ",n))
        stop("Dammit this shouldn't happen")
      }
      
      l_star <- sum(Y_trunc*log(p_star))
      
      
      #This takes forever, so I'm not including it right now
      #adjusted <- dtmvnorm(x_cur_n,mean = x_star_n,sigma = diag(X_kap[,n]),lower = l,upper = u,log = TRUE) -
       # dtmvnorm(x_star_n,mean = x_cur_n,sigma = diag(X_kap[,n]),lower = l,upper = u,log = TRUE)
      
      #Ratio
      ratio <- l_star - l_cur +
        
        #Priors
        dmvnorm(x=x_star_n,sigma = GP_cov(d_scaled,lambda = lam_cur[n], ell = ell_cur[n], nugget = 1E-8),
                log = TRUE) - 
        dmvnorm(x=x_cur_n,sigma = GP_cov(d_scaled,lambda = lam_cur[n], ell = ell_cur[n], nugget = 1E-8),
                log = TRUE)  #+
        
        #Proposals
        #Arrggg Botev doesn't include pdfs
    #    adjusted
      
#       if(n %in% c(1,10)){
#         print(paste('xst',x_star_n))
#         print(paste('xcur',x_cur_n))
#         print(paste('l',l))
#         print(paste('u',u))
#       }
      
      if(length(ratio)>1){stop(paste("wtf n=",n))}
      if(runif(1)<exp(ratio)){
        x_acc[n] <- x_acc[n] + 1
        l_cur <- l_star
        X_cur[,n] <- x_star_n
      }
      
      #############
      #Update ell_n
      #############
      z_ell <- rnorm(1)
      ell_st <- ell_cur[n]*exp(delta_ell[n]*z_ell)
      adjust_ratio <- delta_ell[n]*z_ell
      
      #MH Ratio
        
      ratio <- 
        #"Likelihood" - X prior
        dmvnorm(X_cur[,n],sigma = GP_cov(d_scaled,lam_cur[n],ell_st,nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[,n],sigma = GP_cov(d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
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
        dmvnorm(X_cur[,n],sigma = GP_cov(d_scaled,lam_st,ell_cur[n],nugget = 1E-8),log = TRUE) -
        dmvnorm(X_cur[,n],sigma = GP_cov(d_scaled,lam_cur[n],ell_cur[n],nugget = 1E-8), log = TRUE) +
        
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

hope <- hier_gp_mh(iters = 1E3,verbose = TRUE, burnin_prop = 0.3,
                   X_kap = matrix(nrow = I, byrow = FALSE,data = c(
                     rep(1E-7,I),#1
                     rep(1E-5,I),#2
                     rep(1E-6,I), #3
                     rep(1E-4,I), #4
                     rep(1E-7,I),#5
                     rep(1E-6,I),#6
                     rep(1E-5,I),#7
                     rep(1E-4,I),#8
                     rep(1E-3,I))) #9
                   )

hope$x_acc
hope$lam_acc
hope$ell_acc

save(hope_i, file = "run_1k.Rdata")

stop("We're done here")

blah = 1:5
save(blah,file = "didnt_stop.Rdata")
#summaryRprof('Hope_against_hope.out')
X <- hope$X
i <- 4
X_i <- X[,i,]

plot_traces(i,save_pics = FALSE)


plot_dens_i(X_i)
plot(as.numeric(colnames(Y_trunc)),Y_trunc[i,])



gam=1
T_out=seq(0,t_star,length=75)



plot_traces <- function(i,save_pics = FALSE){
  
  for(n in 1:N){
    if(n<ceiling(N/2)){
      var_type = 'Z'
      index = n
    }
    else{
      var_type = 'W'
      index = n-floor(N/2)
    }
    if(save_pics) pdf(paste0(path,var_type,i,suffix,'.pdf'))
    plot(X[,i,n],type = 'l',ylab = bquote(.(var_type)[.(index)]),
         main = bquote('Trace Plot of '~.(var_type)[.(index)]),
         xlab = 'Iteration')
    if(save_pics)dev.off()
  }
}


est_dens_i <- function(x_mat, r){
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




